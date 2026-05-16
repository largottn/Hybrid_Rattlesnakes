#!/bin/bash
#SBATCH --job-name=mie_pipeline
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --output=logs/mie_pipeline.%j.out
#SBATCH --error=logs/mie_pipeline.%j.err

# ============================================================================
# Step 4 - Mendelian Inheritance Error (MIE) detection pipeline
# ----------------------------------------------------------------------------
# Full trio-based MIE detection pipeline using GATK GVCF-mode joint
# genotyping. Designed for a single nuclear family: one dam, one sire,
# and one or more F1 offspring. The pipeline produces, for each trio,
# a list of sites at which the offspring's genotype is inconsistent with
# Mendelian transmission from the parents.
#
# Pipeline overview:
#   1.  Per-sample HaplotypeCaller in GVCF mode (-ERC GVCF)
#   2.  GenomicsDBImport combines per-sample gVCFs
#   3.  GenotypeGVCFs produces a joint multi-sample VCF
#   4.  Restriction to biallelic autosomal SNPs
#   5.  Hard filtering (GATK best practices)
#   6.  Cluster filter (>=3 SNPs in a 10 bp window)
#   7.  Indel proximity filter (SNPs within 50 bp of any indel)
#   8.  Masking of repeat (RepeatMasker) and CNV (CNVnator) regions
#   9.  Per-trio extraction, DP/GQ filtering, and bcftools +mendelian2
#   10. Diagnostic: top MIE genotype patterns per offspring
#   11. Diagnostic: parental genotype composition at callable sites
#
# This is a long-running pipeline (multiple hours per sample for Step 1
# at ~15x coverage on a snake-sized genome). Resume-after-failure is
# partial — Step 1 skips per-sample gVCFs that already exist; later
# steps do not check for existing outputs and will overwrite them.
#
# Inputs:   indexed reference FASTA (samtools faidx + GATK dict),
#           per-sample BAMs in ${BAM_DIR}, repeat + CNV mask BEDs,
#           autosome list file (one chromosome name per line)
# Outputs:  per-sample gVCFs, joint VCF, filtered VCFs at each step,
#           per-trio MIE VCFs and position TSVs
#
# Usage:    Set the paths and SAMPLES below, then submit with sbatch.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and samples ----
WORKDIR="${WORKDIR:-./mie_output}"
BAM_DIR="${BAM_DIR:-./alignments/final_bams}"
REF="${REF:-./reference/genome.fna}"
AUTOSOMES_FILE="${AUTOSOMES_FILE:-./reference/autosomes.list}"
CNV_MASK="${CNV_MASK:-./cnvnator/cnv_mask.bed}"
REPEAT_MASK="${REPEAT_MASK:-./repeatmasker/repeat_mask.bed}"

# Trio structure: one sire, one dam, one or more offspring.
SIRE="${SIRE:-sire_sample}"
DAM="${DAM:-dam_sample}"
OFFSPRING=(child1 child2 child3)

# All five samples participate in joint genotyping.
SAMPLES=("${SIRE}" "${DAM}" "${OFFSPRING[@]}")

# Filtering thresholds
MIN_DP="${MIN_DP:-7}"
MIN_GQ="${MIN_GQ:-30}"
CLUSTER_SIZE="${CLUSTER_SIZE:-3}"
CLUSTER_WINDOW="${CLUSTER_WINDOW:-10}"
INDEL_PROXIMITY_BP="${INDEL_PROXIMITY_BP:-50}"

# ---- Setup ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "${WORKDIR}/gvcfs" "${WORKDIR}/joint" "${WORKDIR}/filtered" logs

# ---- Sanity checks ----
for tool in gatk bcftools vcftools bedtools tabix awk; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done
for f in "${REF}" "${AUTOSOMES_FILE}" "${CNV_MASK}" "${REPEAT_MASK}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done
if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: reference not indexed. Run: samtools faidx ${REF}" >&2
    exit 1
fi
REF_DICT="${REF%.fna}.dict"
if [[ ! -f "${REF_DICT}" && ! -f "${REF%.fa}.dict" && ! -f "${REF%.fasta}.dict" ]]; then
    echo "ERROR: reference dictionary missing. Run: gatk CreateSequenceDictionary -R ${REF}" >&2
    exit 1
fi

# Build a comma-separated chromosome string for bcftools -r
AUTOSOMES_CSV=$(paste -sd, "${AUTOSOMES_FILE}")

# ============================================================================
# STEP 1: HaplotypeCaller in GVCF mode (per sample)
# ============================================================================
echo "=========================================="
echo "STEP 1: HaplotypeCaller in GVCF mode"
echo "=========================================="

for SAMPLE in "${SAMPLES[@]}"; do
    GVCF_OUT="${WORKDIR}/gvcfs/${SAMPLE}.g.vcf.gz"

    if [[ -s "${GVCF_OUT}" ]]; then
        echo "  ${SAMPLE}: gVCF already exists, skipping."
        continue
    fi

    echo "  Calling ${SAMPLE}..."
    gatk HaplotypeCaller \
        -R "${REF}" \
        -I "${BAM_DIR}/${SAMPLE}.bam" \
        -O "${GVCF_OUT}" \
        -ERC GVCF \
        -L "${AUTOSOMES_FILE}" \
        --native-pair-hmm-threads "${THREADS}"

    if [[ ! -s "${GVCF_OUT}" ]]; then
        echo "ERROR: HaplotypeCaller produced empty output for ${SAMPLE}" >&2
        exit 1
    fi
    echo "  ${SAMPLE}: done."
done

# Verify all gVCFs exist before proceeding
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ ! -s "${WORKDIR}/gvcfs/${SAMPLE}.g.vcf.gz" ]]; then
        echo "ERROR: missing gVCF for ${SAMPLE}" >&2
        exit 1
    fi
done
echo "All gVCFs generated and verified."
echo ""

# ============================================================================
# STEP 2: GenomicsDBImport
# ============================================================================
echo "=========================================="
echo "STEP 2: GenomicsDBImport"
echo "=========================================="

GENOMICSDB="${WORKDIR}/joint/genomicsdb"
SAMPLE_MAP="${WORKDIR}/joint/sample_map.txt"

: > "${SAMPLE_MAP}"
for SAMPLE in "${SAMPLES[@]}"; do
    printf '%s\t%s\n' "${SAMPLE}" "${WORKDIR}/gvcfs/${SAMPLE}.g.vcf.gz" >> "${SAMPLE_MAP}"
done

# GenomicsDB workspace cannot already exist; remove if present.
if [[ -d "${GENOMICSDB}" ]]; then
    echo "  Removing existing GenomicsDB workspace..."
    rm -rf "${GENOMICSDB}"
fi

gatk GenomicsDBImport \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${GENOMICSDB}" \
    -L "${AUTOSOMES_FILE}" \
    --reader-threads "${THREADS}"

echo "GenomicsDBImport done."
echo ""

# ============================================================================
# STEP 3: GenotypeGVCFs
# ============================================================================
echo "=========================================="
echo "STEP 3: GenotypeGVCFs"
echo "=========================================="

JOINT_VCF="${WORKDIR}/joint/all_samples.joint.vcf.gz"

gatk GenotypeGVCFs \
    -R "${REF}" \
    -V "gendb://${GENOMICSDB}" \
    -O "${JOINT_VCF}"

bcftools index "${JOINT_VCF}"

echo "  Joint genotyped sites: $(bcftools view -H "${JOINT_VCF}" | wc -l)"
echo ""

# ============================================================================
# STEP 4: Biallelic autosomal SNPs
# ============================================================================
echo "=========================================="
echo "STEP 4: Biallelic autosomal SNPs"
echo "=========================================="

BIALLELIC_VCF="${WORKDIR}/filtered/biallelic_snps.vcf.gz"

bcftools view \
    -v snps -m2 -M2 \
    -r "${AUTOSOMES_CSV}" \
    "${JOINT_VCF}" \
    -Oz -o "${BIALLELIC_VCF}"
bcftools index "${BIALLELIC_VCF}"

echo "  Biallelic autosomal SNPs: $(bcftools view -H "${BIALLELIC_VCF}" | wc -l)"
echo ""

# ============================================================================
# STEP 5: Hard filtering (GATK best practices)
# ============================================================================
echo "=========================================="
echo "STEP 5: Hard filtering"
echo "=========================================="

HARDFILTERED_VCF="${WORKDIR}/filtered/hard_filtered.vcf.gz"

bcftools filter \
    -e 'QUAL<30 || QD<2.0 || FS>60.0 || MQ<40.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0' \
    "${BIALLELIC_VCF}" \
    -Oz -o "${HARDFILTERED_VCF}"
bcftools index "${HARDFILTERED_VCF}"

echo "  After hard filtering: $(bcftools view -H "${HARDFILTERED_VCF}" | wc -l)"
echo ""

# ============================================================================
# STEP 6: Cluster filter (>=3 SNPs in a 10 bp window)
# ============================================================================
echo "=========================================="
echo "STEP 6: Cluster filter"
echo "=========================================="

CLUSTER_TAGGED="${WORKDIR}/filtered/cluster_tagged.vcf.gz"
CLUSTER_VCF="${WORKDIR}/filtered/no_clusters.vcf.gz"

gatk VariantFiltration \
    -R "${REF}" \
    -V "${HARDFILTERED_VCF}" \
    --cluster-size "${CLUSTER_SIZE}" \
    --cluster-window-size "${CLUSTER_WINDOW}" \
    -O "${CLUSTER_TAGGED}"

# Retain sites passing the cluster filter (PASS or unannotated)
bcftools view -f .,PASS \
    "${CLUSTER_TAGGED}" \
    -Oz -o "${CLUSTER_VCF}"
bcftools index "${CLUSTER_VCF}"

echo "  After cluster filter: $(bcftools view -H "${CLUSTER_VCF}" | wc -l)"
echo ""

# ============================================================================
# STEP 7: Indel proximity filter (remove SNPs within INDEL_PROXIMITY_BP of any indel)
# ============================================================================
echo "=========================================="
echo "STEP 7: Indel proximity filter (${INDEL_PROXIMITY_BP} bp)"
echo "=========================================="

INDEL_BED="${WORKDIR}/filtered/indel_proximity.bed"
PROXIMITY_VCF="${WORKDIR}/filtered/no_indel_prox.vcf.gz"

bcftools view -v indels "${JOINT_VCF}" \
    | bcftools query -f '%CHROM\t%POS\n' \
    | awk -v p="${INDEL_PROXIMITY_BP}" \
        '{OFS="\t"; s=$2-p; if(s<0) s=0; print $1, s, $2+p}' \
    > "${INDEL_BED}"

bcftools view \
    -T "^${INDEL_BED}" \
    "${CLUSTER_VCF}" \
    -Oz -o "${PROXIMITY_VCF}"
bcftools index "${PROXIMITY_VCF}"

echo "  After indel proximity filter: $(bcftools view -H "${PROXIMITY_VCF}" | wc -l)"
echo ""

# ============================================================================
# STEP 8: CNV + Repeat masking
# ============================================================================
echo "=========================================="
echo "STEP 8: CNV + Repeat masking"
echo "=========================================="

COMBINED_MASK="${WORKDIR}/filtered/combined_mask.bed"
MASKED_VCF="${WORKDIR}/filtered/masked.vcf.gz"

cat "${CNV_MASK}" "${REPEAT_MASK}" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${COMBINED_MASK}"

bcftools view \
    -T "^${COMBINED_MASK}" \
    "${PROXIMITY_VCF}" \
    -Oz -o "${MASKED_VCF}"
bcftools index "${MASKED_VCF}"

echo "  After CNV + repeat masking: $(bcftools view -H "${MASKED_VCF}" | wc -l)"
echo ""

# ============================================================================
# STEP 9: Per-trio MIE analysis
# ============================================================================
echo "=========================================="
echo "STEP 9: Per-trio MIE analysis (DP >= ${MIN_DP}, GQ >= ${MIN_GQ})"
echo "=========================================="

for KID in "${OFFSPRING[@]}"; do
    echo "--- Processing trio: ${KID} ---"

    TRIO_DIR="${WORKDIR}/${KID}_trio"
    mkdir -p "${TRIO_DIR}"

    TRIO_VCF="${TRIO_DIR}/trio.vcf.gz"
    NOMISSING_VCF="${TRIO_DIR}/trio.nomissing.vcf.gz"
    FILTERED_VCF="${TRIO_DIR}/trio.filtered.vcf.gz"
    PED_FILE="${TRIO_DIR}/trio.ped"
    MIE_POS="${TRIO_DIR}/${KID}.mie.positions.tsv"
    MIE_VCF="${TRIO_DIR}/${KID}.mie.vcf.gz"

    # Extract trio samples
    bcftools view \
        -s "${DAM},${SIRE},${KID}" \
        "${MASKED_VCF}" \
        -Oz -o "${TRIO_VCF}"
    bcftools index "${TRIO_VCF}"

    TRIO_COUNT=$(bcftools view -H "${TRIO_VCF}" | wc -l)
    echo "  Trio sites: ${TRIO_COUNT}"

    # Remove sites with any missing genotype
    bcftools view -i 'F_MISSING=0' \
        "${TRIO_VCF}" \
        -Oz -o "${NOMISSING_VCF}"
    bcftools index "${NOMISSING_VCF}"

    NOMISSING_COUNT=$(bcftools view -H "${NOMISSING_VCF}" | wc -l)
    echo "  Sites with no missing data: ${NOMISSING_COUNT}"

    # DP and GQ thresholds
    bcftools view \
        -i "MIN(FMT/DP)>=${MIN_DP} && MIN(FMT/GQ)>=${MIN_GQ}" \
        "${NOMISSING_VCF}" \
        -Oz -o "${FILTERED_VCF}"
    bcftools index "${FILTERED_VCF}"

    FILTERED_COUNT=$(bcftools view -H "${FILTERED_VCF}" | wc -l)
    echo "  Callable sites after DP/GQ filter: ${FILTERED_COUNT}"

    # PED file: FAM ID, sample ID, father, mother, sex (1=M/2=F), phenotype (0)
    cat > "${PED_FILE}" << EOF
FAM1	${DAM}	0	0	2	0
FAM1	${SIRE}	0	0	1	0
FAM1	${KID}	${SIRE}	${DAM}	0	0
EOF

    # MIE counts (mode c = summary counts)
    bcftools +mendelian2 \
        "${FILTERED_VCF}" \
        -P "${PED_FILE}" \
        -m c 2>&1 | tee "${TRIO_DIR}/${KID}.mendelian2.log"

    # MIE positions (mode e = error sites)
    bcftools +mendelian2 \
        "${FILTERED_VCF}" \
        -P "${PED_FILE}" \
        -m e \
        | bcftools query -f '%CHROM\t%POS\n' \
        > "${MIE_POS}"

    # MIE VCF (mode e, VCF output)
    bcftools +mendelian2 \
        "${FILTERED_VCF}" \
        -P "${PED_FILE}" \
        -m e \
        -Oz -o "${MIE_VCF}"
    bcftools index "${MIE_VCF}"

    MIE_COUNT=$(wc -l < "${MIE_POS}")

    # Independent validation with vcftools --mendel
    vcftools \
        --gzvcf "${FILTERED_VCF}" \
        --mendel "${PED_FILE}" \
        --out "${TRIO_DIR}/${KID}.vcftools_mendel" 2>&1

    echo ""
    echo "  --- SUMMARY for ${KID} ---"
    echo "  Trio sites:                 ${TRIO_COUNT}"
    echo "  No missing data:            ${NOMISSING_COUNT}"
    echo "  After DP/GQ filter:         ${FILTERED_COUNT}"
    echo "  MIE count (mendelian2):     ${MIE_COUNT}"
    if [[ ${FILTERED_COUNT} -gt 0 ]]; then
        RATE=$(awk "BEGIN {printf \"%.4f\", 100*${MIE_COUNT}/${FILTERED_COUNT}}")
        echo "  MIE rate:                   ${RATE}%"
    fi
    echo ""
done

# ============================================================================
# STEP 10 (diagnostic): top MIE genotype patterns per offspring
# ============================================================================
echo "=========================================="
echo "STEP 10 (diagnostic): MIE genotype patterns"
echo "=========================================="

for KID in "${OFFSPRING[@]}"; do
    TRIO_DIR="${WORKDIR}/${KID}_trio"
    echo "--- ${KID} ---"
    TOTAL=$(wc -l < "${TRIO_DIR}/${KID}.mie.positions.tsv")
    echo "  Total MIEs: ${TOTAL}"
    if [[ ${TOTAL} -gt 0 ]]; then
        echo "  Top genotype patterns (dam sire offspring):"
        bcftools query -f '[%GT\t]\n' "${TRIO_DIR}/${KID}.mie.vcf.gz" \
            | awk '{print $1, $2, $3}' | sort | uniq -c | sort -rn | head -10
    fi
    echo ""
done

# ============================================================================
# STEP 11 (diagnostic): callable site parental genotype composition
# ============================================================================
echo "=========================================="
echo "STEP 11 (diagnostic): callable site composition"
echo "=========================================="

for KID in "${OFFSPRING[@]}"; do
    TRIO_DIR="${WORKDIR}/${KID}_trio"
    echo "--- ${KID} ---"
    echo "  Parental genotype combinations (dam x sire) at callable sites (top 10):"
    bcftools query -f '[%GT\t]\n' "${TRIO_DIR}/trio.filtered.vcf.gz" \
        | awk '{print $1, $2}' | sort | uniq -c | sort -rn | head -10
    echo ""
done

echo "=========================================="
echo "MIE PIPELINE COMPLETE"
echo "=========================================="
