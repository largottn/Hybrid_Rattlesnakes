#!/bin/bash
#SBATCH --job-name=dnm_refine
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-2
#SBATCH --output=logs/dnm_refine.%a.%j.out
#SBATCH --error=logs/dnm_refine.%a.%j.err
#SBATCH --mem=32G

# ============================================================================
# Step 5b - DNM pipeline: hard filtering and Genotype Refinement workflow
# ----------------------------------------------------------------------------
# Takes the raw trio VCFs from Step 5a (multi-sample HaplotypeCaller) and:
#   1. Hard-filters SNVs and indels with separate thresholds
#   2. Removes variant clusters (>= CLUSTER_SIZE within CLUSTER_WINDOW bp)
#   3. Runs GATK CalculateGenotypePosteriors with a trio PED file
#      (--skip-population-priors, since we have no population-level data)
#   4. Marks low-GQ genotypes via VariantFiltration
#   5. Annotates with PossibleDeNovo via VariantAnnotator
#
# Output is a per-trio annotated VCF carrying the hiConfDeNovo / loConfDeNovo
# INFO fields, ready for Step 5c (read-level filtering).
#
# Methodology adapted from Zhang et al. 2020 (Nat Commun) and GATK Best
# Practices for de novo mutation detection.
#
# Inputs:   per-trio raw VCFs from Step 5a in ${WORK_DIR}/raw_calls/
# Outputs:  per-trio annotated VCFs in ${WORK_DIR}/refined/
#
# Usage:    Set REF, WORK_DIR, SIRE, DAM, and OFFSPRING below.
#           Adjust #SBATCH --array=0-N to match (number of offspring - 1).
#           Set CHILD_SEX_MAP entries for each offspring (1=male, 2=female,
#           0=unknown). Submit with sbatch.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and samples ----
REF="${REF:-./reference/genome.fna}"
WORK_DIR="${WORK_DIR:-./dnm_output}"

SIRE="${SIRE:-sire_sample}"
DAM="${DAM:-dam_sample}"
OFFSPRING=(child1 child2 child3)

# Sex of each offspring for the PED file: 1=male, 2=female, 0=unknown.
# Required because parental sex matters for downstream sex-chromosome
# interpretation (autosomal DNMs are unaffected).
declare -A CHILD_SEX_MAP=(
    [child1]=0
    [child2]=0
    [child3]=0
)

# ---- Hard-filter thresholds (Zhang et al. 2020) ----
# SNVs: QD<2.0, FS>60, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8.0, QUAL<100
# Indels: QD<2.0, FS>200, ReadPosRankSum<-20, QUAL<100
# Cluster filter: >= CLUSTER_SIZE variants within CLUSTER_WINDOW bp
CLUSTER_SIZE="${CLUSTER_SIZE:-5}"
CLUSTER_WINDOW="${CLUSTER_WINDOW:-100}"
GQ_MIN="${GQ_MIN:-30}"

# ---- Setup ----
mkdir -p "${WORK_DIR}/raw_calls" "${WORK_DIR}/refined" "${WORK_DIR}/ped" logs

# ---- Sanity checks ----
for tool in gatk bcftools; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done
if [[ ! -f "${REF}" || ! -f "${REF}.fai" ]]; then
    echo "ERROR: reference not found or not indexed: ${REF}" >&2
    exit 1
fi

CHILD="${OFFSPRING[${SLURM_ARRAY_TASK_ID:-0}]}"
CHILD_SEX="${CHILD_SEX_MAP[${CHILD}]:-0}"

RAW_VCF="${WORK_DIR}/raw_calls/trio_${CHILD}_raw.vcf.gz"
if [[ ! -s "${RAW_VCF}" ]]; then
    echo "ERROR: raw trio VCF not found: ${RAW_VCF}" >&2
    echo "  Run Step 5a (DNM trio HaplotypeCaller) first." >&2
    exit 1
fi

# ---- PED file for this trio ----
# PED format: FamilyID IndividualID PaternalID MaternalID Sex Phenotype
# Sex: 1=male, 2=female, 0=unknown
# Phenotype: 0=unknown, 1=unaffected, 2=affected
PED_FILE="${WORK_DIR}/ped/trio_${CHILD}.ped"
cat > "${PED_FILE}" << EOF
trio_${CHILD}	${SIRE}	0	0	1	1
trio_${CHILD}	${DAM}	0	0	2	1
trio_${CHILD}	${CHILD}	${SIRE}	${DAM}	${CHILD_SEX}	1
EOF

echo "============================================"
echo "  Trio:     ${SIRE} + ${DAM} -> ${CHILD}"
echo "  PED file: ${PED_FILE}"
echo "  Start:    $(date)"
echo "============================================"

# ---- Convenience: short alias to avoid repeating -Xmx flag ----
GATK="gatk --java-options -Xmx24g"

# Step 1: select and hard-filter SNVs
echo "[$(date)] Selecting and hard-filtering SNVs..."
${GATK} SelectVariants \
    -R "${REF}" -V "${RAW_VCF}" \
    --select-type-to-include SNP \
    -O "${WORK_DIR}/refined/trio_${CHILD}_raw_snps.vcf.gz"

${GATK} VariantFiltration \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_raw_snps.vcf.gz" \
    --filter-expression "QD < 2.0"             --filter-name "QD2" \
    --filter-expression "FS > 60.0"            --filter-name "FS60" \
    --filter-expression "MQ < 40.0"            --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"    --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "QUAL < 100.0"         --filter-name "QUAL100" \
    -O "${WORK_DIR}/refined/trio_${CHILD}_filtered_snps.vcf.gz"

# Step 2: select and hard-filter indels (different thresholds)
echo "[$(date)] Selecting and hard-filtering indels..."
${GATK} SelectVariants \
    -R "${REF}" -V "${RAW_VCF}" \
    --select-type-to-include INDEL \
    -O "${WORK_DIR}/refined/trio_${CHILD}_raw_indels.vcf.gz"

${GATK} VariantFiltration \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_raw_indels.vcf.gz" \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 200.0"            --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --filter-expression "QUAL < 100.0"          --filter-name "QUAL100" \
    -O "${WORK_DIR}/refined/trio_${CHILD}_filtered_indels.vcf.gz"

# Step 3: merge SNVs + indels, keep only PASS
echo "[$(date)] Merging filtered SNVs and indels..."
${GATK} MergeVcfs \
    -I "${WORK_DIR}/refined/trio_${CHILD}_filtered_snps.vcf.gz" \
    -I "${WORK_DIR}/refined/trio_${CHILD}_filtered_indels.vcf.gz" \
    -O "${WORK_DIR}/refined/trio_${CHILD}_filtered_merged.vcf.gz"

${GATK} SelectVariants \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_filtered_merged.vcf.gz" \
    --exclude-filtered \
    -O "${WORK_DIR}/refined/trio_${CHILD}_pass.vcf.gz"

# Step 4: cluster filter
echo "[$(date)] Filtering variant clusters (${CLUSTER_SIZE} in ${CLUSTER_WINDOW} bp)..."
${GATK} VariantFiltration \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_pass.vcf.gz" \
    --cluster-size "${CLUSTER_SIZE}" \
    --cluster-window-size "${CLUSTER_WINDOW}" \
    -O "${WORK_DIR}/refined/trio_${CHILD}_cluster_flagged.vcf.gz"

${GATK} SelectVariants \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_cluster_flagged.vcf.gz" \
    --exclude-filtered \
    -O "${WORK_DIR}/refined/trio_${CHILD}_hardfiltered.vcf.gz"

# Step 5: CalculateGenotypePosteriors
echo "[$(date)] Calculating genotype posteriors..."
${GATK} CalculateGenotypePosteriors \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_hardfiltered.vcf.gz" \
    --pedigree "${PED_FILE}" \
    --skip-population-priors \
    -O "${WORK_DIR}/refined/trio_${CHILD}_posteriors.vcf.gz"

# Step 6: mark low-GQ genotypes
echo "[$(date)] Marking low-GQ genotypes (GQ < ${GQ_MIN})..."
${GATK} VariantFiltration \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_posteriors.vcf.gz" \
    --genotype-filter-expression "GQ < ${GQ_MIN}" \
    --genotype-filter-name "lowGQ" \
    -O "${WORK_DIR}/refined/trio_${CHILD}_gq_filtered.vcf.gz"

# Step 7: annotate PossibleDeNovo
echo "[$(date)] Annotating PossibleDeNovo..."
${GATK} VariantAnnotator \
    -R "${REF}" \
    -V "${WORK_DIR}/refined/trio_${CHILD}_gq_filtered.vcf.gz" \
    --pedigree "${PED_FILE}" \
    -A PossibleDeNovo \
    -O "${WORK_DIR}/refined/trio_${CHILD}_denovo_annotated.vcf.gz"

# ---- Summary ----
echo ""
echo "============================================"
echo "Pipeline summary for trio: ${CHILD}"
echo "============================================"

HICONF=$(bcftools query \
    -f '%INFO/hiConfDeNovo\n' \
    "${WORK_DIR}/refined/trio_${CHILD}_denovo_annotated.vcf.gz" \
    | grep -cv '^\.$' || true)

LOCONF=$(bcftools query \
    -f '%INFO/loConfDeNovo\n' \
    "${WORK_DIR}/refined/trio_${CHILD}_denovo_annotated.vcf.gz" \
    | grep -cv '^\.$' || true)

echo "  High-confidence de novo candidates: ${HICONF}"
echo "  Low-confidence de novo candidates:  ${LOCONF}"
echo "  Annotated VCF: ${WORK_DIR}/refined/trio_${CHILD}_denovo_annotated.vcf.gz"
echo ""
echo "  End: $(date)"
echo "============================================"
