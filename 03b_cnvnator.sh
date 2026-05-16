#!/bin/bash
#SBATCH --job-name=cnvnator
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/cnvnator.%j.out
#SBATCH --error=logs/cnvnator.%j.err
#SBATCH --mem=32G

# ============================================================================
# Step 3b - CNVnator: detect copy number variants for variant masking
# ----------------------------------------------------------------------------
# Runs CNVnator on each sample's deduplicated BAM to detect copy number
# variants (CNVs) across autosomal chromosomes. Filters calls by statistical
# significance and minimum size, producing a per-sample BED of CNV regions
# suitable for masking in downstream variant analyses.
#
# Two preparatory artifacts are produced the first time the script runs
# and reused on subsequent runs:
#   * Per-chromosome FASTAs in ${CHROM_DIR}/  (required by `cnvnator -his`)
#   * A `.split_done` lockfile in ${CHROM_DIR}/
#
# Inputs:   reference FASTA (indexed), autosome list, per-sample BAMs
# Outputs:  per-sample CNVnator ROOT file, raw call text, and filtered BED
#
# Filter:   e-value < 0.05 (both natural and Gaussian) AND size >= 1 kb
#
# Bin size: tuned to coverage. The default (300 bp) is appropriate for
#           ~15x coverage; ~10x → 500 bp; ~30x+ → 100 bp (see CNVnator docs).
#
# Note:     This script produces per-sample CNV BEDs only. To use them
#           as a single mask across all samples (recommended for trio
#           analyses), merge them with `bedtools merge` in a downstream step.
#
# Usage:    Set REF, AUTOSOMES, BAM_DIR, OUT_DIR, CHROM_DIR, SAMPLES,
#           and BIN_SIZE below, then submit with sbatch.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and parameters ----
REF="${REF:-./reference/genome.fna}"
AUTOSOMES="${AUTOSOMES:-./reference/autosomes.txt}"
BAM_DIR="${BAM_DIR:-./alignments/final_bams}"
OUT_DIR="${OUT_DIR:-./cnvnator}"
CHROM_DIR="${CHROM_DIR:-./reference/chroms}"

BIN_SIZE="${BIN_SIZE:-300}"

# Sample IDs to process. Expects ${BAM_DIR}/${SAMPLE}.bam for each.
SAMPLES=(sample1 sample2 sample3 sample4)

# ---- Setup ----
mkdir -p "${OUT_DIR}" "${CHROM_DIR}" logs

# ---- Sanity checks ----
for tool in cnvnator samtools awk; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done
for f in "${REF}" "${AUTOSOMES}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done

# ---- Step 0: split reference into per-chromosome FASTAs (once) ----
LOCK_FILE="${CHROM_DIR}/.split_done"
if [[ ! -f "${LOCK_FILE}" ]]; then
    echo "Splitting reference into per-chromosome FASTAs (one-time setup)..."
    while IFS= read -r chrom; do
        [[ -z "${chrom}" ]] && continue
        samtools faidx "${REF}" "${chrom}" > "${CHROM_DIR}/${chrom}.fa"
    done < "${AUTOSOMES}"
    touch "${LOCK_FILE}"
    echo "  Per-chromosome FASTAs written to ${CHROM_DIR}/"
else
    echo "Per-chromosome FASTAs already exist, skipping split."
fi

# Build a space-separated chromosome list for cnvnator's -chrom flag.
# shellcheck disable=SC2002
CHROM_STRING=$(tr '\n' ' ' < "${AUTOSOMES}")

# ---- Process each sample ----
for SAMPLE in "${SAMPLES[@]}"; do
    echo "========================================"
    echo "Processing sample: ${SAMPLE}"
    echo "========================================"

    BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"
    ROOT_FILE="${OUT_DIR}/${SAMPLE}.root"

    if [[ ! -f "${BAM_FILE}" ]]; then
        echo "ERROR: BAM file not found: ${BAM_FILE}. Skipping." >&2
        continue
    fi

    # CNVnator's ROOT file is rebuilt from scratch each run to avoid
    # corruption from prior failed jobs.
    rm -f "${ROOT_FILE}"

    echo "Step 1: extracting reads..."
    # shellcheck disable=SC2086
    cnvnator -root "${ROOT_FILE}" -tree "${BAM_FILE}" -chrom ${CHROM_STRING}

    echo "Step 2: generating histograms (bin size = ${BIN_SIZE})..."
    cnvnator -root "${ROOT_FILE}" -his "${BIN_SIZE}" -d "${CHROM_DIR}"

    echo "Step 3: calculating statistics..."
    cnvnator -root "${ROOT_FILE}" -stat "${BIN_SIZE}"

    echo "Step 4: partitioning..."
    cnvnator -root "${ROOT_FILE}" -partition "${BIN_SIZE}"

    echo "Step 5: calling CNVs..."
    cnvnator -root "${ROOT_FILE}" -call "${BIN_SIZE}" \
        > "${OUT_DIR}/${SAMPLE}.cnvs.raw.txt"

    echo "Step 6: filtering CNVs (e-value < 0.05, size >= 1 kb)..."
    awk -F'\t' '{
        split($2, coords, "[:-]")
        chrom = coords[1]
        start = coords[2]
        end   = coords[3]
        size  = $3
        eval1 = $5
        eval2 = $6
        if (eval1 < 0.05 && eval2 < 0.05 && size >= 1000) {
            print chrom "\t" start "\t" end "\t" $1 "\t" size
        }
    }' "${OUT_DIR}/${SAMPLE}.cnvs.raw.txt" \
        > "${OUT_DIR}/${SAMPLE}.cnvs.filtered.bed"

    RAW=$(wc -l < "${OUT_DIR}/${SAMPLE}.cnvs.raw.txt")
    FILT=$(wc -l < "${OUT_DIR}/${SAMPLE}.cnvs.filtered.bed")
    echo "  Raw CNV calls:      ${RAW}"
    echo "  Filtered CNV calls: ${FILT}"
    echo "Done with ${SAMPLE}."
done

echo ""
echo "============================================"
echo "All samples processed."
echo "  Per-sample filtered CNV BEDs in: ${OUT_DIR}/"
echo "  To produce a single cross-sample mask:"
echo "    cat ${OUT_DIR}/*.cnvs.filtered.bed \\"
echo "      | sort -k1,1 -k2,2n \\"
echo "      | bedtools merge -i - > ${OUT_DIR}/cnv_mask.bed"
echo "============================================"
