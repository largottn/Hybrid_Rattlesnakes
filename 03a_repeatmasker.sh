#!/bin/bash
#SBATCH --job-name=repeatmasker
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/repeatmasker.%j.out
#SBATCH --error=logs/repeatmasker.%j.err
#SBATCH --mem=32G

# ============================================================================
# Step 3a - RepeatMasker: identify repetitive regions for variant masking
# ----------------------------------------------------------------------------
# Runs RepeatMasker on the autosomal portion of the reference genome to
# identify repetitive and low-complexity regions, then converts the output
# to a sorted, merged BED file suitable for variant-call masking. Sex
# chromosomes and unplaced scaffolds are excluded up front since downstream
# analyses are restricted to autosomes.
#
# Inputs:   reference FASTA (indexed with samtools faidx), autosome list
# Outputs:  RepeatMasker .out and .gff under ${OUT_DIR}/, plus a merged
#           BED file (${OUT_DIR}/repeat_mask.bed) for downstream masking
#
# Usage:    Set REF, AUTOSOMES, OUT_DIR, and SPECIES below, then submit
#           with sbatch. Run once per reference genome.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and parameters ----
REF="${REF:-./reference/genome.fna}"
AUTOSOMES="${AUTOSOMES:-./reference/autosomes.txt}"
OUT_DIR="${OUT_DIR:-./repeatmasker}"

# RepeatMasker species. Tries the primary first; if unavailable in the
# installed RepBase/Dfam library, falls back to FALLBACK_SPECIES.
SPECIES="${SPECIES:-crotalus viridis}"
FALLBACK_SPECIES="${FALLBACK_SPECIES:-serpentes}"

# ---- Setup ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "${OUT_DIR}" logs

# ---- Sanity checks ----
for tool in samtools RepeatMasker bedtools awk; do
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

# ---- Step 1: Extract autosomes from the reference ----
AUTO_REF="${OUT_DIR}/autosomes_only.fna"

if [[ -f "${AUTO_REF}" ]]; then
    echo "Autosome-only FASTA already exists, skipping extraction."
else
    echo "Extracting autosomes from reference..."
    # shellcheck disable=SC2046
    samtools faidx "${REF}" $(tr '\n' ' ' < "${AUTOSOMES}") > "${AUTO_REF}"
fi

# ---- Step 2: Run RepeatMasker ----
echo "Running RepeatMasker (${SPECIES}) on ${AUTO_REF}..."
echo "Start time: $(date)"

if ! RepeatMasker \
        -species "${SPECIES}" \
        -pa "${THREADS}" \
        -xsmall \
        -gff \
        -dir "${OUT_DIR}" \
        "${AUTO_REF}"; then
    echo "RepeatMasker with '${SPECIES}' failed. Retrying with '${FALLBACK_SPECIES}'..."
    RepeatMasker \
        -species "${FALLBACK_SPECIES}" \
        -pa "${THREADS}" \
        -xsmall \
        -gff \
        -dir "${OUT_DIR}" \
        "${AUTO_REF}"
fi

echo "RepeatMasker finished: $(date)"

# ---- Step 3: Convert RepeatMasker output to merged BED ----
RM_OUT="${OUT_DIR}/$(basename "${AUTO_REF}").out"
MASK_BED="${OUT_DIR}/repeat_mask.bed"

if [[ ! -f "${RM_OUT}" ]]; then
    echo "ERROR: RepeatMasker output file not found: ${RM_OUT}" >&2
    exit 1
fi

echo "Converting RepeatMasker output to merged BED..."

# .out format: 3 header lines, then space-padded columns.
# Columns 5,6,7 = query, start, end (1-based, inclusive).
# Convert to 0-based BED, sort, and merge overlapping intervals.
tail -n +4 "${RM_OUT}" \
    | awk '{print $5"\t"$6-1"\t"$7"\t"$11"\t"$10}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${MASK_BED}"

TOTAL_REGIONS=$(wc -l < "${MASK_BED}")
TOTAL_BP=$(awk '{sum += $3-$2} END {print sum+0}' "${MASK_BED}")

echo ""
echo "============================================"
echo "RepeatMasker complete."
echo "  Merged repeat regions: ${TOTAL_REGIONS}"
echo "  Total masked bp:       ${TOTAL_BP}"
echo "  Mask BED:              ${MASK_BED}"
echo "============================================"
