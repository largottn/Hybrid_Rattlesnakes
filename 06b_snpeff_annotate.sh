#!/bin/bash
#SBATCH --job-name=snpeff_annotate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-2
#SBATCH --output=logs/snpeff_annotate.%a.%j.out
#SBATCH --error=logs/snpeff_annotate.%a.%j.err
#SBATCH --mem=8G

# ============================================================================
# Step 6b - SnpEff annotation of per-trio de novo mutations
# ----------------------------------------------------------------------------
# For each offspring, converts the final DNM TSV (from Step 5c) into a
# minimal VCF, annotates it with the custom SnpEff database built in
# Step 6a, and parses the annotated VCF into a human-readable summary.
#
# This script assumes Step 6a (database build) has already completed.
# Submit Step 6a first, wait for it to finish, then submit this script.
#
# The per-offspring SLURM array means one task runs per trio in parallel,
# all reading from the same (read-only) database.
#
# Inputs:   per-trio DNM TSVs from Step 5c in ${DNM_RESULTS_DIR}/
#           SnpEff database + config from Step 6a in ${SNPEFF_DIR}/
#           sibling Python scripts: dnm_tsv_to_vcf.py, parse_snpeff_results.py
# Outputs:  per-trio annotated VCFs and parsed TSVs in ${OUT_DIR}/
#
# Usage:    Set OFFSPRING, paths, and GENOME_NAME below.
#           Adjust #SBATCH --array=0-N to match (number of offspring - 1).
#           Submit with sbatch after Step 6a has completed.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and samples ----
SNPEFF_DIR="${SNPEFF_DIR:-./snpeff}"
GENOME_NAME="${GENOME_NAME:-custom_genome}"

DNM_RESULTS_DIR="${DNM_RESULTS_DIR:-./dnm_output/results}"
OUT_DIR="${OUT_DIR:-./dnm_output/snpeff_annotated}"

OFFSPRING=(child1 child2 child3)

# Path to the sibling Python helpers
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TSV_TO_VCF="${TSV_TO_VCF:-${SCRIPT_DIR}/dnm_tsv_to_vcf.py}"
PARSE_RESULTS="${PARSE_RESULTS:-${SCRIPT_DIR}/parse_snpeff_results.py}"

# ---- Setup ----
mkdir -p "${OUT_DIR}" logs

# ---- Sanity checks ----
for tool in snpEff python3; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done

CONFIG_FILE="${SNPEFF_DIR}/snpEff.config"
DB_FILE="${SNPEFF_DIR}/data/${GENOME_NAME}/snpEffectPredictor.bin"
for f in "${CONFIG_FILE}" "${DB_FILE}" "${TSV_TO_VCF}" "${PARSE_RESULTS}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required file not found: ${f}" >&2
        echo "  (Did Step 6a complete?)" >&2
        exit 1
    fi
done

CHILD="${OFFSPRING[${SLURM_ARRAY_TASK_ID:-0}]}"
INPUT_TSV="${DNM_RESULTS_DIR}/trio_${CHILD}_final_denovo.tsv"

if [[ ! -f "${INPUT_TSV}" ]]; then
    echo "ERROR: DNM TSV not found for ${CHILD}: ${INPUT_TSV}" >&2
    exit 1
fi

echo "============================================"
echo "  Offspring: ${CHILD}"
echo "  Input TSV: ${INPUT_TSV}"
echo "  Output:    ${OUT_DIR}"
echo "  Start:     $(date)"
echo "============================================"

# ---- Short-circuit when there are no variants to annotate ----
N_VARIANTS=$(($(wc -l < "${INPUT_TSV}") - 1))   # subtract header
if [[ "${N_VARIANTS}" -le 0 ]]; then
    echo "No DNMs to annotate for ${CHILD}. Exiting cleanly."
    exit 0
fi
echo "Annotating ${N_VARIANTS} DNM(s)."

# ---- Step 1: convert DNM TSV to minimal VCF for SnpEff input ----
VCF_FOR_SNPEFF="${OUT_DIR}/trio_${CHILD}_denovo_for_snpeff.vcf"
python3 "${TSV_TO_VCF}" \
    --input  "${INPUT_TSV}" \
    --output "${VCF_FOR_SNPEFF}"

# ---- Step 2: run SnpEff annotation ----
ANNOTATED_VCF="${OUT_DIR}/trio_${CHILD}_denovo_annotated.vcf"
STATS_HTML="${OUT_DIR}/trio_${CHILD}_snpeff_stats.html"

snpEff ann \
    -c "${CONFIG_FILE}" \
    -stats "${STATS_HTML}" \
    -v \
    "${GENOME_NAME}" \
    "${VCF_FOR_SNPEFF}" \
    > "${ANNOTATED_VCF}"

# ---- Step 3: parse annotated VCF into human-readable summary ----
SUMMARY_TSV="${OUT_DIR}/trio_${CHILD}_denovo_annotation_summary.tsv"
python3 "${PARSE_RESULTS}" \
    --input  "${ANNOTATED_VCF}" \
    --output "${SUMMARY_TSV}"

echo ""
echo "============================================"
echo "  ${CHILD} annotation complete."
echo "  Annotated VCF: ${ANNOTATED_VCF}"
echo "  Summary TSV:   ${SUMMARY_TSV}"
echo "  Stats HTML:    ${STATS_HTML}"
echo "  End:           $(date)"
echo "============================================"
