#!/bin/bash
#SBATCH --job-name=dnm_readlevel
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-2
#SBATCH --output=logs/dnm_readlevel.%a.%j.out
#SBATCH --error=logs/dnm_readlevel.%a.%j.err
#SBATCH --mem=8G

# ============================================================================
# Step 5c - DNM pipeline: read-level filtering of de novo candidates
# ----------------------------------------------------------------------------
# Takes the annotated trio VCF from Step 5b and applies the read-level
# filters from Zhang et al. 2020 to the hiConfDeNovo candidates:
#
#   1. Minimum depth (DP) >= MIN_DP in all three trio members
#   2. Parental allele filter:
#        - Hom-ref parents must have zero alt-allele reads
#        - Hom-alt parents must have zero ref-allele reads
#        - Heterozygous parents disqualify the site
#   3. Offspring allele balance in [AB_MIN, AB_MAX]
#   4. Restriction to autosomes
#
# The actual filtering logic lives in filter_dnm_candidates.py (called
# from this script). Separating Python from bash makes the logic
# auditable and individually testable.
#
# Note: Zhang et al. 2020 used DP >= 20 with their 30-40x coverage. For
# lower coverage data, adjust MIN_DP. Each halving of DP roughly doubles
# false-positive risk, so changes here are worth sensitivity-testing.
#
# Inputs:   annotated VCF from Step 5b in ${WORK_DIR}/refined/
#           autosome list (one chromosome name per line)
#           sibling Python script: filter_dnm_candidates.py
# Outputs:  per-trio TSV of final DNMs in ${WORK_DIR}/results/
#
# Usage:    Set REF-related paths, SIRE, DAM, OFFSPRING, and thresholds
#           below. Submit with sbatch.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and samples ----
WORK_DIR="${WORK_DIR:-./dnm_output}"
AUTOSOMES_FILE="${AUTOSOMES_FILE:-./reference/autosomes.list}"

SIRE="${SIRE:-sire_sample}"
DAM="${DAM:-dam_sample}"
OFFSPRING=(child1 child2 child3)

# Read-level filter thresholds
MIN_DP="${MIN_DP:-7}"      # Zhang et al. used 20 with ~30-40x; lower for lower coverage
AB_MIN="${AB_MIN:-0.3}"
AB_MAX="${AB_MAX:-0.7}"

# Path to the Python filter (expected as a sibling of this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER_PY="${FILTER_PY:-${SCRIPT_DIR}/filter_dnm_candidates.py}"

# ---- Setup ----
mkdir -p "${WORK_DIR}/results" "${WORK_DIR}/extracted" logs

# ---- Sanity checks ----
for tool in bcftools python3; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done
for f in "${AUTOSOMES_FILE}" "${FILTER_PY}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done

CHILD="${OFFSPRING[${SLURM_ARRAY_TASK_ID:-0}]}"
DENOVO_VCF="${WORK_DIR}/refined/trio_${CHILD}_denovo_annotated.vcf.gz"

if [[ ! -s "${DENOVO_VCF}" ]]; then
    echo "ERROR: annotated VCF not found: ${DENOVO_VCF}" >&2
    echo "  Run Step 5b (Genotype Refinement) first." >&2
    exit 1
fi

echo "============================================"
echo "  Trio:       ${SIRE} + ${DAM} -> ${CHILD}"
echo "  Input VCF:  ${DENOVO_VCF}"
echo "  DP min:     ${MIN_DP}"
echo "  AB range:   [${AB_MIN}, ${AB_MAX}]"
echo "  Start:      $(date)"
echo "============================================"

# Show sample order in the VCF (used by the Python filter to disambiguate columns)
echo ""
echo "Sample order in VCF:"
bcftools query -l "${DENOVO_VCF}"
echo ""

# ---- Extract hiConfDeNovo candidates with per-sample fields ----
AUTOSOMES_CSV=$(paste -sd, "${AUTOSOMES_FILE}")
EXTRACTED="${WORK_DIR}/extracted/trio_${CHILD}_denovo_raw_extract.tsv"

bcftools query \
    -r "${AUTOSOMES_CSV}" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/hiConfDeNovo\t[%SAMPLE\t%GT\t%DP\t%AD\t]\n' \
    -i 'INFO/hiConfDeNovo!="."' \
    "${DENOVO_VCF}" \
    > "${EXTRACTED}"

N_CANDIDATES=$(wc -l < "${EXTRACTED}")
echo "Extracted ${N_CANDIDATES} hiConfDeNovo candidates on autosomes"

# ---- Apply read-level filters ----
FINAL_TSV="${WORK_DIR}/results/trio_${CHILD}_final_denovo.tsv"

python3 "${FILTER_PY}" \
    --input "${EXTRACTED}" \
    --output "${FINAL_TSV}" \
    --sire "${SIRE}" \
    --dam "${DAM}" \
    --child "${CHILD}" \
    --min-dp "${MIN_DP}" \
    --ab-min "${AB_MIN}" \
    --ab-max "${AB_MAX}"

echo ""
echo "Final DNMs for ${CHILD} written to: ${FINAL_TSV}"
echo "End: $(date)"
