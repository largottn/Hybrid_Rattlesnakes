#!/bin/bash
#SBATCH --job-name=dnm_trio_call
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-2
#SBATCH --output=logs/dnm_trio_call.%a.%j.out
#SBATCH --error=logs/dnm_trio_call.%a.%j.err
#SBATCH --mem=32G

# ============================================================================
# Step 5a - DNM pipeline: multi-sample trio variant calling
# ----------------------------------------------------------------------------
# For each offspring, runs GATK HaplotypeCaller on the full trio (sire,
# dam, and offspring BAMs) in a single joint call. This is the first step
# of the de novo mutation (DNM) detection pipeline and is distinct from
# the MIE pipeline's per-sample GVCF-mode calling.
#
# Inputs:   reference FASTA (indexed), per-sample BAMs in ${BAM_DIR}
# Outputs:  one raw multi-sample VCF per trio in ${OUT_DIR}/
#           (named trio_${OFFSPRING}_raw.vcf.gz)
#
# Usage:    Set REF, BAM_DIR, OUT_DIR, SIRE, DAM, and OFFSPRING below.
#           Adjust #SBATCH --array=0-N to match (number of offspring - 1).
#           Submit with sbatch.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and samples ----
REF="${REF:-./reference/genome.fna}"
BAM_DIR="${BAM_DIR:-./alignments/final_bams}"
OUT_DIR="${OUT_DIR:-./dnm_output/raw_calls}"

SIRE="${SIRE:-sire_sample}"
DAM="${DAM:-dam_sample}"

# Offspring IDs to process. Length of this array must match the SLURM
# array range above (--array=0-N where N = number of offspring - 1).
OFFSPRING=(child1 child2 child3)

# ---- Setup ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "${OUT_DIR}" logs

# ---- Sanity checks ----
if ! command -v gatk >/dev/null 2>&1; then
    echo "ERROR: gatk not on PATH" >&2
    exit 1
fi
if [[ ! -f "${REF}" || ! -f "${REF}.fai" ]]; then
    echo "ERROR: reference not found or not indexed: ${REF}" >&2
    echo "  Run: samtools faidx ${REF}" >&2
    exit 1
fi

CHILD="${OFFSPRING[${SLURM_ARRAY_TASK_ID:-0}]}"
SIRE_BAM="${BAM_DIR}/${SIRE}.bam"
DAM_BAM="${BAM_DIR}/${DAM}.bam"
CHILD_BAM="${BAM_DIR}/${CHILD}.bam"

for bam in "${SIRE_BAM}" "${DAM_BAM}" "${CHILD_BAM}"; do
    if [[ ! -f "${bam}" ]]; then
        echo "ERROR: BAM not found: ${bam}" >&2
        exit 1
    fi
done

OUT_VCF="${OUT_DIR}/trio_${CHILD}_raw.vcf.gz"

echo "============================================"
echo "  Trio:    ${SIRE} + ${DAM} -> ${CHILD}"
echo "  Output:  ${OUT_VCF}"
echo "  Threads: ${THREADS}"
echo "  Start:   $(date)"
echo "============================================"

if [[ -s "${OUT_VCF}" ]]; then
    echo "Output VCF already exists, skipping."
    exit 0
fi

# ---- Multi-sample trio calling ----
gatk --java-options "-Xmx24g" HaplotypeCaller \
    -R "${REF}" \
    -I "${SIRE_BAM}" \
    -I "${DAM_BAM}" \
    -I "${CHILD_BAM}" \
    -O "${OUT_VCF}" \
    --verbosity ERROR \
    --native-pair-hmm-threads "${THREADS}"

echo ""
echo "============================================"
echo "  Trio ${CHILD} complete."
echo "  End: $(date)"
echo "============================================"
