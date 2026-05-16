#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-3
#SBATCH --output=logs/align.%a.%j.out
#SBATCH --error=logs/align.%a.%j.err
#SBATCH --mem=32G

# ============================================================================
# Step 2 - Alignment: BWA-MEM, sort, lane merge, MarkDuplicates
# ----------------------------------------------------------------------------
# For each sample in the SAMPLES array, this script:
#   1. Aligns trimmed paired-end reads to the reference with BWA-MEM,
#      attaching read group metadata.
#   2. Sorts the alignment with samtools.
#   3. If the sample was sequenced across multiple lanes, merges the
#      per-lane sorted BAMs into a single sample-level BAM.
#   4. Marks duplicates with Picard MarkDuplicates.
#   5. Writes flagstat and idxstats reports.
#
# Multi-lane samples are detected automatically: if no single-file FASTQ
# (${SAMPLE}_R1.P.trim.fq.gz) is found in READS_DIR, the script looks for
# lane-split files matching ${SAMPLE}_L<n>_R1.P.trim.fq.gz and processes
# each lane independently before merging.
#
# Inputs:   reference FASTA (BWA-indexed) and trimmed paired-end FASTQs
# Outputs:  sorted BAM, merged BAM (if multi-lane), markdup BAM,
#           dup_metrics.txt, flagstat.txt, idxstats.txt
#
# Usage:    Set REF, READS_DIR, OUT_DIR, and SAMPLES below.
#           Adjust #SBATCH --array=0-N to match (number of samples - 1).
#           Submit with sbatch.
#
# Notes:    - The reference FASTA must already be indexed (bwa index, samtools faidx).
#           - Resumable: per-sample steps are skipped if their outputs exist.
# ============================================================================

set -eu

# ---- User-configurable paths and samples ----
REF="${REF:-./reference/genome.fna}"
READS_DIR="${READS_DIR:-./trimmed}"
OUT_DIR="${OUT_DIR:-./alignments}"

# Sample IDs to process. The length of this array must match the SLURM
# array range above (--array=0-N where N = number of samples - 1).
SAMPLES=(sample1 sample2 sample3 sample4)

# ---- Setup ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "${OUT_DIR}" logs

SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "============================================"
echo "  Sample:  ${SAMPLE}"
echo "  Task:    ${SLURM_ARRAY_TASK_ID:-0}"
echo "  Threads: ${THREADS}"
echo "  Start:   $(date)"
echo "============================================"

# ---- Sanity checks ----
if [[ ! -f "${REF}" ]]; then
    echo "ERROR: reference FASTA not found: ${REF}" >&2
    exit 1
fi
for tool in bwa samtools picard; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done

# ---- Step 1: Discover input FASTQs (single-file or multi-lane) ----
SINGLE_R1="${READS_DIR}/${SAMPLE}_R1.P.trim.fq.gz"
shopt -s nullglob

if [[ -f "${SINGLE_R1}" ]]; then
    LANE_R1S=("${SINGLE_R1}")
else
    LANE_R1S=("${READS_DIR}/${SAMPLE}"_L*_R1.P.trim.fq.gz)
fi

if [[ ${#LANE_R1S[@]} -eq 0 ]]; then
    echo "ERROR: No trimmed R1 files found for ${SAMPLE} in ${READS_DIR}" >&2
    echo "  Expected: ${SAMPLE}_R1.P.trim.fq.gz  or  ${SAMPLE}_L<n>_R1.P.trim.fq.gz" >&2
    exit 1
fi

echo ""
echo "=== Step 1: Aligning ${#LANE_R1S[@]} read file(s) for ${SAMPLE} ==="

# ---- Step 2: Align each R1 (one alignment per lane) ----
SORTED_BAMS=()
for R1 in "${LANE_R1S[@]}"; do
    R2="${R1/_R1.P./_R2.P.}"
    bn=$(basename "${R1}")

    # Derive a lane label from the filename if present, else treat as single-lane
    if [[ "${bn}" =~ ^${SAMPLE}_(L[0-9]+)_R1\.P\.trim\.fq\.gz$ ]]; then
        LANE="${BASH_REMATCH[1]}"
        RG_ID="${SAMPLE}_${LANE}"
    else
        RG_ID="${SAMPLE}"
    fi

    LANE_BAM="${OUT_DIR}/${RG_ID}.sorted.bam"

    if [[ ! -f "${R2}" ]]; then
        echo "ERROR: mate not found for ${R1}" >&2
        exit 1
    fi

    if [[ -f "${LANE_BAM}" ]]; then
        echo "  ${RG_ID}: sorted BAM exists, skipping."
    else
        echo "  Aligning ${RG_ID}..."
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${RG_ID}\tSM:${SAMPLE}\tLB:${SAMPLE}_lib1\tPL:ILLUMINA\tPU:${RG_ID}" \
            "${REF}" "${R1}" "${R2}" \
        | samtools sort -@ "${THREADS}" -o "${LANE_BAM}" -
        samtools index "${LANE_BAM}"
        echo "  ${RG_ID}: alignment complete."
    fi

    SORTED_BAMS+=("${LANE_BAM}")
done

# ---- Step 3: Merge lanes (if multi-lane) or rename single lane ----
MERGED="${OUT_DIR}/${SAMPLE}.sorted.bam"

if [[ ${#SORTED_BAMS[@]} -gt 1 ]]; then
    echo ""
    echo "=== Step 2: Merging ${#SORTED_BAMS[@]} lanes for ${SAMPLE} ==="
    if [[ -f "${MERGED}" ]]; then
        echo "  Merged BAM exists, skipping."
    else
        samtools merge -@ "${THREADS}" "${MERGED}" "${SORTED_BAMS[@]}"
        samtools index "${MERGED}"
        echo "  Merge complete."
    fi
elif [[ "${SORTED_BAMS[0]}" != "${MERGED}" ]]; then
    # Single lane with a lane-suffixed filename: symlink to the canonical name
    # so the downstream MarkDuplicates step finds a consistent input.
    ln -sf "$(basename "${SORTED_BAMS[0]}")" "${MERGED}"
    ln -sf "$(basename "${SORTED_BAMS[0]}").bai" "${MERGED}.bai"
fi

# ---- Step 4: Mark duplicates ----
echo ""
echo "=== Step 3: Marking duplicates ==="

MARKDUP="${OUT_DIR}/${SAMPLE}.markdup.bam"
METRICS="${OUT_DIR}/${SAMPLE}.dup_metrics.txt"

if [[ -f "${MARKDUP}" ]]; then
    echo "  ${SAMPLE}: markdup BAM exists, skipping."
else
    picard -Xmx24g MarkDuplicates \
        I="${MERGED}" \
        O="${MARKDUP}" \
        M="${METRICS}" \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT
    echo "  Duplicates marked."
fi

# ---- Step 5: Alignment stats ----
echo ""
echo "=== Step 4: Generating stats ==="
samtools flagstat "${MARKDUP}" > "${OUT_DIR}/${SAMPLE}.flagstat.txt"
samtools idxstats "${MARKDUP}" > "${OUT_DIR}/${SAMPLE}.idxstats.txt"

echo ""
echo "============================================"
echo "  ${SAMPLE} complete!"
echo "  Final BAM:    ${MARKDUP}"
echo "  Dup metrics:  ${METRICS}"
echo "  Flagstat:     ${OUT_DIR}/${SAMPLE}.flagstat.txt"
echo "  End:          $(date)"
echo "============================================"
