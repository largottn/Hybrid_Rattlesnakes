#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/trimmomatic.%j.out
#SBATCH --error=logs/trimmomatic.%j.err
#SBATCH --mem=32G

# ============================================================================
# Step 1 - Trimmomatic: paired-end read trimming with TruSeq adapter removal
# ----------------------------------------------------------------------------
# Trims adapter sequences and low-quality bases from raw Illumina paired-end
# reads using Trimmomatic v0.39, producing paired and unpaired output files.
# Adapter sequences are written to a temporary FASTA file at runtime.
#
# Inputs:   raw paired-end FASTQs (R1/R2) in ${RAW_DIR}
# Outputs:  trimmed paired and unpaired FASTQs in ${OUT_DIR}
#
# Usage:    Set RAW_DIR, OUT_DIR, and SAMPLES below, then submit with sbatch.
#           If running outside SLURM, set THREADS manually.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and samples ----
# Directory containing the raw FASTQs. Files are expected to be named
# ${SAMPLE_PREFIX}${SAMPLE}_R1.fastq.gz / _R2.fastq.gz (see below).
RAW_DIR="${RAW_DIR:-./raw_fastq}"

# Where trimmed reads will be written.
OUT_DIR="${OUT_DIR:-./trimmed}"

# Optional prefix that precedes each sample ID in the raw filenames.
# Set to "" if your raw files are just ${SAMPLE}_R1.fastq.gz.
SAMPLE_PREFIX="${SAMPLE_PREFIX:-}"

# Sample IDs to process. Adjust to match your data.
SAMPLES=(sample1 sample2 sample3 sample4)

# ---- TruSeq adapter sequences (Illumina TruSeq DNA/RNA, R1 and R2) ----
READ1_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
READ2_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# ---- Setup ----
THREADS="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "${OUT_DIR}" logs

# Write adapter FASTA to a temp file for ILLUMINACLIP
ADAPTERS_FILE=$(mktemp --suffix=.fa)
cat > "${ADAPTERS_FILE}" << EOF
>PrefixPE/1
${READ1_ADAPTER}
>PrefixPE/2
${READ2_ADAPTER}
EOF
trap 'rm -f "${ADAPTERS_FILE}"' EXIT

# ---- Trim each sample ----
for SAMPLE in "${SAMPLES[@]}"; do
    echo "=== Processing ${SAMPLE} ==="

    RAW_R1="${RAW_DIR}/${SAMPLE_PREFIX}${SAMPLE}_R1.fastq.gz"
    RAW_R2="${RAW_DIR}/${SAMPLE_PREFIX}${SAMPLE}_R2.fastq.gz"

    TRIM_R1="${OUT_DIR}/${SAMPLE}_R1.P.trim.fq.gz"
    TRIM_R2="${OUT_DIR}/${SAMPLE}_R2.P.trim.fq.gz"
    UNPAIRED_R1="${OUT_DIR}/${SAMPLE}_R1.U.trim.fq.gz"
    UNPAIRED_R2="${OUT_DIR}/${SAMPLE}_R2.U.trim.fq.gz"

    if [[ ! -f "${RAW_R1}" || ! -f "${RAW_R2}" ]]; then
        echo "  WARNING: raw FASTQ not found for ${SAMPLE}, skipping"
        continue
    fi

    trimmomatic PE \
        -threads "${THREADS}" \
        "${RAW_R1}" "${RAW_R2}" \
        "${TRIM_R1}" "${UNPAIRED_R1}" \
        "${TRIM_R2}" "${UNPAIRED_R2}" \
        ILLUMINACLIP:"${ADAPTERS_FILE}":2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36

    echo "  ${SAMPLE} complete."
done

echo "Done. Trimmed reads written to: ${OUT_DIR}"
