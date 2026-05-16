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
# Step 1 - Trimmomatic: multi-lane paired-end read trimming
# ----------------------------------------------------------------------------
# Trims adapter sequences and low-quality bases from raw Illumina paired-end
# reads using Trimmomatic v0.39. Designed for samples sequenced across
# multiple lanes, where each lane is delivered as a separate FASTQ pair
# (e.g. ${SAMPLE}_L5_1.fq.gz, ${SAMPLE}_L5_2.fq.gz, ${SAMPLE}_L6_1.fq.gz...).
# Each lane is trimmed independently; lane-level BAMs are merged later
# during the alignment step.
#
# Inputs:   raw lane-split paired-end FASTQs in ${RAW_DIR},
#           named *_L<lane>_1.fq.gz / *_L<lane>_2.fq.gz
# Outputs:  trimmed paired and unpaired FASTQs in ${OUT_DIR}/{paired,unpaired};
#           per-lane logs in ${OUT_DIR}/logs
#
# Usage:    Set RAW_DIR and OUT_DIR below, then submit with sbatch.
#           If running outside SLURM, set THREADS manually.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths ----
RAW_DIR="${RAW_DIR:-./raw_fastq}"
OUT_DIR="${OUT_DIR:-./trimmed}"

# ---- Setup ----
THREADS="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "${OUT_DIR}/paired" "${OUT_DIR}/unpaired" "${OUT_DIR}/logs" logs

# ---- Resolve Trimmomatic command ----
# Prefer the conda wrapper; fall back to the bundled JAR if needed.
if command -v trimmomatic >/dev/null 2>&1; then
    TRIMMER=(trimmomatic)
else
    TRIMJAR="${CONDA_PREFIX:-}/share/trimmomatic-0.39-2/trimmomatic.jar"
    if [[ ! -f "${TRIMJAR}" ]]; then
        echo "ERROR: Trimmomatic not found. Activate a conda env with Trimmomatic installed or set TRIMJAR." >&2
        exit 1
    fi
    TRIMMER=(java -Xmx16g -jar "${TRIMJAR}")
fi

# ---- Adapter file ----
# Uses the TruSeq3-PE.fa bundled with Trimmomatic. Adapter trimming is
# skipped automatically if the file is not found (e.g. non-conda install).
ADAPTER="${CONDA_PREFIX:-}/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
USE_ADAPT=()
if [[ -f "${ADAPTER}" ]]; then
    USE_ADAPT=(ILLUMINACLIP:"${ADAPTER}":2:30:10)
else
    echo "WARNING: TruSeq3-PE.fa adapter file not found; running without ILLUMINACLIP."
fi

# ---- Trim each lane-split read pair ----
shopt -s nullglob
for R1 in "${RAW_DIR}"/*_L[0-9]*_1.fq.gz; do
    base=$(basename "${R1}")
    root="${base%_1.fq.gz}"
    R2="${RAW_DIR}/${root}_2.fq.gz"

    if [[ ! -f "${R2}" ]]; then
        echo "WARNING: mate not found for ${R1}, skipping" >&2
        continue
    fi

    outP1="${OUT_DIR}/paired/${root}_R1.P.trim.fq.gz"
    outU1="${OUT_DIR}/unpaired/${root}_R1.U.trim.fq.gz"
    outP2="${OUT_DIR}/paired/${root}_R2.P.trim.fq.gz"
    outU2="${OUT_DIR}/unpaired/${root}_R2.U.trim.fq.gz"

    echo "=== Trimming ${root} ==="
    "${TRIMMER[@]}" PE -phred33 -threads "${THREADS}" \
        "${R1}" "${R2}" \
        "${outP1}" "${outU1}" \
        "${outP2}" "${outU2}" \
        "${USE_ADAPT[@]}" \
        LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30 \
        2>&1 | tee "${OUT_DIR}/logs/${root}.trim.log"
done

echo "Done. Trimmed reads written to: ${OUT_DIR}/{paired,unpaired}"
