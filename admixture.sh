#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=admixture.%A.out
#SBATCH --error=admixture.%A.err
#SBATCH --mem=32G
#
# admixture.sh
#
# Full ADMIXTURE workflow:
#   1. Recode PLINK .bim chromosomes for ADMIXTURE compatibility
#   2. Run ADMIXTURE across K = 2-5 with cross-validation
#   3. Summarise CV error per K (lowest = best-supported K)
#
# ADMIXTURE ignores physical position and uses only the genotype matrix, but it
# rejects .bim chromosome codes it does not recognise (scaffold names, etc.), so
# step 1 sets the chromosome column to a placeholder ("1").
#
# INPUT : ${IN_PREFIX}.{bed,bim,fam}  (LD-pruned PLINK set; see note)
# OUTPUT: ${OUT_PREFIX}.{bed,bim,fam} (ADMIXTURE-ready inputs)
#         ${BASE_OUT_DIR}/run_K_<K>/  with .Q, .P, and CV.out for each K
#
# NOTE: this assumes the VCF has already been converted to PLINK binary format
#       and LD-pruned upstream (producing population_pruned.*). That conversion
#       and pruning step is not part of this script.
#
# Dependencies: ADMIXTURE (v1.3.0), awk, grep, coreutils
# Submit with: sbatch admixture.sh

# Exit on any error
set -euo pipefail

##1. SETTINGS  (edit these paths for your system)
CORES=$SLURM_CPUS_PER_TASK
IN_DIR="/path/to/plink_files"             # directory containing the PLINK input set
IN_PREFIX="${IN_DIR}/population_pruned"    # LD-pruned PLINK input prefix
OUT_PREFIX="${IN_DIR}/population_admix"     # chromosome-corrected ADMIXTURE input prefix
IN_FILE="${OUT_PREFIX}.bed"
BASE_OUT_DIR="/path/to/admixture_output"   # parent directory for per-K results

## 2. PREPARE ADMIXTURE INPUT
echo "----------------------------------------------------"
echo "Creating ADMIXTURE-compatible files..."
echo "   >> Forcing all chromosomes in .bim file to '1'"

# Rewrite column 1 (chromosome) of the .bim to "1"; leave all other columns intact.
awk 'BEGIN {OFS="\t"} {$1 = 1; print $0}' "${IN_PREFIX}.bim" > "${OUT_PREFIX}.bim"

# The .bed  and .fam  do not need to change.
cp "${IN_PREFIX}.bed" "${OUT_PREFIX}.bed"
cp "${IN_PREFIX}.fam" "${OUT_PREFIX}.fam"

echo "   >> Input files ready: ${OUT_PREFIX}.bed, .bim, .fam"

##3. RUN ADMIXTURE
mkdir -p "$BASE_OUT_DIR"

echo "----------------------------------------------------"
echo "Starting ADMIXTURE analysis loop for K=2 through 5..."
echo "   >> Input file: ${IN_FILE}"

# Loop from K=2 to K=5 (edit the range here to scan more values of K)
for K in {2..5}; do

    # K-specific output directory
    OUT_DIR="${BASE_OUT_DIR}/run_K_${K}"
    mkdir -p "$OUT_DIR"

    echo "---"
    echo "Running K=$K... (Threads: $CORES)"
    echo "   >> Output directory: ${OUT_DIR}"

    # Run inside the output directory so ADMIXTURE writes its .Q/.P files there.
    cd "$OUT_DIR"

    # --cv calculates cross-validation error, which ADMIXTURE prints to stdout.
    # ADMIXTURE writes no log of its own, so tee the output to CV.out in each run
    # directory (while still echoing to the SLURM log) for the summary below.
    admixture --cv -j"$CORES" "$IN_FILE" "$K" 2>&1 | tee "CV.out"

    # Return to the submission directory.
    cd "$SLURM_SUBMIT_DIR"
done

## --- 4. SUMMARISE CROSS-VALIDATION ERROR ---
echo "----------------------------------------------------"
echo "Success! All ADMIXTURE runs are complete."
echo "---"
echo "CV errors (lowest indicates the best-supported K):"
grep -h 'CV error' "${BASE_OUT_DIR}"/run_K_*/CV.out | sort -V
