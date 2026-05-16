#!/bin/bash
#SBATCH --job-name=snpeff_build
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=logs/snpeff_build.%j.out
#SBATCH --error=logs/snpeff_build.%j.err
#SBATCH --mem=8G

# ============================================================================
# Step 6a - SnpEff database build (one-time setup)
# ----------------------------------------------------------------------------
# Builds a custom SnpEff database from a reference FASTA and a MAKER (or
# other) GFF3 annotation, ready for variant annotation in Step 6b.
#
# This step must complete before Step 6b is submitted. Do NOT run as a
# SLURM array — only one process should build the database at a time
# (concurrent builds will race on the data directory).
#
# If your GFF's sequence IDs do not match those in the reference FASTA
# (e.g. MAKER's scaffold names vs. GenBank accessions), provide a
# two-column name-mapping file via SCAFFOLD_NAME_MAP. The mapping is
# applied with a single awk pass to produce a renamed GFF for SnpEff.
#
# The -noCheckCds and -noCheckProtein flags are required for MAKER-style
# annotations where translated CDS sequences may not pass SnpEff's strict
# protein-coding consistency checks.
#
# Inputs:   reference FASTA, GFF3 annotation, optional name-mapping TSV
# Outputs:  SnpEff database in ${SNPEFF_DIR}/data/${GENOME_NAME}/
#           Custom snpEff.config in ${SNPEFF_DIR}/snpEff.config
#
# Usage:    Set REF, GFF, SNPEFF_DIR, GENOME_NAME, and (optionally)
#           SCAFFOLD_NAME_MAP and SPECIES_LABEL below, then submit
#           with sbatch.
# ============================================================================

set -euo pipefail

# ---- User-configurable paths and parameters ----
REF="${REF:-./reference/genome.fna}"
GFF="${GFF:-./reference/annotation.gff}"
SNPEFF_DIR="${SNPEFF_DIR:-./snpeff}"
GENOME_NAME="${GENOME_NAME:-custom_genome}"
SPECIES_LABEL="${SPECIES_LABEL:-Custom genome}"

# Optional: TSV mapping GFF sequence IDs (col 1) -> reference IDs (col 2).
# Leave empty if your GFF already uses the same IDs as your reference.
SCAFFOLD_NAME_MAP="${SCAFFOLD_NAME_MAP:-}"

# ---- Setup ----
SNPEFF_DATA="${SNPEFF_DIR}/data"
GENOME_DIR="${SNPEFF_DATA}/${GENOME_NAME}"
CONFIG_FILE="${SNPEFF_DIR}/snpEff.config"
mkdir -p "${GENOME_DIR}" logs

# ---- Sanity checks ----
for tool in snpEff awk; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "ERROR: required tool not on PATH: ${tool}" >&2
        exit 1
    fi
done
for f in "${REF}" "${GFF}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done

# Skip the build if a complete database already exists.
if [[ -f "${GENOME_DIR}/snpEffectPredictor.bin" ]]; then
    echo "SnpEff database already exists at ${GENOME_DIR}/snpEffectPredictor.bin"
    echo "Delete the file if you want to rebuild."
    exit 0
fi

echo "============================================"
echo "  Building SnpEff database"
echo "  Genome:      ${GENOME_NAME}"
echo "  Species:     ${SPECIES_LABEL}"
echo "  Reference:   ${REF}"
echo "  Annotation:  ${GFF}"
echo "  Output dir:  ${GENOME_DIR}"
echo "  Start:       $(date)"
echo "============================================"

# ---- Stage genome ----
echo "[$(date)] Copying genome FASTA into SnpEff data directory..."
cp "${REF}" "${GENOME_DIR}/sequences.fa"

# ---- Stage GFF (with optional sequence-ID rewriting) ----
if [[ -n "${SCAFFOLD_NAME_MAP}" && -f "${SCAFFOLD_NAME_MAP}" ]]; then
    echo "[$(date)] Rewriting GFF sequence IDs using ${SCAFFOLD_NAME_MAP}..."
    awk -F'\t' 'BEGIN{OFS="\t"}
                NR==FNR{map[$1]=$2; next}
                /^#/{print; next}
                {if($1 in map) $1=map[$1]; print}' \
        "${SCAFFOLD_NAME_MAP}" "${GFF}" \
        > "${GENOME_DIR}/genes.gff"
else
    echo "[$(date)] Copying GFF as-is (no sequence-ID rewriting)..."
    cp "${GFF}" "${GENOME_DIR}/genes.gff"
fi

# ---- Write custom snpEff.config ----
# data.dir MUST be an absolute path; SnpEff resolves it from the working
# directory at runtime, which can move unexpectedly inside SLURM jobs.
SNPEFF_DATA_ABS="$(cd "${SNPEFF_DATA}" && pwd)"

cat > "${CONFIG_FILE}" << EOF
data.dir = ${SNPEFF_DATA_ABS}

${GENOME_NAME}.genome : ${SPECIES_LABEL}
EOF

echo "[$(date)] Custom snpEff.config written to ${CONFIG_FILE}"

# ---- Build the database ----
# -noCheckCds / -noCheckProtein: skip translation-consistency checks
# that often fail on MAKER annotations. Required for many non-model genomes.
echo "[$(date)] Running snpEff build..."
snpEff build \
    -gff3 \
    -c "${CONFIG_FILE}" \
    -noCheckCds \
    -noCheckProtein \
    -v \
    "${GENOME_NAME}"

if [[ ! -f "${GENOME_DIR}/snpEffectPredictor.bin" ]]; then
    echo "ERROR: snpEff build did not produce snpEffectPredictor.bin" >&2
    exit 1
fi

echo ""
echo "============================================"
echo "  SnpEff database build complete."
echo "  Database:  ${GENOME_DIR}/snpEffectPredictor.bin"
echo "  Config:    ${CONFIG_FILE}"
echo ""
echo "  Next step: submit 06b_snpeff_annotate.sh"
echo "  End:       $(date)"
echo "============================================"
