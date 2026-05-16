#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Step 0 - Reference preparation: integrate the W chromosome
# ----------------------------------------------------------------------------
# Concatenates a primary reference genome assembly with a separate
# W chromosome assembly, sanity-checks for duplicate sequence names,
# and indexes the resulting combined FASTA with samtools and bwa.
#
# For the rattlesnake hybrid project this was used to add the Crotalus
# viridis W chromosome assembly (Schield et al. 2022; GBE; GenBank
# GCA_024760675.1; BioProject PRJNA853338) onto the autosomes + Z of
# the main CroVir_3.0 reference (Schield et al. 2019; BJLS; GCA_003400415.2),
# producing a single FASTA suitable for whole-genome resequencing of
# female samples.
#
# Downloads required before running:
#   - Reference assembly (e.g. CroVir_3.0): NCBI Datasets or FTP, save
#     the .fna file and point REF_GENOME to it.
#   - W chromosome assembly: NCBI GenBank accession GCA_024760675.1,
#     download as a .fna and point W_FASTA to it.
#   - (Optional) W chromosome annotation: cloned from a public repo at
#     runtime if W_ANNOTATION_REPO is set.
#
# This script is run once, on a login or interactive node. It is not a
# SLURM job — the indexing step is the only slow part (~10-30 minutes
# for bwa index on a snake-sized genome).
#
# Usage:    Set the paths below, then run as a regular bash script.
# ============================================================================

# ---- User-configurable paths ----
REF_GENOME="${REF_GENOME:-./reference/main_assembly.fna}"
W_FASTA="${W_FASTA:-./reference/W_chromosome.fna}"
OUT_DIR="${OUT_DIR:-./reference}"
COMBINED_NAME="${COMBINED_NAME:-combined_reference.fna}"

# Optional: clone a public repo containing the W chromosome annotation.
# Leave empty to skip. The script copies a GFF out of the repo if found.
W_ANNOTATION_REPO="${W_ANNOTATION_REPO:-}"
W_ANNOTATION_PATH_IN_REPO="${W_ANNOTATION_PATH_IN_REPO:-}"

# ---- Setup ----
mkdir -p "${OUT_DIR}"
COMBINED_GENOME="${OUT_DIR}/${COMBINED_NAME}"

# ---- Sanity checks ----
for f in "${REF_GENOME}" "${W_FASTA}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required input not found: ${f}" >&2
        exit 1
    fi
done

# ---- Step 1: (optional) clone the W chromosome annotation repo ----
W_GFF=""
if [[ -n "${W_ANNOTATION_REPO}" ]]; then
    echo "=== Step 1: Downloading W chromosome annotation ==="
    REPO_NAME=$(basename "${W_ANNOTATION_REPO}" .git)
    if [[ ! -d "${REPO_NAME}" ]]; then
        echo "  Cloning ${W_ANNOTATION_REPO}..."
        git clone "${W_ANNOTATION_REPO}"
    else
        echo "  Repository ${REPO_NAME} already present, skipping clone."
    fi

    if [[ -n "${W_ANNOTATION_PATH_IN_REPO}" && -f "${REPO_NAME}/${W_ANNOTATION_PATH_IN_REPO}" ]]; then
        cp "${REPO_NAME}/${W_ANNOTATION_PATH_IN_REPO}" "${OUT_DIR}/"
        W_GFF="${OUT_DIR}/$(basename "${W_ANNOTATION_PATH_IN_REPO}")"
        echo "  Copied W annotation GFF to ${W_GFF}"
    else
        echo "  W chromosome GFF not found at expected path; skipping."
    fi
fi

# ---- Step 2: check for duplicate sequence names ----
echo ""
echo "=== Step 2: Checking for duplicate sequence names ==="

REF_HEADERS=$(grep -c "^>" "${REF_GENOME}")
W_HEADERS=$(grep -c "^>" "${W_FASTA}")
echo "  Reference sequences:    ${REF_HEADERS}"
echo "  W chromosome sequences: ${W_HEADERS}"
echo "  First W chromosome headers:"
grep "^>" "${W_FASTA}" | head -5 | sed 's/^/    /'

TMP=$(mktemp -d)
trap 'rm -rf "${TMP}"' EXIT
grep "^>" "${REF_GENOME}" | sed 's/^>//' | cut -d' ' -f1 | sort > "${TMP}/ref.txt"
grep "^>" "${W_FASTA}"    | sed 's/^>//' | cut -d' ' -f1 | sort > "${TMP}/w.txt"
COLLISIONS=$(comm -12 "${TMP}/ref.txt" "${TMP}/w.txt" | wc -l)

if [[ "${COLLISIONS}" -gt 0 ]]; then
    echo ""
    echo "  WARNING: ${COLLISIONS} duplicate header(s) found between assemblies!"
    echo "  Overlapping names (first 10):"
    comm -12 "${TMP}/ref.txt" "${TMP}/w.txt" | head -10 | sed 's/^/    /'
    echo ""
    echo "  Rename these before proceeding, or downstream tools may"
    echo "  silently pick the wrong sequence."
    exit 1
else
    echo "  No collisions found. Safe to concatenate."
fi

# ---- Step 3: concatenate ----
echo ""
echo "=== Step 3: Creating combined reference ==="
cat "${REF_GENOME}" "${W_FASTA}" > "${COMBINED_GENOME}"

TOTAL_HEADERS=$(grep -c "^>" "${COMBINED_GENOME}")
EXPECTED=$((REF_HEADERS + W_HEADERS))
echo "  Combined sequences: ${TOTAL_HEADERS} (expected: ${EXPECTED})"

if [[ "${TOTAL_HEADERS}" -ne "${EXPECTED}" ]]; then
    echo "  ERROR: sequence count mismatch — check inputs." >&2
    exit 1
fi
echo "  Output: ${COMBINED_GENOME}"

# ---- Step 4: index ----
echo ""
echo "=== Step 4: Indexing combined genome ==="

if command -v samtools >/dev/null 2>&1; then
    echo "  samtools faidx..."
    samtools faidx "${COMBINED_GENOME}"
    echo "  Created ${COMBINED_GENOME}.fai"
else
    echo "  WARNING: samtools not found. Run manually:"
    echo "    samtools faidx ${COMBINED_GENOME}"
fi

if command -v bwa >/dev/null 2>&1; then
    echo "  bwa index (this can take 10-30 minutes for snake-sized genomes)..."
    bwa index "${COMBINED_GENOME}"
    echo "  bwa index complete."
else
    echo "  WARNING: bwa not found. Run manually:"
    echo "    bwa index ${COMBINED_GENOME}"
fi

# ---- Summary ----
echo ""
echo "============================================"
echo "  Reference preparation complete."
echo "============================================"
echo "  Combined genome: ${COMBINED_GENOME}"
echo "  Sequences:       ${TOTAL_HEADERS}"
echo "                   (${REF_HEADERS} from main assembly + ${W_HEADERS} from W)"
if [[ -n "${W_GFF}" ]]; then
    echo "  W annotation:    ${W_GFF}"
fi
echo ""
echo "  Before running GATK, also create a sequence dictionary:"
echo "    gatk CreateSequenceDictionary -R ${COMBINED_GENOME}"
echo "============================================"
