#!/bin/bash
set -eo pipefail

PRIMARY_VCF=${1:-data/majorSet.fix_delins.vcf.gz}
ALT_VCF=${2:-majorSet.fix_delins.extract_anchor.vcf}
BASE_FASTA=${3:-Homo_sapiens_assembly38.fasta}
OUT_FASTA=${4:-VGR.consensus.fasta}
THREADS=${THREADS:-8}

usage() {
    echo "Usage: $0 [primary_variants.vcf.gz] [alt_variants.vcf] [reference.fasta] [output.fasta]"
    echo ""
    echo "Options via environment variables:"
    echo "  THREADS  Number of threads (default: 8)"
    exit 1
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    usage
fi

if [ ! -f "${PRIMARY_VCF}" ]; then
    echo "Error: Primary variant VCF not found: ${PRIMARY_VCF}"
    exit 1
fi

if [ ! -f "${BASE_FASTA}" ]; then
    echo "Error: Reference FASTA not found: ${BASE_FASTA}"
    exit 1
fi

OUT_DIR="$(dirname "${OUT_FASTA}")"
mkdir -p "${OUT_DIR}"
TMP_DIR="$(mktemp -d -p "${OUT_DIR}" tmp_consensus_XXXXXX)"
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "Primary VCF: ${PRIMARY_VCF}"
echo "ALT VCF:     ${ALT_VCF}"
echo "Reference:   ${BASE_FASTA}"
echo "Output:      ${OUT_FASTA}"

# 1. Process ALT contig variants if present
MERGE_INPUTS=()
if [ -f "${ALT_VCF}" ]; then
    echo "Indexing ALT contig variants..."
    if [[ "${ALT_VCF}" == *.gz ]]; then
        ALT_VCF_GZ="${ALT_VCF}"
    else
        ALT_VCF_GZ="${TMP_DIR}/alt_variants.vcf.gz"
        bgzip -@ "${THREADS}" -c "${ALT_VCF}" > "${ALT_VCF_GZ}"
    fi
    tabix -f -p vcf "${ALT_VCF_GZ}"
    MERGE_INPUTS+=("${ALT_VCF_GZ}")
fi

# 2. Exclude random and unplaced contigs
echo "Filtering primary variants..."
PRIMARY_CLEAN="${TMP_DIR}/primary_clean.vcf.gz"
bcftools view \
    --threads "${THREADS}" \
    -t ^$(grep -E "_random|^chrUn_" "${BASE_FASTA}.fai" 2>/dev/null | cut -f1 | paste -sd, - || echo "") \
    "${PRIMARY_VCF}" \
    -Oz -o "${PRIMARY_CLEAN}"
tabix -f -p vcf "${PRIMARY_CLEAN}"

# 3. Concatenate variants
ALL_MAJOR_VCF="${TMP_DIR}/all_major_variants.vcf.gz"
if [ ${#MERGE_INPUTS[@]} -gt 0 ]; then
    echo "Merging primary and ALT variants..."
    bcftools concat \
        --threads "${THREADS}" \
        -a "${PRIMARY_CLEAN}" "${MERGE_INPUTS[@]}" \
        -Oz -o "${ALL_MAJOR_VCF}"
else
    cp "${PRIMARY_CLEAN}" "${ALL_MAJOR_VCF}"
fi
tabix -f -p vcf "${ALL_MAJOR_VCF}"

# 4. Apply variants to reference
echo "Building consensus reference with bcftools..."
bcftools consensus \
    -f "${BASE_FASTA}" \
    "${ALL_MAJOR_VCF}" \
    > "${OUT_FASTA}"

# 5. Index reference
echo "Indexing consensus FASTA..."
samtools faidx "${OUT_FASTA}"
samtools dict -o "${OUT_FASTA%.fasta}.dict" "${OUT_FASTA}" 2>/dev/null || true

# Copy ALT index if available
BASE_DIR="$(dirname "${BASE_FASTA}")"
if [ -f "${BASE_DIR}/Homo_sapiens_assembly38.fasta.64.alt" ]; then
    cp "${BASE_DIR}/Homo_sapiens_assembly38.fasta.64.alt" "${OUT_FASTA}.64.alt"
    cp "${BASE_DIR}/Homo_sapiens_assembly38.fasta.64.alt" "${OUT_FASTA}.alt"
elif [ -f "Homo_sapiens_assembly38.fasta.64.alt" ]; then
    cp "Homo_sapiens_assembly38.fasta.64.alt" "${OUT_FASTA}.64.alt"
    cp "Homo_sapiens_assembly38.fasta.64.alt" "${OUT_FASTA}.alt"
fi

echo "Done. Reference built: ${OUT_FASTA}"
