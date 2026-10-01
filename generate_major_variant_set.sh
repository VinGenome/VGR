#!/bin/bash
set -eo pipefail

COHORT_VCF=$1
REF_FASTA=${2:-Homo_sapiens_assembly38.fasta}
OUT_DIR=${3:-data}
THREADS=${THREADS:-8}
AF_CUTOFF=${AF_CUTOFF:-0.50}
HWE_PVAL=${HWE_PVAL:-3.4e-6}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    echo "Usage: $0 <cohort.vcf.gz> [reference.fasta] [output_dir]"
    echo ""
    echo "Options via environment variables:"
    echo "  THREADS    Number of threads (default: 8)"
    echo "  AF_CUTOFF  Major allele frequency threshold (default: 0.50)"
    echo "  HWE_PVAL   HWE p-value cutoff (default: 3.4e-6)"
    exit 1
}

if [ -z "${COHORT_VCF}" ] || [ ! -f "${COHORT_VCF}" ]; then
    echo "Error: Input cohort VCF file not found: '${COHORT_VCF}'"
    usage
fi

if [ ! -f "${REF_FASTA}" ]; then
    echo "Error: Reference FASTA not found: '${REF_FASTA}'"
    exit 1
fi

mkdir -p "${OUT_DIR}"
TMP_DIR="$(mktemp -d -p "${OUT_DIR}" tmp_major_XXXXXX)"
trap 'rm -rf "${TMP_DIR}"' EXIT

FINAL_VCF="${OUT_DIR}/majorSet.fix_delins.vcf.gz"

echo "Input VCF: ${COHORT_VCF}"
echo "Reference: ${REF_FASTA}"
echo "Output:    ${FINAL_VCF}"

# 1. Filter PASS variants
echo "Filtering PASS variants..."
bcftools view \
    --threads "${THREADS}" \
    -f PASS \
    -Oz -o "${TMP_DIR}/01.pass.vcf.gz" \
    "${COHORT_VCF}"
tabix -p vcf "${TMP_DIR}/01.pass.vcf.gz"

# 2. Split multiallelic and normalize indels
echo "Normalizing and splitting multiallelic variants..."
bcftools norm \
    --threads "${THREADS}" \
    -m -any \
    -f "${REF_FASTA}" \
    -Oz -o "${TMP_DIR}/02.norm.vcf.gz" \
    "${TMP_DIR}/01.pass.vcf.gz"
tabix -p vcf "${TMP_DIR}/02.norm.vcf.gz"

# 3. Fill tags (AF, AC, AN, HWE)
echo "Annotating tags (AF, AC, AN, HWE)..."
bcftools plugin fill-tags \
    --threads "${THREADS}" \
    "${TMP_DIR}/02.norm.vcf.gz" \
    -Oz -o "${TMP_DIR}/03.tags.vcf.gz" \
    -- -t AF,AC,AN,HWE
tabix -p vcf "${TMP_DIR}/03.tags.vcf.gz"

# 4. Filter by HWE
echo "Filtering by HWE (p > ${HWE_PVAL})..."
bcftools view \
    --threads "${THREADS}" \
    -i "HWE > ${HWE_PVAL}" \
    -Oz -o "${TMP_DIR}/04.hwe.vcf.gz" \
    "${TMP_DIR}/03.tags.vcf.gz"
tabix -p vcf "${TMP_DIR}/04.hwe.vcf.gz"

# 5. Filter major alleles (AF >= AF_CUTOFF)
echo "Extracting major alleles (AF >= ${AF_CUTOFF})..."
bcftools view \
    --threads "${THREADS}" \
    -i "INFO/AF >= ${AF_CUTOFF}" \
    -Oz -o "${TMP_DIR}/05.major.raw.vcf.gz" \
    "${TMP_DIR}/04.hwe.vcf.gz"
tabix -p vcf "${TMP_DIR}/05.major.raw.vcf.gz"

# 6. Standardize indels and resolve overlaps
echo "Standardizing indels and delins..."
python3 "${SCRIPT_DIR}/prepare_delins.py" \
    -i "${TMP_DIR}/05.major.raw.vcf.gz" \
    -o "${TMP_DIR}/06.major.fix_delins.vcf"

bgzip -@ "${THREADS}" -c "${TMP_DIR}/06.major.fix_delins.vcf" > "${FINAL_VCF}"
tabix -p vcf "${FINAL_VCF}"

COUNT=$(bcftools index -n "${FINAL_VCF}")
echo "Done. Generated ${COUNT} major variants in ${FINAL_VCF}"
