#!/bin/bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COHORT_VCF=""
MAJOR_VCF="${SCRIPT_DIR}/data/majorSet.fix_delins.vcf.gz"
REF_FASTA="Homo_sapiens_assembly38.fasta"
ALT_FILE="${SCRIPT_DIR}/Homo_sapiens_assembly38.fasta.64.alt"
OUT_FASTA="${SCRIPT_DIR}/VGR.consensus.fasta"
THREADS=8
RUN_STEP0=false

usage() {
    cat << EOF
Usage: $(basename "$0") [options]

Build population-specific consensus reference with ALT-contig awareness.

Options:
  --cohort <file>       Input cohort VCF to generate majorSet from scratch
  --major-vcf <file>    Path to precomputed major variants (default: data/majorSet.fix_delins.vcf.gz)
  --ref <file>          Reference FASTA (default: Homo_sapiens_assembly38.fasta)
  --alt <file>          ALT definition file (default: Homo_sapiens_assembly38.fasta.64.alt)
  --output <file>       Output FASTA (default: VGR.consensus.fasta)
  --threads <int>       Number of threads (default: 8)
  -h, --help            Show this help message

Examples:
  # Using precomputed major variants:
  $(basename "$0") --ref /path/to/Homo_sapiens_assembly38.fasta --output VGR.fasta

  # From cohort WGS VCF:
  $(basename "$0") --cohort cohort.vcf.gz --ref Homo_sapiens_assembly38.fasta --output VGR.fasta
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort)
            COHORT_VCF="$2"
            RUN_STEP0=true
            shift 2
            ;;
        --major-vcf)
            MAJOR_VCF="$2"
            shift 2
            ;;
        --ref)
            REF_FASTA="$2"
            shift 2
            ;;
        --alt)
            ALT_FILE="$2"
            shift 2
            ;;
        --output)
            OUT_FASTA="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

export THREADS

if [ ! -f "${REF_FASTA}" ]; then
    echo "Error: Reference FASTA not found: ${REF_FASTA}"
    exit 1
fi

# Step 0: Generate major variants if cohort VCF provided
if [ "${RUN_STEP0}" = true ]; then
    echo "=== Step 0: Generating major variants from cohort VCF ==="
    bash "${SCRIPT_DIR}/generate_major_variant_set.sh" \
        "${COHORT_VCF}" \
        "${REF_FASTA}" \
        "${SCRIPT_DIR}/data"
    MAJOR_VCF="${SCRIPT_DIR}/data/majorSet.fix_delins.vcf.gz"
fi

if [ ! -f "${MAJOR_VCF}" ]; then
    echo "Error: Major variant file not found: ${MAJOR_VCF}"
    echo "Use --major-vcf or provide --cohort to generate it."
    exit 1
fi

# Step 1: Query UCSC API for ALT tiling paths
echo "=== Step 1: Getting ALT contig sequences ==="
python3 "${SCRIPT_DIR}/get_sequences.py" \
    --alt "${ALT_FILE}" \
    --output "${SCRIPT_DIR}/sequence.csv"

# Step 2: Check anchors
echo "=== Step 2: Checking anchored segments ==="
python3 "${SCRIPT_DIR}/check_anchor.py" \
    --input "${SCRIPT_DIR}/sequence.csv" \
    --output "${SCRIPT_DIR}/sequence.isAnchor.csv"

# Step 3: Map coordinates to primary assembly
echo "=== Step 3: Mapping coordinates to primary reference ==="
python3 "${SCRIPT_DIR}/get_coord_onhg38.py" \
    --input "${SCRIPT_DIR}/sequence.isAnchor.csv" \
    --output "${SCRIPT_DIR}/sequence.isAnchor.hg38pos.csv"

# Step 4: Project variants to ALT contigs
echo "=== Step 4: Extracting and projecting variants to ALT contigs ==="
python3 "${SCRIPT_DIR}/extract_major.py" \
    --coords "${SCRIPT_DIR}/sequence.isAnchor.hg38pos.csv" \
    --vcf "${MAJOR_VCF}" \
    --output "${SCRIPT_DIR}/majorSet.fix_delins.extract_anchor.vcf"

# Step 5: Build consensus reference
echo "=== Step 5: Building consensus reference ==="
bash "${SCRIPT_DIR}/build_consensus.sh" \
    "${MAJOR_VCF}" \
    "${SCRIPT_DIR}/majorSet.fix_delins.extract_anchor.vcf" \
    "${REF_FASTA}" \
    "${OUT_FASTA}"

echo "Pipeline finished successfully: ${OUT_FASTA}"
