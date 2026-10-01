# Vietnamese Genome Reference (VGR)

Pipeline to construct population-specific consensus reference genomes (GRCh38-based) with ALT-contig awareness.

Standard human reference assemblies (such as GRCh38) carry alleles that may be minor or non-existent in specific populations, introducing reference allele bias during read alignment and variant calling. VGR incorporates high-frequency major alleles (allele frequency >= 50%) into GRCh38 across both canonical chromosomes (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`) and alternate contigs (`chr*_alt`). Homologous primary chromosome coordinates are projected onto anchored ALT contigs to ensure compatibility with ALT-aware aligners (such as BWA-MEM using `bwa-postalt.js`).

---

## Workflow Overview

```text
Population Cohort VCF (e.g., VN1K joint callset)
   │
   ▼ Step 0: generate_major_variant_set.sh & prepare_delins.py
   │  - Filter PASS variants (GATK VQSR tranche)
   │  - Split multiallelic sites & normalize indels (bcftools norm)
   │  - Annotate AF & filter HWE (p > 3.4e-6)
   │  - Retain major alleles (AF >= 0.50)
   │  - Trim indels/delins & resolve overlapping variant collisions
   │
   ▼ Major variant set (data/majorSet.fix_delins.vcf.gz)
   │
   ├─────────────────────────────────────────────────┐
   ▼ Primary chromosomes                             ▼ ALT contig projection (Steps 1–4)
   Primary major variants                            - get_sequences.py (UCSC tiling path)
   │                                                 - check_anchor.py (anchored contig check)
   │                                                 - get_coord_onhg38.py (homologous coordinates)
   │                                                 - extract_major.py (project to ALT contigs)
   │                                                 │
   │                                                 ▼
   │                                                 Projected ALT major variants
   │                                                 │
   └────────────────────────┬────────────────────────┘
                            ▼ Step 5: build_consensus.sh
                            - Merge variant sets (bcftools concat)
                            - Apply consensus alleles (bcftools consensus)
                            - Index FASTA & link ALT indices (.fai, .dict, .alt)
                            │
                            ▼
                    VGR.consensus.fasta
```

---

## Requirements

- Linux environment
- Python >= 3.8 (`pandas`, `requests`, `pysam`, `tqdm`)
- `bcftools` (>= 1.12)
- `samtools` (>= 1.12)
- `tabix` / `bgzip`
- `k8` (optional, required if running `bwa-postalt.js`)

Dependencies can be installed via Conda:

```bash
conda create -n vgr -c bioconda -c conda-forge python=3.10 pandas requests pysam tqdm bcftools samtools htslib -y
conda activate vgr
```

---

## Step-by-Step Guide

### Step 0: Prepare Population Major Variant Set (Optional)

If starting from a raw joint-called cohort VCF (e.g. from GATK HaplotypeCaller):

```bash
bash generate_major_variant_set.sh cohort.vcf.gz Homo_sapiens_assembly38.fasta data/
```

This step executes the following operations:
1. **Quality filtering**: Retains only `PASS` variants from VQSR (99.0% tranche).
2. **Multiallelic decomposition & normalization**: Runs `bcftools norm -m -any -f reference.fasta` to decompose multiallelic variants into biallelic records and left-align indels against the reference genome.
3. **Hardy-Weinberg Equilibrium (HWE) filtering**: Annotates AC, AN, AF, and HWE p-values (`bcftools plugin fill-tags`). Filters out sites with severe HWE deviation (p <= 3.4e-6) to remove potential genotype calling or sequencing artifacts.
4. **Major allele selection**: Selects variants with allele frequency AF >= 0.50 (where the alternate allele is the majority in the target population).
5. **InDel & DELINS standardization (`prepare_delins.py`)**:
   - Trims shared prefix and suffix bases between REF and ALT alleles to ensure atomic representation.
   - Sorts records by chromosome and coordinate.
   - Resolves overlapping variant loci by retaining non-overlapping sites, preventing allele collision errors during `bcftools consensus`.

*Note: For the Vietnamese population reference, the processed major variant set is already provided at `data/majorSet.fix_delins.vcf.gz`.*

### Step 1: Query ALT Contig Tiling Paths

Extracts ALT contigs from `Homo_sapiens_assembly38.fasta.64.alt` and queries the UCSC Genome Browser API (`gold` track) for fragment accessions and coordinates:

```bash
python3 get_sequences.py \
    --alt Homo_sapiens_assembly38.fasta.64.alt \
    --output sequence.csv
```

### Step 2: Identify Anchored Contig Segments

Cross-references ALT fragment accessions against primary chromosome assembly tracks to identify anchored sequence segments:

```bash
python3 check_anchor.py \
    --input sequence.csv \
    --output sequence.isAnchor.csv
```

### Step 3: Map ALT Segments to Primary GRCh38 Coordinates

Computes homologous primary chromosome coordinate intervals (`anchor_hg38_start`, `anchor_hg38_end`) for each anchored ALT segment:

```bash
python3 get_coord_onhg38.py \
    --input sequence.isAnchor.csv \
    --output sequence.isAnchor.hg38pos.csv
```

### Step 4: Project Primary Major Variants onto ALT Contigs

Extracts primary chromosome major variants within the homologous anchored regions and translates their positions onto the respective ALT contig coordinates:

```bash
python3 extract_major.py \
    --coords sequence.isAnchor.hg38pos.csv \
    --vcf data/majorSet.fix_delins.vcf.gz \
    --output majorSet.fix_delins.extract_anchor.vcf
```

### Step 5: Build Final Consensus Genome

Combines primary and ALT variants, applies them to the GRCh38 reference FASTA, and creates index files (`.fai`, `.dict`, `.alt`):

```bash
bash build_consensus.sh \
    data/majorSet.fix_delins.vcf.gz \
    majorSet.fix_delins.extract_anchor.vcf \
    Homo_sapiens_assembly38.fasta \
    VGR.consensus.fasta
```

---

## One-Command Execution

The complete workflow can also be executed with `run_pipeline.sh`:

```bash
# Using precomputed major variants:
bash run_pipeline.sh \
    --ref Homo_sapiens_assembly38.fasta \
    --output VGR.consensus.fasta

# Starting from raw cohort VCF:
bash run_pipeline.sh \
    --cohort cohort.vcf.gz \
    --ref Homo_sapiens_assembly38.fasta \
    --output VGR.consensus.fasta
```

---

## Downstream Read Alignment (ALT-Aware BWA-MEM)

To align short reads against the consensus reference with ALT-contig awareness:

```bash
# 1. Index consensus reference (if not already indexed)
bwa index VGR.consensus.fasta

# 2. Align reads and process ALT contigs with bwa-postalt.js
bwa mem -t 16 -K 100000000 -Y VGR.consensus.fasta read_1.fastq.gz read_2.fastq.gz | \
    k8 bwa-postalt.js -p Homo_sapiens_assembly38.fasta.64.alt > aligned.sam
```

---

## Repository Files

- `generate_major_variant_set.sh`: Shell pipeline for cohort filtering, HWE calculation, and major variant extraction.
- `prepare_delins.py`: Normalizes InDels/DELINS and resolves overlapping variant positions.
- `get_sequences.py`: Queries UCSC Genome Browser API for ALT contig tiling paths.
- `check_anchor.py`: Identifies anchored ALT contig segments.
- `get_coord_onhg38.py`: Maps ALT segments to homologous primary chromosome coordinates.
- `extract_major.py`: Remaps primary major variants onto corresponding ALT contigs.
- `build_consensus.sh`: Combines variant callsets and builds the consensus reference FASTA with indices.
- `run_pipeline.sh`: End-to-end automated runner.
- `bwa-postalt.js`: BWA-KIT post-processing script for ALT-aware read alignment.
- `data/majorSet.fix_delins.vcf.gz`: Curated Vietnamese population major variant callset.
- `Homo_sapiens_assembly38.fasta.64.alt`: Alternate contig definition file.

---

## Citation

VinBigdata / VinGenome Bioinformatics Team.
