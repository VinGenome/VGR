import sys
import os
import argparse
import pandas as pd
from pysam import VariantFile
from tqdm import tqdm


def find_vcf_file(path=None):
    if path and os.path.exists(path):
        return path
    for default_path in ["data/majorSet.fix_delins.vcf.gz", "majorSet.fix_delins.vcf.gz"]:
        if os.path.exists(default_path):
            return default_path
    raise FileNotFoundError("Could not find input major variants VCF. Please provide --vcf.")


def main():
    parser = argparse.ArgumentParser(description="Project primary major variants onto anchored ALT contigs.")
    parser.add_argument("-c", "--coords", default="sequence.isAnchor.hg38pos.csv", help="Coordinate mapping file")
    parser.add_argument("-v", "--vcf", default=None, help="Input major variants VCF")
    parser.add_argument("-o", "--output", default="majorSet.fix_delins.extract_anchor.vcf", help="Output VCF")
    args = parser.parse_args()

    if not os.path.exists(args.coords):
        sys.exit(f"Error: {args.coords} not found")

    vcf_file = find_vcf_file(args.vcf)
    anchor_df = pd.read_csv(args.coords, sep="\t")

    bcf_in = VariantFile(vcf_file)
    bcf_out = VariantFile(args.output, "w", header=bcf_in.header)

    extracted_count = 0
    for anchor in tqdm(anchor_df.itertuples(), total=len(anchor_df)):
        try:
            records = bcf_in.fetch(anchor.hg38_chrom, int(anchor.anchor_hg38_start), int(anchor.anchor_hg38_end))
        except ValueError:
            continue

        for rec in records:
            changed_rec = rec.copy()
            changed_rec.contig = anchor.name
            if (anchor.hg38_fragLength < anchor.fragLength) and (anchor.fragStart < anchor.hg38_fragStart):
                changed_rec.pos = rec.pos - int(anchor.anchor_hg38_start) + int(anchor.chromStart) + (int(anchor.hg38_fragStart) - int(anchor.fragStart))
            else:
                changed_rec.pos = rec.pos - int(anchor.anchor_hg38_start) + int(anchor.chromStart)

            bcf_out.write(changed_rec)
            extracted_count += 1

    bcf_out.close()
    bcf_in.close()
    print(f"Extracted {extracted_count} variants for ALT contigs to {args.output}")


if __name__ == "__main__":
    main()