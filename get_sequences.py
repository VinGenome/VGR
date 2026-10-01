import sys
import os
import argparse
import requests
import pandas as pd
from tqdm import tqdm


def read_hg38alt(alt_file="Homo_sapiens_assembly38.fasta.64.alt"):
    alt_contigs_list = []
    with open(alt_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) > 0 and "_alt" in cols[0]:
                alt_contigs_list.append(cols[0])
    return alt_contigs_list


def get_contig_seq(chr_name, max_retries=3):
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gold;chrom={chr_name}"
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                return resp.json()
        except requests.RequestException:
            if attempt == max_retries - 1:
                raise
    return {}


def parse_info(index, contig_json):
    seqs = []
    if "chrom" not in contig_json or "gold" not in contig_json:
        return seqs
    name = contig_json["chrom"]
    length = contig_json["end"]
    for seq in contig_json["gold"]:
        seqs.append([
            index + 1,
            name,
            length,
            seq["chromStart"],
            seq["chromEnd"],
            seq["frag"],
            seq["fragStart"],
            seq["fragEnd"],
            seq.get("strand", "+"),
            seq["fragEnd"] - seq["fragStart"],
        ])
    return seqs


def main():
    parser = argparse.ArgumentParser(description="Query UCSC API for ALT contig tiling path.")
    parser.add_argument("-a", "--alt", default="Homo_sapiens_assembly38.fasta.64.alt", help="Path to alt contig file")
    parser.add_argument("-o", "--output", default="sequence.csv", help="Output path (default: sequence.csv)")
    args = parser.parse_args()

    if not os.path.exists(args.alt):
        sys.exit(f"Error: {args.alt} not found")

    alt_contigs = read_hg38alt(args.alt)
    print(f"Querying UCSC gold track for {len(alt_contigs)} ALT contigs...")

    all_seq_list = []
    for index, chr_name in tqdm(enumerate(alt_contigs), total=len(alt_contigs)):
        contig_json = get_contig_seq(chr_name)
        seqs = parse_info(index, contig_json)
        all_seq_list.extend(seqs)

    df = pd.DataFrame(all_seq_list)
    df.to_csv(args.output, index=False, header=False, sep="\t")
    print(f"Saved {len(all_seq_list)} sequences to {args.output}")


if __name__ == "__main__":
    main()