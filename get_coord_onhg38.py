import sys
import os
import argparse
import requests
import pandas as pd
from tqdm import tqdm


def read_sequences(file_path="sequence.isAnchor.csv"):
    return pd.read_csv(file_path, sep="\t")


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


def update_mainChr_assemblyList(seq, main_chr, assembly_df):
    update_main_chr = main_chr
    update_assembly_df = assembly_df
    if update_main_chr not in seq.name:
        update_main_chr = seq.name.split("_")[0]
        data = get_contig_seq(update_main_chr)
        if "gold" in data:
            update_assembly_df = pd.DataFrame(data["gold"])
        else:
            update_assembly_df = pd.DataFrame()
    return update_main_chr, update_assembly_df


def calculate_region_on_hg38(seq, selected_seq_hg38):
    if (seq.fragStart - selected_seq_hg38.fragStart) >= 0:
        anchor_hg38_start = selected_seq_hg38.chromStart + (seq.fragStart - selected_seq_hg38.fragStart)
    else:
        anchor_hg38_start = selected_seq_hg38.chromStart

    if selected_seq_hg38.fragEnd - seq.fragEnd >= 0:
        anchor_hg38_end = selected_seq_hg38.chromEnd - (selected_seq_hg38.fragEnd - seq.fragEnd)
    else:
        anchor_hg38_end = selected_seq_hg38.chromEnd
    return anchor_hg38_start, anchor_hg38_end


def write_list(lst, output_file):
    header = [
        "pdindex", "index", "name", "chrLength", "chromStart", "chromEnd",
        "frag", "fragStart", "fragEnd", "strand", "fragLength", "isAnchor",
        "hg38_chrom", "hg38_chromStart", "hg38_chromEnd", "hg38_fragStart",
        "hg38_fragEnd", "hg38_fragLength", "anchor_hg38_start", "anchor_hg38_end"
    ]
    df = pd.DataFrame(lst, columns=header)
    df.drop(["pdindex"], axis=1).to_csv(output_file, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Map anchored ALT contig regions to primary hg38 coordinates.")
    parser.add_argument("-i", "--input", default="sequence.isAnchor.csv", help="Input file")
    parser.add_argument("-o", "--output", default="sequence.isAnchor.hg38pos.csv", help="Output file")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Error: {args.input} not found")

    df = read_sequences(args.input)
    df = df[df["anchor"] == 1]
    print(f"Mapping {len(df)} anchored segments to primary chromosomes...")

    seq_full_info_list = []
    main_chr = "chr0"
    assembly_df = pd.DataFrame()

    for seq in tqdm(df.itertuples(), total=len(df)):
        main_chr, assembly_df = update_mainChr_assemblyList(seq, main_chr, assembly_df)
        if assembly_df.empty:
            continue

        seq_hg38 = assembly_df[
            (assembly_df["frag"] == seq.frag) &
            (assembly_df["fragEnd"] > int(seq.fragStart)) &
            (assembly_df["fragStart"] < int(seq.fragEnd)) &
            (assembly_df["fragStart"] <= int(seq.fragStart)) &
            (assembly_df["fragEnd"] >= int(seq.fragEnd))
        ]
        if len(seq_hg38) != 1:
            seq_hg38 = assembly_df[
                (assembly_df["frag"] == seq.frag) &
                (assembly_df["fragEnd"] > int(seq.fragStart)) &
                (assembly_df["fragStart"] < int(seq.fragEnd)) &
                ((assembly_df["fragStart"] <= int(seq.fragStart)) | (assembly_df["fragEnd"] >= int(seq.fragEnd)))
            ]

        if len(seq_hg38) == 0:
            continue

        max_isec_base = 0
        selected_seq_hg38 = seq_hg38.iloc[0]
        for s in seq_hg38.itertuples():
            anchor_hg38_start, anchor_hg38_end = calculate_region_on_hg38(seq, s)
            if anchor_hg38_end - anchor_hg38_start > max_isec_base:
                selected_seq_hg38 = s
        anchor_hg38_start, anchor_hg38_end = calculate_region_on_hg38(seq, selected_seq_hg38)

        update_seq = list(seq)
        update_seq.extend([
            selected_seq_hg38.chrom,
            selected_seq_hg38.chromStart,
            selected_seq_hg38.chromEnd,
            selected_seq_hg38.fragStart,
            selected_seq_hg38.fragEnd,
            selected_seq_hg38.fragEnd - selected_seq_hg38.fragStart,
            anchor_hg38_start,
            anchor_hg38_end
        ])
        seq_full_info_list.append(update_seq)

    write_list(seq_full_info_list, args.output)
    print(f"Saved {len(seq_full_info_list)} mapped regions to {args.output}")


if __name__ == "__main__":
    main()