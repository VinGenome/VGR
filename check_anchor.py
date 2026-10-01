import sys
import os
import argparse
import requests
import pandas as pd
from tqdm import tqdm


def read_sequences(file_path="sequence.csv"):
    header = [
        "index", "name", "chrLength", "chromStart", "chromEnd",
        "frag", "fragStart", "fragEnd", "fragStrand", "fragLength"
    ]
    return pd.read_csv(file_path, names=header, sep="\t")


def get_anchors_from_file(file_name):
    anchors_list = []
    with open(file_name, "r") as f:
        lines = f.readlines()

    if len(lines) < 6:
        return anchors_list

    number_of_anchors = int((len(lines) - 1) / 6) - 1
    chr_name = lines[0].replace("\n", "").strip()

    for i in range(number_of_anchors):
        try:
            frag = lines[(i + 2) * 6].replace("\n", "").split(": ")[-1].replace('"', "").strip()
            chromStart, chromEnd = lines[(i + 1) * 6 + 1].replace("\n", "").split()[-1].split("..")
            anchors_list.append([chr_name, frag, int(chromStart), int(chromEnd)])
        except (IndexError, ValueError):
            continue

    return anchors_list


def get_all_anchors_from_dir(crawl_dir="crawl_data"):
    all_anchors_list = []
    if os.path.exists(crawl_dir) and os.path.isdir(crawl_dir):
        for file in os.listdir(crawl_dir):
            all_anchors_list.extend(get_anchors_from_file(os.path.join(crawl_dir, file)))
    return all_anchors_list


_PRIMARY_GOLD_CACHE = {}

def get_primary_gold_frags(primary_chr):
    if primary_chr in _PRIMARY_GOLD_CACHE:
        return _PRIMARY_GOLD_CACHE[primary_chr]
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=gold;chrom={primary_chr}"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            data = resp.json().get("gold", [])
            frags = {item["frag"] for item in data if "frag" in item}
            _PRIMARY_GOLD_CACHE[primary_chr] = frags
            return frags
    except requests.RequestException:
        pass
    _PRIMARY_GOLD_CACHE[primary_chr] = set()
    return _PRIMARY_GOLD_CACHE[primary_chr]


def determine_anchors_dynamic(df):
    is_anchor = []
    print("Checking anchor status via UCSC assembly tracks...")
    for seq in tqdm(df.itertuples(), total=len(df)):
        primary_chr = seq.name.split("_")[0]
        primary_frags = get_primary_gold_frags(primary_chr)
        is_anchor.append(1 if seq.frag in primary_frags else 0)
    return is_anchor


def main():
    parser = argparse.ArgumentParser(description="Identify anchored fragments in ALT contigs.")
    parser.add_argument("-i", "--input", default="sequence.csv", help="Input sequence.csv")
    parser.add_argument("-o", "--output", default="sequence.isAnchor.csv", help="Output file")
    parser.add_argument("-c", "--crawl-dir", default="crawl_data", help="Directory of crawled NCBI files")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Error: {args.input} not found")

    df = read_sequences(args.input)
    all_anchors = get_all_anchors_from_dir(args.crawl_dir)

    if len(all_anchors) > 0:
        anchor_df = pd.DataFrame(all_anchors, columns=["chrom", "frag", "chromStart", "chromEnd"])
        is_anchor = []
        for seq in df.itertuples():
            df_check = anchor_df[
                (anchor_df['chrom'] == seq.name) &
                (anchor_df["frag"] == seq.frag) &
                (anchor_df['chromStart'] == seq.chromStart + 1) &
                (anchor_df['chromEnd'] == seq.chromEnd)
            ]
            is_anchor.append(1 if len(df_check) > 0 else 0)
    else:
        is_anchor = determine_anchors_dynamic(df)

    df["anchor"] = is_anchor
    df.to_csv(args.output, index=False, sep="\t")
    print(f"Saved {args.output} (anchored: {sum(is_anchor)}/{len(df)})")


if __name__ == "__main__":
    main()