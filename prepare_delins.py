#!/usr/bin/env python3
import sys
import argparse
import gzip
from typing import Tuple


def trim_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Trim common suffix and prefix bases from REF and ALT."""
    ref = ref.upper()
    alt = alt.upper()

    # Trim common suffix
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # Trim common prefix
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    return pos, ref, alt


def open_vcf(path: str, mode: str = "rt"):
    if path == "-" or path is None:
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def process_vcf(input_vcf: str, output_vcf: str, resolve_overlaps: bool = True):
    in_f = open_vcf(input_vcf, "rt")
    out_f = open_vcf(output_vcf, "wt")

    header_lines = []
    records = []

    for line in in_f:
        if line.startswith("#"):
            header_lines.append(line)
            continue

        parts = line.strip().split("\t")
        if len(parts) < 8:
            continue

        chrom = parts[0]
        pos = int(parts[1])
        var_id = parts[2]
        ref = parts[3]
        alt = parts[4]
        qual = parts[5]
        filt = parts[6]
        info = parts[7]

        new_pos, new_ref, new_alt = trim_alleles(pos, ref, alt)
        end_pos = new_pos + len(new_ref) - 1

        records.append({
            "chrom": chrom,
            "pos": new_pos,
            "end": end_pos,
            "id": var_id,
            "ref": new_ref,
            "alt": new_alt,
            "qual": qual,
            "filt": filt,
            "info": info
        })

    if in_f != sys.stdin:
        in_f.close()

    # Sort records
    records.sort(key=lambda r: (r["chrom"], r["pos"], r["end"]))

    for h in header_lines:
        out_f.write(h)

    written_count = 0
    skipped_count = 0
    last_chrom = None
    last_end = -1

    for r in records:
        if resolve_overlaps:
            if r["chrom"] == last_chrom and r["pos"] <= last_end:
                skipped_count += 1
                continue

        out_f.write(f"{r['chrom']}\t{r['pos']}\t{r['id']}\t{r['ref']}\t{r['alt']}\t{r['qual']}\t{r['filt']}\t{r['info']}\n")
        written_count += 1
        last_chrom = r["chrom"]
        last_end = r["end"]

    if out_f != sys.stdout:
        out_f.close()

    sys.stderr.write(f"Processed: {written_count} variants retained, {skipped_count} overlapping variants removed.\n")


def main():
    parser = argparse.ArgumentParser(description="Standardize indels/delins and remove overlapping variants for reference consensus.")
    parser.add_argument("-i", "--input", default="-", help="Input VCF (default: stdin)")
    parser.add_argument("-o", "--output", default="-", help="Output VCF (default: stdout)")
    parser.add_argument("--keep-overlaps", action="store_true", help="Do not drop overlapping variants")

    args = parser.parse_args()
    process_vcf(args.input, args.output, resolve_overlaps=not args.keep_overlaps)


if __name__ == "__main__":
    main()
