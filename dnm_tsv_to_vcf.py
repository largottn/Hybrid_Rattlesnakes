#!/usr/bin/env python3
"""
Convert a de novo mutation TSV (output of the read-level filter, step 5c)
to a minimal VCF suitable as input to SnpEff.

The input TSV is expected to have a header row whose columns include at
least CHROM, POS, REF, ALT. Optional per-trio depth/AB columns are
copied into the VCF INFO field if present, so they survive annotation.

Optional INFO carry-throughs (added when columns are present):
    CHILD_AB, CHILD_DP, SIRE_DP, DAM_DP

Usage:
    python3 dnm_tsv_to_vcf.py \\
        --input  trio_child1_final_denovo.tsv \\
        --output trio_child1_for_snpeff.vcf
"""

import argparse
import sys

REQUIRED_COLS = ("CHROM", "POS", "REF", "ALT")
INFO_COLS = ("CHILD_AB", "CHILD_DP", "SIRE_DP", "DAM_DP")

INFO_HEADER_LINES = {
    "CHILD_AB": '##INFO=<ID=CHILD_AB,Number=1,Type=Float,Description="Allele balance in offspring">',
    "CHILD_DP": '##INFO=<ID=CHILD_DP,Number=1,Type=Integer,Description="Read depth in offspring">',
    "SIRE_DP":  '##INFO=<ID=SIRE_DP,Number=1,Type=Integer,Description="Read depth in sire">',
    "DAM_DP":   '##INFO=<ID=DAM_DP,Number=1,Type=Integer,Description="Read depth in dam">',
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",  required=True, help="DNM TSV (with header)")
    p.add_argument("--output", required=True, help="Output VCF (uncompressed)")
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.input) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    missing = [c for c in REQUIRED_COLS if c not in header]
    if missing:
        sys.exit(f"ERROR: required column(s) missing from {args.input}: {missing}")

    col_idx = {c: header.index(c) for c in header}
    info_cols_present = [c for c in INFO_COLS if c in col_idx]

    with open(args.output, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=DeNovoMutationPipeline\n")
        for c in info_cols_present:
            out.write(INFO_HEADER_LINES[c] + "\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for r in rows:
            chrom = r[col_idx["CHROM"]]
            pos   = r[col_idx["POS"]]
            ref   = r[col_idx["REF"]]
            alt   = r[col_idx["ALT"]]
            info  = ";".join(f"{c}={r[col_idx[c]]}" for c in info_cols_present) or "."
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")

    print(f"Wrote {len(rows)} variant(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
