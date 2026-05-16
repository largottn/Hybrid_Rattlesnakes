#!/usr/bin/env python3
"""
Parse a SnpEff-annotated VCF into a tidy summary TSV.

For each variant, extracts the highest-impact annotation from the ANN
INFO field plus any carry-through INFO fields (CHILD_AB, CHILD_DP, etc.)
added by dnm_tsv_to_vcf.py.

ANN field format (per SnpEff docs):
    Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID |
    Feature_Type | Feature_ID | Transcript_BioType | Rank |
    HGVS.c | HGVS.p | ...

SnpEff orders multiple annotations per variant by impact (HIGH > MODERATE
> LOW > MODIFIER), so the first comma-separated entry is the "top" hit.

Usage:
    python3 parse_snpeff_results.py \\
        --input  trio_child1_denovo_annotated.vcf \\
        --output trio_child1_denovo_summary.tsv
"""

import argparse
import sys
from collections import Counter


OUTPUT_COLS = [
    "CHROM", "POS", "REF", "ALT",
    "EFFECT", "IMPACT", "GENE", "GENE_ID", "BIOTYPE",
    "HGVS_C", "HGVS_P", "CHILD_AB", "CHILD_DP",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",  required=True, help="SnpEff-annotated VCF")
    p.add_argument("--output", required=True, help="Output summary TSV")
    return p.parse_args()


def parse_info(info_str):
    """Split a VCF INFO field into a dict. Flag-style entries map to ''."""
    out = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        elif item:
            out[item] = ""
    return out


def top_annotation(ann_str):
    """
    Return the first (highest-impact) annotation from an ANN string.

    Returns a dict with effect/impact/gene/gene_id/biotype/hgvs_c/hgvs_p,
    or all dots if the ANN field is missing or malformed.
    """
    blank = {k: "." for k in ("effect", "impact", "gene", "gene_id",
                              "biotype", "hgvs_c", "hgvs_p")}
    if not ann_str:
        return blank

    first = ann_str.split(",")[0].split("|")
    # Indices per SnpEff ANN spec
    def get(i):
        return first[i] if i < len(first) and first[i] else "."

    return {
        "effect":  get(1),
        "impact":  get(2),
        "gene":    get(3),
        "gene_id": get(4),
        "biotype": get(7),
        "hgvs_c":  get(9),
        "hgvs_p":  get(10),
    }


def main():
    args = parse_args()

    variants = []
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            info = parse_info(f[7])
            ann = top_annotation(info.get("ANN", ""))

            variants.append({
                "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt,
                "EFFECT":   ann["effect"],
                "IMPACT":   ann["impact"],
                "GENE":     ann["gene"],
                "GENE_ID":  ann["gene_id"],
                "BIOTYPE":  ann["biotype"],
                "HGVS_C":   ann["hgvs_c"],
                "HGVS_P":   ann["hgvs_p"],
                "CHILD_AB": info.get("CHILD_AB", "."),
                "CHILD_DP": info.get("CHILD_DP", "."),
            })

    with open(args.output, "w") as out:
        out.write("\t".join(OUTPUT_COLS) + "\n")
        for v in variants:
            out.write("\t".join(v[c] for c in OUTPUT_COLS) + "\n")

    # ---- Summary to stderr ----
    bar = "=" * 60
    print(f"\n{bar}", file=sys.stderr)
    print(f"SnpEff annotation summary ({len(variants)} variant(s))", file=sys.stderr)
    print(bar, file=sys.stderr)
    for v in variants:
        print(f"  {v['CHROM']}:{v['POS']} {v['REF']}>{v['ALT']}", file=sys.stderr)
        print(f"    Effect:   {v['EFFECT']} ({v['IMPACT']})", file=sys.stderr)
        print(f"    Gene:     {v['GENE']} ({v['GENE_ID']})", file=sys.stderr)
        if v["HGVS_P"] != ".":
            print(f"    Protein:  {v['HGVS_P']}", file=sys.stderr)
        print("", file=sys.stderr)

    if variants:
        impacts = Counter(v["IMPACT"] for v in variants)
        print("Impact distribution:", file=sys.stderr)
        for impact, count in sorted(impacts.items()):
            print(f"  {impact}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
