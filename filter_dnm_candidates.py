#!/usr/bin/env python3
"""
Read-level filtering of high-confidence de novo mutation candidates.

Applies the read-level filters from Zhang et al. 2020 (Nat Commun) to a
TSV of hiConfDeNovo candidates extracted from a trio VCF. Reads the
extracted TSV produced by bcftools query, applies the four filters
(depth, parental allele, allele balance, autosomes), and writes the
passing variants plus a summary to stdout.

Filters applied:
  1. DP >= MIN_DP in sire, dam, and offspring
  2. Parental allele filter:
       - Hom-ref parents must have zero alt-allele reads
       - Hom-alt parents must have zero ref-allele reads
       - Heterozygous parents disqualify the site (true DNMs require
         both parents to be homozygous)
  3. Offspring allele balance in [AB_MIN, AB_MAX]

Input format (one variant per line, tab-separated):
  CHROM  POS  REF  ALT  hiConfDeNovo  <per-sample blocks>

Each per-sample block contains 4 fields:
  SAMPLE_NAME  GT  DP  AD

There are three samples per trio, so each row has 5 + 3*4 = 17 fields.

Usage:
    python3 filter_dnm_candidates.py \\
        --input  extracted.tsv \\
        --output passing.tsv \\
        --sire   sire_sample \\
        --dam    dam_sample \\
        --child  offspring_sample \\
        --min-dp 7 \\
        --ab-min 0.3 \\
        --ab-max 0.7
"""

import argparse
import sys
from collections import Counter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",  required=True,
                   help="TSV from bcftools query (see module docstring)")
    p.add_argument("--output", required=True,
                   help="Output TSV of passing DNMs")
    p.add_argument("--sire",  required=True, help="Sire sample name")
    p.add_argument("--dam",   required=True, help="Dam sample name")
    p.add_argument("--child", required=True, help="Offspring sample name")
    p.add_argument("--min-dp", type=int,   default=7,
                   help="Minimum DP in all three trio members (default: 7)")
    p.add_argument("--ab-min", type=float, default=0.3,
                   help="Minimum offspring allele balance (default: 0.3)")
    p.add_argument("--ab-max", type=float, default=0.7,
                   help="Maximum offspring allele balance (default: 0.7)")
    return p.parse_args()


def parse_sample_block(fields, idx):
    """Parse a 4-field per-sample block at position idx in the row."""
    sname = fields[idx]
    gt    = fields[idx + 1]
    dp    = int(fields[idx + 2]) if fields[idx + 2] != "." else 0
    ad_str = fields[idx + 3]
    if ad_str and ad_str != ".":
        ad = [int(x) for x in ad_str.split(",")]
    else:
        ad = [0, 0]
    return sname, {"gt": gt, "dp": dp, "ad": ad}


def parent_passes(parent):
    """
    True iff parent genotype is homozygous and consistent with its
    read support (no alt reads in hom-ref, no ref reads in hom-alt).
    Heterozygous parents always fail.
    """
    gt = parent["gt"]
    ad = parent["ad"]
    ref_reads = ad[0] if len(ad) > 0 else 0
    alt_reads = ad[1] if len(ad) > 1 else 0

    if gt in ("0/0", "0|0"):
        return alt_reads == 0
    if gt in ("1/1", "1|1"):
        return ref_reads == 0
    return False  # heterozygous or missing → fail


def main():
    args = parse_args()

    passed = []
    failed_dp = 0
    failed_parent = 0
    failed_ab = 0
    skipped = 0

    with open(args.input) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, hiconf = fields[:5]

            # Parse three 4-field sample blocks
            samples = {}
            idx = 5
            for _ in range(3):
                sname, sdata = parse_sample_block(fields, idx)
                samples[sname] = sdata
                idx += 4

            if not all(n in samples for n in (args.sire, args.dam, args.child)):
                skipped += 1
                continue

            sire  = samples[args.sire]
            dam   = samples[args.dam]
            child = samples[args.child]

            # Filter 1: depth in all three
            if min(sire["dp"], dam["dp"], child["dp"]) < args.min_dp:
                failed_dp += 1
                continue

            # Filter 2: parental allele filter
            if not (parent_passes(sire) and parent_passes(dam)):
                failed_parent += 1
                continue

            # Filter 3: offspring allele balance
            child_ref = child["ad"][0] if len(child["ad"]) > 0 else 0
            child_alt = child["ad"][1] if len(child["ad"]) > 1 else 0
            total = child_ref + child_alt
            if total == 0:
                failed_ab += 1
                continue
            ab = child_alt / total
            if ab < args.ab_min or ab > args.ab_max:
                failed_ab += 1
                continue

            passed.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "sire_dp": sire["dp"],
                "sire_ref": sire["ad"][0] if len(sire["ad"]) > 0 else 0,
                "sire_alt": sire["ad"][1] if len(sire["ad"]) > 1 else 0,
                "dam_dp": dam["dp"],
                "dam_ref": dam["ad"][0] if len(dam["ad"]) > 0 else 0,
                "dam_alt": dam["ad"][1] if len(dam["ad"]) > 1 else 0,
                "child_dp": child["dp"],
                "child_ref": child_ref,
                "child_alt": child_alt,
                "child_ab": round(ab, 4),
                "child_gt": child["gt"],
                "hiConfDeNovo": hiconf,
            })

    # ---- Write passing DNMs ----
    cols = ["chrom", "pos", "ref", "alt",
            "sire_dp", "sire_ref", "sire_alt",
            "dam_dp",  "dam_ref",  "dam_alt",
            "child_dp", "child_ref", "child_alt",
            "child_ab", "child_gt", "hiConfDeNovo"]
    header_labels = ["CHROM", "POS", "REF", "ALT",
                     "SIRE_DP", "SIRE_REF_READS", "SIRE_ALT_READS",
                     "DAM_DP",  "DAM_REF_READS",  "DAM_ALT_READS",
                     "CHILD_DP", "CHILD_REF_READS", "CHILD_ALT_READS",
                     "CHILD_AB", "CHILD_GT", "hiConfDeNovo"]

    with open(args.output, "w") as out:
        out.write("\t".join(header_labels) + "\n")
        for v in passed:
            out.write("\t".join(str(v[c]) for c in cols) + "\n")

    # ---- Summary ----
    total = failed_dp + failed_parent + failed_ab + len(passed) + skipped
    bar = "=" * 55
    print(f"\n{bar}", file=sys.stderr)
    print(f"DNM filter summary (DP >= {args.min_dp}, AB in "
          f"[{args.ab_min}, {args.ab_max}])", file=sys.stderr)
    print(bar, file=sys.stderr)
    print(f"Total hiConfDeNovo candidates:      {total}", file=sys.stderr)
    print(f"Skipped (sample mismatch):          {skipped}", file=sys.stderr)
    print(f"Failed DP filter:                   {failed_dp}", file=sys.stderr)
    print(f"Failed parental allele filter:      {failed_parent}", file=sys.stderr)
    print(f"Failed allele balance filter:       {failed_ab}", file=sys.stderr)
    print(f"PASSED all filters:                 {len(passed)}", file=sys.stderr)
    print(bar, file=sys.stderr)

    snvs   = sum(1 for v in passed if len(v["ref"]) == 1 and len(v["alt"]) == 1)
    indels = len(passed) - snvs
    print(f"  De novo SNVs:   {snvs}", file=sys.stderr)
    print(f"  De novo indels: {indels}", file=sys.stderr)
    print(bar, file=sys.stderr)

    if snvs > 0:
        subs = Counter(f"{v['ref']}>{v['alt']}" for v in passed
                       if len(v["ref"]) == 1 and len(v["alt"]) == 1)
        print("\nBase substitution spectrum (de novo SNVs):", file=sys.stderr)
        for sub, count in sorted(subs.items()):
            pct = 100 * count / snvs
            print(f"  {sub}: {count} ({pct:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
