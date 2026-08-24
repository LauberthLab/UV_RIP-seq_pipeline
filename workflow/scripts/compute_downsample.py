#!/usr/bin/env python3
"""Compute per-sample downsampling fractions so libraries share a common depth.

Reproduces, automatically, what was done by hand in UVRIP_downsample_BW.sh
(e.g. seed 42, fractions 0.60172 / 0.18306 / 0.18339 ...). It counts mapped,
primary, deduplicated reads in each nodup BAM and picks a target depth, then
emits the fraction and the ready-to-use `samtools view -s` argument.

Usage:
  compute_downsample.py \
      --bam-dir DIR --meta metadata.tsv \
      --mode {min,target} --target-reads N \
      --seed 42 --group-by {all,type} \
      --out downsample_plan.tsv

Output columns (TSV):
  sample  bam  mapped_reads  target  fraction  s_arg  need_downsample
"""
import argparse, csv, os, subprocess, sys


def count_reads(bam):
    """Primary, mapped, non-duplicate, non-supplementary read count.
    -F 3332 excludes: unmapped(4) + secondary(256) + dup(1024) + supplementary(2048)."""
    out = subprocess.check_output(["samtools", "view", "-c", "-F", "3332", bam])
    return int(out.strip())


def load_meta(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    for r in rows:
        for k in list(r.keys()):
            r[k] = (r[k] or "").strip()
    return rows


def s_arg(seed, frac):
    """Build the samtools `-s SEED.FRACTION` argument.
    frac is kept fraction in (0,1]; formatted to 5 decimals, decimal digits
    appended after the integer seed (e.g. seed=42, frac=0.18306 -> '42.18306')."""
    if frac >= 1.0:
        return None  # no downsampling
    digits = f"{frac:.5f}".split(".")[1]
    return f"{int(seed)}.{digits}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-dir", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--mode", choices=["min", "target"], default="min")
    ap.add_argument("--target-reads", type=int, default=20_000_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--group-by", choices=["all", "type"], default="all")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = load_meta(args.meta)

    # Count reads per sample.
    for r in rows:
        bam = os.path.join(args.bam_dir, f"nodup_{r['sample']}.bam")
        if not os.path.exists(bam):
            sys.stderr.write(f"ERROR: missing BAM {bam} — run the process step first.\n")
            sys.exit(1)
        r["_bam"] = bam
        r["_reads"] = count_reads(bam)

    # Determine the target depth per grouping.
    def group_key(r):
        return r["type"] if args.group_by == "type" else "ALL"

    targets = {}
    if args.mode == "min":
        for r in rows:
            k = group_key(r)
            targets[k] = min(targets.get(k, r["_reads"]), r["_reads"])
    else:  # target
        for r in rows:
            targets[group_key(r)] = args.target_reads

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "bam", "mapped_reads", "target",
                    "fraction", "s_arg", "need_downsample"])
        for r in rows:
            tgt = targets[group_key(r)]
            # Manual override wins.
            override = r.get("downsample_frac", "")
            if override:
                frac = float(override)
            else:
                frac = min(1.0, tgt / r["_reads"]) if r["_reads"] else 1.0
            arg = s_arg(args.seed, frac)
            need = "true" if arg is not None else "false"
            w.writerow([r["sample"], r["_bam"], r["_reads"], tgt,
                        f"{frac:.5f}", arg if arg else "NA", need])

    sys.stderr.write(f"Wrote downsample plan: {args.out}\n")


if __name__ == "__main__":
    main()
