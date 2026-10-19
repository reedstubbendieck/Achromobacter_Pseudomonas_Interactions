#!/usr/bin/env python3

import argparse
import csv
import shutil
import subprocess
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]


def newest_mix_clustering(bigscape_dir: Path) -> Path:
    hits = sorted(bigscape_dir.glob("output_files/*_c0.3/mix/mix_clustering_c0.3.tsv"))
    if not hits:
        raise SystemExit(f"No mix clustering table under {bigscape_dir}/output_files/*_c0.3/mix/")
    return hits[-1]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--family", default="FAM_00001",
                   help="BiG-SCAPE family to plot (default: FAM_00001 = alcaligin)")
    p.add_argument("--antismash-dir", type=Path, default=PROJ / "derivedData/antismash",
                   help="antiSMASH output root (default: derivedData/antismash)")
    p.add_argument("--bigscape-dir", type=Path, default=PROJ / "derivedData/bigscape",
                   help="BiG-SCAPE output dir (default: derivedData/bigscape)")
    p.add_argument("--reference", type=Path,
                   default=PROJ / "references/BGCs/BGC0002464.gbk",
                   help="MIBiG reference .gbk to include (default: BGC0002464, alcaligin)")
    p.add_argument("--output-dir", type=Path, default=PROJ / "derivedData/clinker",
                   help="clinker output dir (default: derivedData/clinker)")
    p.add_argument("--conda-env", default="clinker",
                   help="conda environment with clinker installed (default: clinker)")
    args = p.parse_args()

    clust = newest_mix_clustering(args.bigscape_dir)
    print(f"Using clustering: {clust}")

    # Family members that are genome regions (skip MIBiG 'BGC*' entries)
    regions = []
    with open(clust, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["Family"] == args.family and not row["GBK"].startswith("BGC"):
                regions.append(row["GBK"])
    print(f"{args.family} genome regions: {len(regions)}")

    input_dir = args.output_dir / "input_BGCs"
    if input_dir.exists():
        shutil.rmtree(input_dir)
    input_dir.mkdir(parents=True)

    for region in regions:
        matches = list(args.antismash_dir.rglob(f"{region}.gbk"))
        if not matches:
            raise SystemExit(f"Could not find {region}.gbk under {args.antismash_dir}")
        shutil.copy(matches[0], input_dir)

    if not args.reference.is_file():
        raise SystemExit(f"Reference BGC not found: {args.reference}")
    shutil.copy(args.reference, input_dir)

    gbks = sorted(input_dir.glob("*.gbk"))
    print(f"Input BGCs ({len(gbks)} files) in {input_dir}")

    html_out = args.output_dir / f"{args.family}_clinker.html"
    aln_out = args.output_dir / f"{args.family}_alignments.txt"
    cmd = [
        "conda", "run", "--no-capture-output", "-n", args.conda_env,
        "clinker", *[str(g) for g in gbks],
        "-p", str(html_out),
        "-o", str(aln_out),
        "-f",
    ]
    subprocess.run(cmd, check=True)

    print(f"Done. Plot: {html_out}")


if __name__ == "__main__":
    main()
