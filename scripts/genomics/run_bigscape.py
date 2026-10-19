#!/usr/bin/env python3

import argparse
import os
import subprocess
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pfam", type=Path, required=True,
                   help="path to Pfam-A.hmm (ships with antiSMASH databases)")
    p.add_argument("--input-dir", type=Path, default=PROJ / "derivedData/antismash",
                   help="antiSMASH output root (default: derivedData/antismash)")
    p.add_argument("--output-dir", type=Path, default=PROJ / "derivedData/bigscape",
                   help="BiG-SCAPE output dir (default: derivedData/bigscape)")
    p.add_argument("--conda-env", default="bigscape",
                   help="conda environment with BiG-SCAPE installed (default: bigscape)")
    p.add_argument("--mibig-version", default="4.0",
                   help="MIBiG version to include (default: 4.0)")
    p.add_argument("--cores", type=int, default=os.cpu_count(),
                   help="CPU cores (default: all available)")
    p.add_argument("--label", default="achromobacter", help="run label (default: achromobacter)")
    args = p.parse_args()

    if not args.input_dir.is_dir():
        raise SystemExit(f"antiSMASH output not found: {args.input_dir}\nRun run_antismash.py first.")
    if not args.pfam.is_file():
        raise SystemExit(f"Pfam HMM not found: {args.pfam}")

    gbks = list(args.input_dir.rglob("*.region*.gbk"))
    print(f"Found {len(gbks)} region .gbk files across all samples")
    if not gbks:
        raise SystemExit("No region .gbk files found. Check antiSMASH completed successfully.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Running BiG-SCAPE cluster with MIBiG {args.mibig_version} "
          f"on {len(gbks)} BGCs using {args.cores} cores...")

    cmd = [
        "conda", "run", "--no-capture-output", "-n", args.conda_env, "bigscape", "cluster",
        "--input-dir", str(args.input_dir),
        "--input-mode", "recursive",
        "--output-dir", str(args.output_dir),
        "--label", args.label,
        "--pfam-path", str(args.pfam),
        "--mibig-version", args.mibig_version,
        "--cores", str(args.cores),
        "--classify", "category",
        "--mix",
        "--include-singletons",
    ]
    subprocess.run(cmd, check=True)

    print(f"BiG-SCAPE complete. Results in: {args.output_dir}")
    print(f"Open the visualisation: {args.output_dir / 'index.html'}")


if __name__ == "__main__":
    main()
