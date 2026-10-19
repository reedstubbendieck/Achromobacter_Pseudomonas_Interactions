#!/usr/bin/env python3

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]


def run_sample(fna: Path, out_dir: Path, log_dir: Path, conda_env: str, cpus: int) -> str:
    sample = fna.stem
    sample_out = out_dir / sample

    if (sample_out / "index.html").exists():
        return f"[SKIP] {sample} - already complete"

    log_file = log_dir / f"{sample}.log"
    cmd = [
        "conda", "run", "-n", conda_env, "antismash",
        "--taxon", "bacteria",
        "--genefinding-tool", "prodigal",
        "--cpus", str(cpus),
        "--output-dir", str(sample_out),
        "--output-basename", sample,
        "--cb-general",
        "--cb-knownclusters",
        "--cb-subclusters",
        str(fna),
    ]
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode == 0:
        return f"[DONE] {sample}"
    return f"[FAIL] {sample} - check {log_file}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome-dir", type=Path, default=PROJ / "rawData/genomes",
                   help="directory of *.fna genomes (default: rawData/genomes)")
    p.add_argument("--output-dir", type=Path, default=PROJ / "derivedData/antismash",
                   help="antiSMASH output root (default: derivedData/antismash)")
    p.add_argument("--conda-env", default="antismash",
                   help="conda environment with antiSMASH installed (default: antismash)")
    p.add_argument("--cpus", type=int, default=4, help="CPUs per antiSMASH job (default: 4)")
    p.add_argument("--parallel", type=int, default=4,
                   help="number of genomes to process at once (default: 4)")
    args = p.parse_args()

    log_dir = args.output_dir.parent / "antismash_logs"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    genomes = sorted(args.genome_dir.glob("*.fna"))
    if not genomes:
        raise SystemExit(f"No *.fna genomes found in {args.genome_dir}")
    print(f"Running antiSMASH on {len(genomes)} genomes "
          f"({args.parallel} at a time, {args.cpus} CPUs each)...")

    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        for msg in pool.map(lambda g: run_sample(g, args.output_dir, log_dir,
                                                 args.conda_env, args.cpus), genomes):
            print(msg)

    print(f"antiSMASH complete. Results in: {args.output_dir}")


if __name__ == "__main__":
    main()
