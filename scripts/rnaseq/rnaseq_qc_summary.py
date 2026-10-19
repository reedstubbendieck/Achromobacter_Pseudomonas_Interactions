#!/usr/bin/env python3

import csv
import json
import os
import sys
from pathlib import Path

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
WORK = Path(os.environ.get(
    "RNASEQ_WORK",
    "/mnt/d/Data_Backup/20261005_Pa_RNA-seq/derivedData/rnaseq"))
OUT = PROJ / "derivedData" / "rnaseq" / "qc_summary.tsv"
OUT_CLASS = PROJ / "derivedData" / "rnaseq" / "qc_feature_class_fractions.tsv"

COLS = [
    "sample", "condition", "replicate", "library_id",
    "raw_read_pairs", "trimmed_read_pairs", "pct_reads_kept",
    "pct_q30_raw", "pct_q30_trimmed", "mean_len_r1_trimmed",
    "pct_adapter_trimmed", "pct_duplication",
    "kallisto_processed", "kallisto_pseudoaligned", "pct_pseudoaligned",
    "kallisto_unique", "pct_unique", "frag_length_mean",
]


def frag_length(name):
    path = WORK / "logs" / f"{name}.kallisto.log"
    if not path.exists():
        return None
    for line in open(path):
        if "estimated average fragment length" in line:
            return float(line.rsplit(":", 1)[1].strip())
    return None


def feature_map():
    path = WORK / "ref" / "PAO1_gene_annotation.tsv"
    with open(path) as fh:
        return {r["target_id"]: r["feature"] for r in csv.DictReader(fh, delimiter="\t")}


def class_fractions(name, feat):
    path = WORK / "kallisto" / name / "abundance.tsv"
    if not path.exists():
        return None
    totals, total = {}, 0.0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cls = feat.get(r["target_id"], "unknown")
            c = float(r["est_counts"])
            totals[cls] = totals.get(cls, 0.0) + c
            total += c
    if total <= 0:
        return None
    return {k: 100.0 * v / total for k, v in totals.items()}, total


def main():
    sample_sheet = PROJ / "rawData" / "rnaseq" / "samples.tsv"
    with open(sample_sheet) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    rows, problems = [], []
    for s in samples:
        name = s["sample"]
        row = {
            "sample": name,
            "condition": s["condition"],
            "replicate": s["replicate"],
            "library_id": s["library_id"],
        }

        fp = WORK / "fastp" / f"{name}.fastp.json"
        if fp.exists():
            d = json.load(open(fp))
            before = d["summary"]["before_filtering"]
            after = d["summary"]["after_filtering"]
            row["raw_read_pairs"] = before["total_reads"] // 2
            row["trimmed_read_pairs"] = after["total_reads"] // 2
            row["pct_reads_kept"] = round(
                100.0 * after["total_reads"] / before["total_reads"], 2)
            row["pct_q30_raw"] = round(100.0 * before["q30_rate"], 2)
            row["pct_q30_trimmed"] = round(100.0 * after["q30_rate"], 2)
            row["mean_len_r1_trimmed"] = after["read1_mean_length"]
            ad = d.get("adapter_cutting", {})
            if ad and before["total_reads"]:
                row["pct_adapter_trimmed"] = round(
                    100.0 * ad.get("adapter_trimmed_reads", 0) / before["total_reads"], 2)
            row["pct_duplication"] = round(
                100.0 * d.get("duplication", {}).get("rate", float("nan")), 2)
        else:
            problems.append(f"{name}: no fastp json")

        kp = WORK / "kallisto" / name / "run_info.json"
        if kp.exists():
            d = json.load(open(kp))
            row["kallisto_processed"] = d["n_processed"]
            row["kallisto_pseudoaligned"] = d["n_pseudoaligned"]
            row["pct_pseudoaligned"] = round(
                100.0 * d["n_pseudoaligned"] / d["n_processed"], 2)
            row["kallisto_unique"] = d.get("n_unique")
            if d.get("n_unique") is not None:
                row["pct_unique"] = round(
                    100.0 * d["n_unique"] / d["n_processed"], 2)
            row["frag_length_mean"] = frag_length(name)
        else:
            problems.append(f"{name}: no kallisto run_info")

        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", restval="NA")
        w.writeheader()
        w.writerows(rows)

    feat = feature_map()
    classes = sorted(set(feat.values()))
    class_rows = []
    for s in samples:
        got = class_fractions(s["sample"], feat)
        if got is None:
            continue
        frac, total = got
        cr = {"sample": s["sample"], "condition": s["condition"],
              "assigned_counts": round(total)}
        for c in classes:
            cr["pct_" + c] = round(frac.get(c, 0.0), 3)
        class_rows.append(cr)

    if class_rows:
        cols2 = ["sample", "condition", "assigned_counts"] + ["pct_" + c for c in classes]
        with open(OUT_CLASS, "w", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=cols2, delimiter="\t", restval="NA")
            w.writeheader()
            w.writerows(class_rows)

    show = ["sample", "condition", "raw_read_pairs", "trimmed_read_pairs",
            "pct_reads_kept", "pct_q30_trimmed", "pct_pseudoaligned",
            "pct_unique", "frag_length_mean"]
    print("\t".join(show))
    for r in rows:
        print("\t".join(str(r.get(c, "NA")) for c in show))
    print(f"\nwritten -> {OUT}")

    if class_rows:
        print("\nshare of assigned reads by feature class (%):")
        show2 = ["sample", "condition"] + ["pct_" + c for c in classes]
        print("\t".join(s.replace("pct_", "") for s in show2))
        for r in class_rows:
            print("\t".join(str(r.get(c, "NA")) for c in show2))
        print(f"\nwritten -> {OUT_CLASS}")
    if problems:
        print("\nincomplete:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
