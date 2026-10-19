#!/usr/bin/env python3

import gzip
import re
import sys
from pathlib import Path

FFN = Path(sys.argv[1])
OUT_FASTA = Path(sys.argv[2])
OUT_ANNOT = Path(sys.argv[3])

COORD_RE = re.compile(r"^(?P<replicon>[^:]+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)$")


def opener(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def parse_header(line):
    """>PGD134012 cds chromosome:483-2027(+) name=dnaA;locus_tag=PA0001;..."""
    head = line[1:].rstrip("\n")
    parts = head.split(" ", 3)
    pgd, feature, coords = parts[0], parts[1], parts[2]
    tail = parts[3] if len(parts) > 3 else ""

    fields = {}
    for kv in tail.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            fields[k.strip()] = v.strip()

    m = COORD_RE.match(coords)
    rec = {
        "pgd_id": pgd,
        "feature": feature,
        "locus_tag": fields.get("locus_tag", ""),
        "gene_name": fields.get("name", ""),
        "product": fields.get("product", ""),
        "replicon": fields.get("replicon_accession", ""),
        "start": m.group("start") if m else "",
        "end": m.group("end") if m else "",
        "strand": m.group("strand") if m else "",
    }
    return rec


def main():
    records = []           # (rec, [seq lines]) in file order
    rec, seq = None, []
    with opener(FFN) as fh:
        for line in fh:
            if line.startswith(">"):
                if rec is not None:
                    records.append((rec, "".join(seq)))
                rec, seq = parse_header(line), []
            else:
                seq.append(line.strip())
    if rec is not None:
        records.append((rec, "".join(seq)))

    seen = {}
    kept, dropped, conflicts = [], 0, []
    for rec, sequence in records:
        key = rec["locus_tag"] or rec["pgd_id"]
        if not key:
            raise SystemExit(f"record with no usable ID: {rec}")
        if key in seen:
            dropped += 1
            if seen[key] != sequence:
                conflicts.append(key)
            continue
        seen[key] = sequence
        rec["target_id"] = key
        rec["length"] = len(sequence)
        kept.append((rec, sequence))

    if conflicts:
        raise SystemExit(
            "duplicate IDs carry DIFFERENT sequences, refusing to pick one: "
            + ", ".join(sorted(conflicts))
        )

    with open(OUT_FASTA, "w", newline="\n") as fh:
        for rec, sequence in kept:
            fh.write(f">{rec['target_id']}\n")
            for i in range(0, len(sequence), 60):
                fh.write(sequence[i:i + 60] + "\n")

    cols = ["target_id", "locus_tag", "pgd_id", "gene_name", "product",
            "feature", "replicon", "start", "end", "strand", "length"]
    with open(OUT_ANNOT, "w", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for rec, _ in kept:
            fh.write("\t".join(str(rec[c]) for c in cols) + "\n")

    by_feature = {}
    for rec, _ in kept:
        by_feature[rec["feature"]] = by_feature.get(rec["feature"], 0) + 1

    print(f"input records      : {len(records)}")
    print(f"dropped duplicates : {dropped} (identical sequence, first kept)")
    print(f"written targets    : {len(kept)}")
    for feat in sorted(by_feature):
        print(f"  {feat:<8} {by_feature[feat]}")
    print(f"fasta  -> {OUT_FASTA}")
    print(f"annot  -> {OUT_ANNOT}")


if __name__ == "__main__":
    main()
