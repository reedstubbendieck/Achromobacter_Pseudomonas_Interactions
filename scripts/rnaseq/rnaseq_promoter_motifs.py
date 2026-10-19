#!/usr/bin/env python3

import csv
import os
import sys

PROJ = os.environ.get(
    "PROJ", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
WORK = os.environ.get(
    "RNASEQ_WORK",
    "/mnt/d/Data_Backup/20261005_Pa_RNA-seq/derivedData/rnaseq")
OUTDIR = os.path.join(WORK, "motifs")

UPSTREAM = 300     # bp upstream of the start codon
DOWNSTREAM = 50    # bp into the gene, so the +1 region is included
MIN_GAP = 100      # intergenic gap required to call a gene an operon leader
ALPHA = 0.05

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_genome(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
    return "".join(chunks)


def main():
    genome = load_genome(os.path.join(PROJ, "references", "PAO1",
                                      "Pseudomonas_aeruginosa_PAO1_107.fna"))
    print(f"genome length: {len(genome)}")

    with open(os.path.join(WORK, "ref", "PAO1_gene_annotation.tsv")) as fh:
        ann = list(csv.DictReader(fh, delimiter="\t"))
    for a in ann:
        a["start"] = int(a["start"])
        a["end"] = int(a["end"])

    with open(os.path.join(PROJ, "derivedData", "rnaseq",
                           "deseq2_results_coculture_vs_monoculture.tsv")) as fh:
        res = {r["target_id"]: r for r in csv.DictReader(fh, delimiter="\t")
               if r["padj"] not in ("", "NA")}

    by_pos = sorted(ann, key=lambda a: a["start"])

    # nearest upstream neighbour on the same strand, in transcription order
    leaders = set()
    for i, a in enumerate(by_pos):
        if a["strand"] == "+":
            prev = [b for b in by_pos[max(0, i - 6):i] if b["end"] < a["start"]]
            gap = a["start"] - max((b["end"] for b in prev), default=0)
        else:
            nxt = [b for b in by_pos[i + 1:i + 7] if b["start"] > a["end"]]
            gap = min((b["start"] for b in nxt), default=len(genome)) - a["end"]
        if gap >= MIN_GAP:
            leaders.add(a["target_id"])
    print(f"operon leaders (gap >= {MIN_GAP} bp): {len(leaders)} of {len(ann)} genes")

    def promoter(a):
        if a["strand"] == "+":
            lo = a["start"] - 1 - UPSTREAM
            hi = a["start"] - 1 + DOWNSTREAM
            if lo < 0 or hi > len(genome):
                return None
            return genome[lo:hi]
        lo = a["end"] - DOWNSTREAM
        hi = a["end"] + UPSTREAM
        if lo < 0 or hi > len(genome):
            return None
        return revcomp(genome[lo:hi])

    os.makedirs(OUTDIR, exist_ok=True)
    sets = {"induced": [], "repressed": [], "background": []}
    for a in ann:
        tid = a["target_id"]
        if tid not in res or tid not in leaders:
            continue
        seq = promoter(a)
        if not seq or len(seq) < UPSTREAM:
            continue
        r = res[tid]
        q = float(r["padj"])
        lfc = float(r["log2FoldChange_shrunk"])
        if q < ALPHA and lfc > 0:
            key = "induced"
        elif q < ALPHA and lfc < 0:
            key = "repressed"
        else:
            key = "background"
        name = a["gene_name"] or tid
        sets[key].append((f"{tid}|{name}", seq))

    for key, items in sets.items():
        path = os.path.join(OUTDIR, f"promoters_{key}.fa")
        with open(path, "w", newline="\n") as fh:
            for h, s in items:
                fh.write(f">{h}\n")
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + "\n")
        print(f"{key:<11} {len(items):>5} promoters -> {path}")

    if not sets["induced"]:
        sys.exit("no induced leader promoters extracted")


if __name__ == "__main__":
    main()
