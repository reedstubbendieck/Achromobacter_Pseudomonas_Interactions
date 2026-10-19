#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import parasail

PROJ = Path(__file__).resolve().parents[2]

# Column order: precursor (ornithine -> putrescine), then alcaligin biosynthesis,
# regulation, and transport.
GENE_ORDER = ["ArgA", "ArgB", "ArgC", "ArgD1", "ArgD2", "ArgE", "Odc",
              "AlcE", "AlcA", "AlcB", "AlcC", "AlcR", "AlcS", "FauA", "AlcD"]


def read_fasta(path):
    name, seq, out = None, [], {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        out[name] = "".join(seq)
    return out


def parse_proteome(gbk_path):
    """Return list of (locus_tag, protein_sequence) from a GenBank file."""
    proteins, locus, parts, in_trans = [], None, None, False
    with open(gbk_path) as fh:
        for line in fh:
            s = line.strip()
            if in_trans:
                if s.endswith('"'):
                    parts.append(s[:-1])
                    proteins.append((locus, "".join(parts)))
                    in_trans, parts = False, None
                else:
                    parts.append(s)
            elif s.startswith('/locus_tag="'):
                locus = s[len('/locus_tag="'):].rstrip('"')
            elif s.startswith('/translation="'):
                rest = s[len('/translation="'):]
                if rest.endswith('"'):
                    proteins.append((locus, rest[:-1]))
                else:
                    in_trans, parts = True, [rest]
    return proteins


def genomes_from_matrix(matrix_tsv):
    with open(matrix_tsv, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return [r[0] for r in rows[1:]]


def best_hit(query, proteome, matrix, gap_open, gap_extend):
    """Best Smith-Waterman hit: (locus, score, pident, qcov, target_len)."""
    qlen = len(query)
    min_len = 0.3 * qlen
    best_score, best_locus, best_target, best_tlen = -1, None, "", 0
    for locus, target in proteome:
        if len(target) < min_len:
            continue
        r = parasail.sw_striped_16(query, target, gap_open, gap_extend, matrix)
        if r.score > best_score:
            best_score, best_locus, best_target, best_tlen = r.score, locus, target, len(target)

    # Realign the winner with traceback for identity / coverage
    rt = parasail.sw_trace_striped_16(query, best_target, gap_open, gap_extend, matrix)
    tb = rt.traceback
    q, t = tb.query, tb.ref
    aln_len = len(q)
    ident = sum(1 for a, b in zip(q, t) if a == b and a != "-")
    qres = sum(1 for a in q if a != "-")
    pident = 100.0 * ident / aln_len if aln_len else 0.0
    qcov = 100.0 * qres / qlen if qlen else 0.0
    return best_locus, best_score, pident, qcov, best_tlen


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refs-dir", type=Path,
                   default=PROJ / "references/bordetella",
                   help="directory of reference *.faa files (default: references/bordetella)")
    p.add_argument("--matrix", type=Path,
                   default=PROJ / "derivedData/bigscape/absence_presence_clean.tsv")
    p.add_argument("--antismash-dir", type=Path, default=PROJ / "derivedData/antismash")
    p.add_argument("--output-dir", type=Path, default=PROJ / "derivedData/alcaligin_pathway_genes")
    p.add_argument("--min-identity", type=float, default=30.0,
                   help="percent identity to call a gene present (default: 30)")
    p.add_argument("--min-coverage", type=float, default=70.0,
                   help="percent query coverage to call a gene present (default: 70)")
    p.add_argument("--gap-open", type=int, default=11)
    p.add_argument("--gap-extend", type=int, default=1)
    args = p.parse_args()

    matrix = parasail.blosum62
    refs = {}
    for faa in sorted(args.refs_dir.glob("*.faa")):
        refs.update(read_fasta(faa))
    genes = [g for g in GENE_ORDER if g in refs] + [g for g in refs if g not in GENE_ORDER]
    genomes = genomes_from_matrix(args.matrix)
    print(f"{len(genomes)} genomes x {len(genes)} reference genes")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    hits = []
    for genome in genomes:
        gbk = args.antismash_dir / genome / f"{genome}.gbk"
        if not gbk.is_file():
            raise SystemExit(f"Missing proteome GBK: {gbk}")
        proteome = parse_proteome(gbk)
        for gene in genes:
            locus, score, pident, qcov, tlen = best_hit(
                refs[gene], proteome, matrix, args.gap_open, args.gap_extend)
            present = pident >= args.min_identity and qcov >= args.min_coverage
            hits.append(dict(genome=genome, gene=gene, locus_tag=locus,
                             target_len=tlen, score=score,
                             pident=round(pident, 1), qcov=round(qcov, 1),
                             present=int(present)))
        n_present = sum(h["present"] for h in hits if h["genome"] == genome)
        print(f"  {genome}: {n_present}/{len(genes)} genes present")

    tidy = args.output_dir / "pathway_gene_hits.tsv"
    with open(tidy, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(hits[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(hits)

    wide = args.output_dir / "pathway_gene_identity.tsv"
    by = {(h["genome"], h["gene"]): h["pident"] for h in hits}
    with open(wide, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome"] + genes)
        for genome in genomes:
            w.writerow([genome] + [by[(genome, g)] for g in genes])

    print(f"\nWrote {tidy}")
    print(f"Wrote {wide}")


if __name__ == "__main__":
    main()
