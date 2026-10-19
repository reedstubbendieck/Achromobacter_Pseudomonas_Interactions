#!/usr/bin/env bash

set -euo pipefail

PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
RAW="${RNASEQ_RAW:-/mnt/d/Data_Backup/20261005_Pa_RNA-seq}"
WORK="${RNASEQ_WORK:-$RAW/derivedData/rnaseq}"
HERE="$PROJ/scripts"

FFN="$PROJ/references/PAO1/Pseudomonas_aeruginosa_PAO1_107.ffn"
KALLISTO="${KALLISTO:-$HOME/miniforge3/envs/kallisto/bin/kallisto}"

export PROJ RNASEQ_RAW="$RAW" RNASEQ_WORK="$WORK"

mkdir -p "$WORK"/{ref,fastp,kallisto,logs,qc}

echo "=== 1. reference ==="
if [ ! -s "$WORK/ref/PAO1_transcripts.fna" ]; then
  python3 "$HERE/rnaseq/rnaseq_prepare_reference.py" "$FFN" \
    "$WORK/ref/PAO1_transcripts.fna" "$WORK/ref/PAO1_gene_annotation.tsv"
else
  echo "reference already built, skipping"
fi

if [ ! -s "$WORK/ref/PAO1_transcripts.idx" ]; then
  "$KALLISTO" index -i "$WORK/ref/PAO1_transcripts.idx" -k 31 --threads 8 \
    "$WORK/ref/PAO1_transcripts.fna"
else
  echo "index already built, skipping"
fi

echo "=== 2. fastp trimming ==="
bash "$HERE/rnaseq/rnaseq_run_fastp.sh"

echo "=== 3. strandedness check ==="
bash "$HERE/rnaseq/rnaseq_detect_strandedness.sh" M1 400000 \
  | tee "$WORK/qc/strandedness.txt"

echo "=== 4. kallisto quantification ==="
bash "$HERE/rnaseq/rnaseq_run_kallisto.sh"

echo "=== 5. QC summary ==="
python3 "$HERE/rnaseq/rnaseq_qc_summary.py"

echo "=== 6. operon-leader promoters and the Fur box scan ==="
res="$PROJ/derivedData/rnaseq/deseq2_results_coculture_vs_monoculture.tsv"
if [ ! -s "$res" ]; then
  echo "missing $res - knit Achromobacter.Rmd first" >&2
  exit 1
fi

bash "$HERE/rnaseq/rnaseq_fur_box_motifs.sh"

echo "=== done ==="
echo "tables  : $PROJ/derivedData/rnaseq"
echo "scratch : $WORK"
echo "now knit Achromobacter.Rmd"
