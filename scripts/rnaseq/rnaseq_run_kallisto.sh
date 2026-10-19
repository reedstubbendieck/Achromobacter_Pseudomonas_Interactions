#!/usr/bin/env bash

set -euo pipefail

PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
RAW="${RNASEQ_RAW:-/mnt/d/Data_Backup/20261005_Pa_RNA-seq}"
WORK="${RNASEQ_WORK:-$RAW/derivedData/rnaseq}"
KALLISTO="$HOME/miniforge3/envs/kallisto/bin/kallisto"
IDX="$WORK/ref/PAO1_transcripts.idx"
THREADS=16
BOOTSTRAPS=100
STRAND="--rf-stranded"

mkdir -p "$WORK/kallisto" "$WORK/logs"

tail -n +2 "$PROJ/rawData/rnaseq/samples.tsv" | while IFS=$'\t' read -r sample cond rep lib pfx; do
  [ -n "${sample:-}" ] || continue
  out="$WORK/kallisto/$sample"

  if [ -s "$out/abundance.tsv" ] && [ -s "$out/run_info.json" ]; then
    echo "[$(date +%H:%M:%S)] $sample already quantified, skipping"
    continue
  fi

  echo "[$(date +%H:%M:%S)] kallisto quant $sample ($cond rep$rep)"
  mkdir -p "$out"
  "$KALLISTO" quant \
    -i "$IDX" \
    -o "$out" \
    --threads "$THREADS" \
    --bootstrap-samples "$BOOTSTRAPS" \
    $STRAND \
    "$WORK/fastp/${sample}_R1.trim.fastq.gz" \
    "$WORK/fastp/${sample}_R2.trim.fastq.gz" \
    > "$WORK/logs/${sample}.kallisto.log" 2>&1
done

echo "[$(date +%H:%M:%S)] kallisto done"
