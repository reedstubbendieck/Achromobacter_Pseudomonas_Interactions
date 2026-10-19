#!/usr/bin/env bash

set -euo pipefail

PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
RAW="${RNASEQ_RAW:-/mnt/d/Data_Backup/20261005_Pa_RNA-seq}"
WORK="${RNASEQ_WORK:-$RAW/derivedData/rnaseq}"
FASTP="$HOME/miniconda3/envs/fastp/bin/fastp"
THREADS=16

mkdir -p "$WORK/fastp" "$WORK/logs"

tail -n +2 "$PROJ/rawData/rnaseq/samples.tsv" | while IFS=$'\t' read -r sample cond rep lib pfx; do
  [ -n "${sample:-}" ] || continue
  out1="$WORK/fastp/${sample}_R1.trim.fastq.gz"
  out2="$WORK/fastp/${sample}_R2.trim.fastq.gz"

  if [ -s "$out1" ] && [ -s "$out2" ] && [ -s "$WORK/fastp/${sample}.fastp.json" ]; then
    echo "[$(date +%H:%M:%S)] $sample already trimmed, skipping"
    continue
  fi

  in1="$RAW/${sample}_R1.fastq.gz"
  in2="$RAW/${sample}_R2.fastq.gz"
  for f in "$in1" "$in2"; do
    [ -s "$f" ] || { echo "missing raw read: $f" >&2; exit 1; }
  done

  echo "[$(date +%H:%M:%S)] fastp $sample ($cond rep$rep, library $lib)"
  "$FASTP" \
    --in1 "$in1" \
    --in2 "$in2" \
    --out1 "$out1" \
    --out2 "$out2" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --cut_front --cut_tail --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 30 \
    --n_base_limit 5 \
    --length_required 31 \
    --thread "$THREADS" \
    --json "$WORK/fastp/${sample}.fastp.json" \
    --html "$WORK/fastp/${sample}.fastp.html" \
    --report_title "fastp: $sample ($cond rep$rep)" \
    2> "$WORK/logs/${sample}.fastp.log"
done

echo "[$(date +%H:%M:%S)] fastp done"
