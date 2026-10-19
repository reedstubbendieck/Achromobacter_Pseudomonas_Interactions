#!/usr/bin/env bash
set -euo pipefail

RAW="${RNASEQ_RAW:-/mnt/d/Data_Backup/20261005_Pa_RNA-seq}"
WORK="${RNASEQ_WORK:-$RAW/derivedData/rnaseq}"
KALLISTO="$HOME/miniforge3/envs/kallisto/bin/kallisto"
IDX="$WORK/ref/PAO1_transcripts.idx"
SAMPLE="${1:-M1}"
NREADS="${2:-400000}"

TMP="$WORK/qc/strand_${SAMPLE}"
mkdir -p "$TMP"

set +o pipefail
zcat "$WORK/fastp/${SAMPLE}_R1.trim.fastq.gz" | head -n $((NREADS * 4)) > "$TMP/sub_R1.fastq"
zcat "$WORK/fastp/${SAMPLE}_R2.trim.fastq.gz" | head -n $((NREADS * 4)) > "$TMP/sub_R2.fastq"
set -o pipefail

for mode in unstranded fr-stranded rf-stranded; do
  flag=""
  [ "$mode" = "fr-stranded" ] && flag="--fr-stranded"
  [ "$mode" = "rf-stranded" ] && flag="--rf-stranded"
  "$KALLISTO" quant -i "$IDX" -o "$TMP/$mode" --threads 8 $flag \
    "$TMP/sub_R1.fastq" "$TMP/sub_R2.fastq" > "$TMP/$mode.log" 2>&1
done

echo "sample: $SAMPLE   subsample: $NREADS pairs"
printf "%-14s %12s %14s %8s\n" mode processed pseudoaligned pct
for mode in unstranded fr-stranded rf-stranded; do
  python3 -c '
import json, sys
d = json.load(open(sys.argv[1]))
print("%-14s %12d %14d %7.2f%%" % (
    sys.argv[2], d["n_processed"], d["n_pseudoaligned"],
    100.0 * d["n_pseudoaligned"] / d["n_processed"]))' "$TMP/$mode/run_info.json" "$mode"
done
rm -f "$TMP/sub_R1.fastq" "$TMP/sub_R2.fastq"
