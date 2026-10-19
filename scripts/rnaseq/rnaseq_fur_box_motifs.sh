#!/usr/bin/env bash

set -euo pipefail

PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
WORK="${RNASEQ_WORK:-/mnt/d/Data_Backup/20261005_Pa_RNA-seq/derivedData/rnaseq}"
M="$WORK/motifs"
OUT="$PROJ/derivedData/rnaseq"
export PATH="$HOME/miniforge3/envs/meme/bin:$PATH"
RSCRIPT="$HOME/miniforge3/envs/rnaseq-de/bin/Rscript"
THRESH=1e-4

python3 "$PROJ/scripts/rnaseq/rnaseq_promoter_motifs.py"
cd "$M"

echo "=== STREME: de novo, induced vs background ==="
streme --p promoters_induced.fa --n promoters_background.fa --oc streme_induced \
       --dna --nmotifs 5 --minw 8 --maxw 22 > streme.log 2>&1
grep -E "^MOTIF" streme_induced/streme.txt

echo
echo "=== FIMO: published Fur box consensus ==="
iupac2meme -dna GATAATGATAATCATTATC > furbox.meme
declare -A HIT TOT
for s in induced repressed background; do
  fimo --verbosity 1 --thresh "$THRESH" --text furbox.meme "promoters_$s.fa" \
       > "fimo_$s.txt" 2>/dev/null
  fimo --verbosity 1 --thresh 0.5 --text furbox.meme "promoters_$s.fa" \
       > "fimo_scores_$s.txt" 2>/dev/null
  HIT[$s]=$(awk 'NR>1 && NF>3{print $3}' "fimo_$s.txt" | sort -u | wc -l)
  TOT[$s]=$(grep -c "^>" "promoters_$s.fa")
  printf "%-11s %4d / %4d promoters with a Fur box\n" "$s" "${HIT[$s]}" "${TOT[$s]}"
done

echo
echo "=== FIMO: P. aeruginosa Fur box PWM built from CollecTF sites ==="
awk -F, 'NR > 1 && length($5) == 19 {printf ">%s_%s\n%s\n", $1, NR - 1, $5}' \
  "$PROJ/references/PAO1/fur_boxes_collectf.csv" > collectf_sites.fa
echo "CollecTF sites: $(grep -c '^>' collectf_sites.fa)"
meme collectf_sites.fa -dna -mod oops -nmotifs 1 -w 19 -revcomp \
     -oc meme_collectf > meme_collectf.log 2>&1
grep -E "^MOTIF" meme_collectf/meme.txt | sed 's/^/  /'

declare -A PAHIT
for s in induced repressed background; do
  fimo --verbosity 1 --thresh "$THRESH" --text meme_collectf/meme.txt \
       "promoters_$s.fa" > "fimo_pa_$s.txt" 2>/dev/null
  PAHIT[$s]=$(awk 'NR>1 && NF>3{print $3}' "fimo_pa_$s.txt" | sort -u | wc -l)
  printf "%-11s %4d / %4d promoters with a Fur box (PA PWM)\n" \
         "$s" "${PAHIT[$s]}" "${TOT[$s]}"
done

cat > fur_fisher.R <<RRR
one <- function(h, t, hb, tb, alt, tag) {
  x <- c(h, t - h); y <- c(hb, tb - hb)
  f <- fisher.test(matrix(c(x, y), 2), alternative = alt)
  cat(sprintf("%-28s %5.1f%% vs %4.1f%%  OR = %6.2f  p = %.3g\n",
              tag, 100*x[1]/sum(x), 100*y[1]/sum(y), f\$estimate, f\$p.value))
}
cat("published E. coli consensus\n")
one(${HIT[induced]}, ${TOT[induced]}, ${HIT[background]}, ${TOT[background]},
    "greater", "  coculture vs unchanged")
one(${HIT[repressed]}, ${TOT[repressed]}, ${HIT[background]}, ${TOT[background]},
    "two.sided", "  monoculture vs unchanged")
cat("\nCollecTF P. aeruginosa PWM\n")
one(${PAHIT[induced]}, ${TOT[induced]}, ${PAHIT[background]}, ${TOT[background]},
    "greater", "  coculture vs unchanged")
one(${PAHIT[repressed]}, ${TOT[repressed]}, ${PAHIT[background]}, ${TOT[background]},
    "two.sided", "  monoculture vs unchanged")
RRR

{
  echo "Fur box in coculture-induced promoters - PAO1 RNA-seq"
  echo "===================================================="
  echo "promoters: 300 bp upstream to 50 bp into the gene, operon leaders only"
  echo "motif: published Fur box consensus GATAATGATAATCATTATC (iupac2meme PWM)"
  echo "FIMO p-value threshold: $THRESH"
  echo
  for s in induced repressed background; do
    printf "%-11s %4d / %4d promoters with a Fur box\n" "$s" "${HIT[$s]}" "${TOT[$s]}"
  done
  echo
  echo "rescanned with a PWM built by meme from the $(grep -c '^>' collectf_sites.fa) CollecTF-curated"
  echo "P. aeruginosa Fur boxes (references/PAO1/fur_boxes_collectf.csv):"
  for s in induced repressed background; do
    printf "%-11s %4d / %4d promoters with a Fur box\n" "$s" "${PAHIT[$s]}" "${TOT[$s]}"
  done
  echo
  "$RSCRIPT" fur_fisher.R
  echo
  echo "STREME de novo motifs (induced vs background):"
  grep -E "^MOTIF" streme_induced/streme.txt | sed 's/^/  /'
} > "$OUT/fur_box_summary.txt"


echo
cat "$OUT/fur_box_summary.txt"
echo
echo "written -> $OUT/fur_box_summary.txt"
