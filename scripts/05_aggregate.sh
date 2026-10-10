#!/usr/bin/env bash
# Builds the summary tables and the MAFFT alignments (reference gene + every sample consensus).
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"
mkdir -p "$SUMMARY_DIR/alignments"
python3 "$HERE/aggregate.py" --stage tables --manifest "$MANIFEST" --samples "$SAMPLE_DIR" --summary "$SUMMARY_DIR" \
        --crypto-fa "$CRYPTO_FA" --genes-tsv "$GENES_TSV" --min-called "$GENE_MIN_CALLED"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
for f in "$SUMMARY_DIR"/alignments/*.unaligned.fa; do
  [[ -e "$f" ]] || continue
  n="$(grep -c '^>' "$f" || true)"
  out="${f%.unaligned.fa}.aln.fa"
  if (( n >= 2 )); then
    log "MAFFT $(basename "$f") ($n sequences)"
    sx mafft mafft --auto --quiet --thread "$THREADS" "$f" >"$out"
  else
    log "skip $(basename "$f"): no sample has a usable consensus for this gene"
    rm -f "$out"
  fi
done
python3 "$HERE/aggregate.py" --stage align --summary "$SUMMARY_DIR"
log "done. Key files: $SUMMARY_DIR/{sample_summary,gene_summary,gene_variants,gene_alignment_identity}.tsv and alignments/*.aln.fa"
