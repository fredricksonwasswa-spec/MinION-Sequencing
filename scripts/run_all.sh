#!/usr/bin/env bash
# Usage:  bash run_all.sh            submit everything to SLURM (index -> per-sample array -> aggregate)
#         bash run_all.sh --local    run sequentially in this shell (small tests / interactive node)
#         bash run_all.sh --aggregate-only
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"
MODE="${1:-slurm}"
mkdir -p "$LOG_DIR" "$OUT_DIR"

for f in "$CRYPTO_FA.fai" "$HUMAN_FA.fai" "$GENES_BED" "$GENES_TSV"; do
  [[ -s "$f" ]] || die "missing $f - run 01_prepare_references.sh first"
done
for t in nanofilt seqkit minimap2 samtools bcftools mosdepth flye medaka mafft; do
  [[ -s "$(sif "$t")" ]] || die "missing container $(sif "$t") - run 00_pull_containers.sh first"
done

bash "$HERE/02_make_manifest.sh"
N=$(( $(wc -l <"$MANIFEST") - 1 ))

case "$MODE" in
  --aggregate-only) bash "$HERE/05_aggregate.sh"; exit 0 ;;
  --local)
    bash "$HERE/01b_index_human.sh"
    for i in $(seq 1 "$N"); do bash "$HERE/03_process_sample.sh" "$i"; done
    bash "$HERE/05_aggregate.sh"
    ;;
  slurm|--slurm)
    command -v sbatch >/dev/null || die "sbatch not found; use --local"
    dep=""
    if [[ ! -s "$HUMAN_MMI" ]]; then
      jidx=$(sbatch --parsable -p "$PARTITION" -c 8 --mem=32G -t 04:00:00 -J hg_index -o "$LOG_DIR/hg_index_%j.out" \
             --wrap "bash $HERE/01b_index_human.sh")
      dep="--dependency=afterok:$jidx"; log "submitted human index job $jidx"
    fi
    jarr=$(sbatch --parsable -p "$PARTITION" -c "$TASK_CPUS" --mem="$TASK_MEM" -t "$TASK_TIME" $dep \
           -J fred_sample --array="1-${N}%${ARRAY_CONCURRENCY}" -o "$LOG_DIR/fred_sample_%A_%a.out" \
           --wrap "bash $HERE/03_process_sample.sh")
    log "submitted per-sample array $jarr ($N tasks, max ${ARRAY_CONCURRENCY} at once)"
    jagg=$(sbatch --parsable -p "$PARTITION" -c 4 --mem=16G -t 04:00:00 --dependency=afterany:"$jarr" \
           -J fred_aggregate -o "$LOG_DIR/fred_aggregate_%j.out" --wrap "bash $HERE/05_aggregate.sh")
    log "submitted aggregation job $jagg (runs when the array finishes)"
    echo "Monitor:  squeue -u \$USER     Logs: $LOG_DIR/"
    echo "Re-run failed samples:  FORCE=1 sbatch ... or  bash $HERE/03_process_sample.sh <row>"
    ;;
  *) die "unknown option $MODE" ;;
esac
