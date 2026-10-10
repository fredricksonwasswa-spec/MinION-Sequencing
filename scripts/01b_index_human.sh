#!/usr/bin/env bash
# Builds the minimap2 index of the human genome (needs ~12 GB RAM). Submitted by run_all.sh.
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"
[[ -s "$HUMAN_FA" ]] || die "missing $HUMAN_FA - run 01_prepare_references.sh first"
if [[ -s "$HUMAN_MMI" ]]; then log "index exists: $HUMAN_MMI"; exit 0; fi
sx minimap2 minimap2 -x map-ont -t "${SLURM_CPUS_PER_TASK:-8}" -I 12G -d "$HUMAN_MMI.tmp" "$HUMAN_FA"
mv "$HUMAN_MMI.tmp" "$HUMAN_MMI"
log "wrote $HUMAN_MMI"
