#!/usr/bin/env bash
# Run ONCE on a node with internet (login node). Pulls every tool as a Singularity/Apptainer image.
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"
mkdir -p "$CONT_DIR" "$LOG_DIR"
for t in nanofilt seqkit minimap2 samtools bcftools mosdepth flye medaka mafft; do
  pull_tool "$t"
done
log "all containers ready in $CONT_DIR"
ls -lh "$CONT_DIR"/*.sif
