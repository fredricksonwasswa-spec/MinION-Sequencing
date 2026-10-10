#!/usr/bin/env bash
# Shared helpers (sourced after config.sh)

log()  { echo "[$(date '+%F %T')] $*" >&2; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR: $*"; exit 1; }

CONTAINER_CMD="${CONTAINER_CMD:-$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null || true)}"
export CONTAINER_CMD

need_runtime() { [[ -n "$CONTAINER_CMD" ]] || die "neither apptainer nor singularity found in PATH (try: module load singularity)"; }

# sif <tool>  -> path of the image for that tool
sif() {
  local v="V_${1^^}"
  [[ -n "${!v:-}" ]] || die "no version set for tool '$1' (V_${1^^} in config.sh)"
  echo "$CONT_DIR/${1}_${!v}.sif"
}

# sx <tool> <command...>  -> run a command inside that tool's container
sx() {
  local tool="$1"; shift
  need_runtime
  "$CONTAINER_CMD" exec -B "$BASE" "$(sif "$tool")" "$@"
}

# pull_tool <tool>  -> download the container once; retries; build tag auto-resolved from quay.io
pull_tool() {
  need_runtime
  local tool="$1" v="V_${1^^}" u="URI_${1^^}" out uri tmp try
  out="$(sif "$tool")"
  if [[ -s "$out" ]]; then log "have $out"; return 0; fi
  uri="${!u:-}"
  if [[ -z "$uri" ]]; then
    uri="$(python3 "$SCRIPTS/resolve_biocontainer.py" "$tool" "${!v}")" \
      || die "cannot resolve a quay.io tag for $tool ${!v}. Set URI_${tool^^} in config.sh, or copy a .sif to $out"
  fi
  tmp="$CONT_DIR/.${tool}.partial.sif"
  for try in 1 2 3; do
    log "pulling $tool  <-  $uri  (attempt $try/3)"
    rm -f "$tmp"
    if "$CONTAINER_CMD" pull "$tmp" "$uri"; then mv "$tmp" "$out"; log "saved $out"; return 0; fi
    sleep 20
  done
  die "failed to pull $tool. Download it on a node with internet and place it at $out"
}

# fq_counts <fastq[.gz]>  -> "reads<TAB>bases"   (robust to empty files; no container needed)
fq_counts() {
  if [[ ! -s "$1" ]]; then printf '0\t0\n'; return 0; fi
  zcat -f -- "$1" | awk 'NR%4==2{n++; b+=length($0)} END{printf "%d\t%d\n", n+0, b+0}'
}
