#!/usr/bin/env bash
# Scans EVERY $BASE/demultiplexed_output* folder and writes outputs/manifest.tsv (one row per barcode fastq).
# Optional $BASE/sample_map.tsv (tab-separated: run, barcode, sample_id) renames barcodes to sample IDs (e.g. MC81).
# Re-run it any time after editing sample_map.tsv; results are keyed by run+barcode so nothing needs recomputing.
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"
mkdir -p "$OUT_DIR"
tmp="$(mktemp)"; mapf="$(mktemp)"; trap 'rm -f "$tmp" "$mapf"' EXIT

while IFS= read -r d; do
  dname="$(basename "$d")"
  lbl="$(printf '%s' "$dname" | sed -E 's/^[Dd][Ee][Mm][Uu][Ll][Tt][Ii][Pp][Ll][Ee][Xx][Ee][Dd]_[Oo][Uu][Tt][Pp][Uu][Tt]_?//' \
         | tr 'A-Z' 'a-z' | sed -E 's/[^a-z0-9]+/_/g; s/^_+|_+$//g')"
  [[ -n "$lbl" ]] || lbl="$(printf '%s' "$dname" | tr 'A-Z' 'a-z' | sed -E 's/[^a-z0-9]+/_/g')"
  run="run_${lbl}"
  while IFS= read -r f; do
    bc="$(basename "$f" | sed -E 's/\.(fastq|fq)(\.gz)?$//; s/^SQK-[A-Za-z0-9-]+_//; s/[^A-Za-z0-9]+/_/g')"
    printf '%s\t%s\t%s\t%s\t%s\n' "$run" "$bc" "$f" "$(stat -c %s "$f")" "$dname" >>"$tmp"
  done < <(find "$d" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | sort)
done < <(find "$BASE" -maxdepth 1 -mindepth 1 -type d -name "$RUN_GLOB" | sort)

[[ -s "$tmp" ]] || die "no fastq files found under $BASE/$RUN_GLOB"
[[ -s "$SAMPLE_MAP" ]] && tr -d '\r' <"$SAMPLE_MAP" >"$mapf" || : >"$mapf"

{
  printf 'sample_key\trun\tbarcode\tfastq\tsample_id\tbytes\n'
  awk -F'\t' -v OFS='\t' '
    function nb(x){ x=tolower(x); gsub(/[^a-z0-9]+/,"_",x); if (x ~ /^[0-9]+$/) x=sprintf("barcode%02d", x+0); return x }
    function nr(x){ x=tolower(x); sub(/^demultiplexed_output_?/,"",x); gsub(/[^a-z0-9]+/,"_",x); sub(/^_+/,"",x); sub(/_+$/,"",x); return "run_" x }
    FILENAME==ARGV[1] { if (tolower($1)!="run" && $3!="") { m[nr($1) SUBSEP nb($2)]=$3 } next }
    {
      key=$1 "__" $2; sid=""
      if ((nr($1) SUBSEP nb($2)) in m) sid=m[nr($1) SUBSEP nb($2)]
      else if ((nr($5) SUBSEP nb($2)) in m) sid=m[nr($5) SUBSEP nb($2)]
      if (sid=="") sid=key
      gsub(/[ \t]+/,"_",sid)
      print key, $1, $2, $3, sid, $4
    }' "$mapf" "$tmp"
} >"$MANIFEST"

n=$(( $(wc -l <"$MANIFEST") - 1 ))
log "manifest: $n fastq files from $(cut -f2 "$MANIFEST" | tail -n +2 | sort -u | wc -l) run folders -> $MANIFEST"
if [[ ! -s "$SAMPLE_MAP" ]]; then
  tpl="$BASE/sample_map_TEMPLATE.tsv"
  { printf 'run\tbarcode\tsample_id\n'; tail -n +2 "$MANIFEST" | awk -F'\t' -v OFS='\t' '{print $2,$3,""}'; } >"$tpl"
  log "no sample_map.tsv yet. Fill the 3rd column of $tpl (e.g. MC81), save it as $SAMPLE_MAP and re-run this script."
fi
cut -f2 "$MANIFEST" | tail -n +2 | sort | uniq -c | sed 's/^/   files per run: /' >&2
