#!/usr/bin/env bash
# Run ONCE on the login node (needs internet only if references are not already in $REF_DIR).
# 1. finds (or downloads) human GRCh38 + C. neoformans H99 genome + annotation, 2. faidx, 3. resolves virulence-gene coordinates.
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"
mkdir -p "$REF_DIR" "$LOG_DIR"

# stage_ref <target> <user-supplied-src> <url> <label> <name-regex> <ext-regex>
stage_ref() {
  local target="$1" src="$2" url="$3" label="$4" nre="$5" ere="$6" dl
  if [[ -s "$target" ]]; then log "$label: already present ($target)"; return 0; fi
  if [[ -z "$src" ]]; then
    src="$(find "$REF_DIR" -maxdepth 2 -type f -regextype posix-extended \
           -iregex ".*/[^/]*(${nre})[^/]*\.(${ere})(\.gz)?$" ! -path "$target" 2>/dev/null | head -n1 || true)"
  fi
  if [[ -n "$src" ]]; then
    log "$label: using $src"
    case "$src" in *.gz) gzip -cd "$src" > "$target.tmp" ;; *) cp "$src" "$target.tmp" ;; esac
  else
    dl="$REF_DIR/$(basename "$url")"
    log "$label: downloading $url"
    wget -c --tries=5 --timeout=60 -O "$dl" "$url" \
      || die "$label download failed. Download it manually, put it in $REF_DIR (any name containing '${nre//|/ or }') and re-run."
    gzip -t "$dl" || die "$dl is corrupt; delete it and re-run"
    gzip -cd "$dl" > "$target.tmp"
  fi
  [[ -s "$target.tmp" ]] || die "$label: staged file is empty"
  mv "$target.tmp" "$target"
}

stage_ref "$HUMAN_FA"   "$HUMAN_FA_SRC"    "$HUMAN_URL"      "human genome"        'GRCh38|hg38|CHM13|T2T|human' 'fa|fna|fasta'
stage_ref "$CRYPTO_FA"  "$CRYPTO_FA_SRC"   "$CRYPTO_FA_URL"  "C. neoformans H99"   'H99|CNA3|neoformans|grubii'  'fa|fna|fasta'
stage_ref "$CRYPTO_GFF" "$CRYPTO_GFF_SRC"  "$CRYPTO_GFF_URL" "H99 annotation"      'H99|CNA3|neoformans|grubii'  'gff|gff3'

for f in "$HUMAN_FA" "$CRYPTO_FA"; do
  grep -q '^>' "$f" || die "$f is not a FASTA file"
done
[[ -s "$CRYPTO_GFF" ]] || die "no GFF at $CRYPTO_GFF"

log "indexing FASTA files (samtools faidx)"
sx samtools samtools faidx "$HUMAN_FA"
sx samtools samtools faidx "$CRYPTO_FA"

log "resolving virulence-gene coordinates"
python3 "$HERE/resolve_virulence_genes.py" --gff "$CRYPTO_GFF" --fai "$CRYPTO_FA.fai" \
        --manual "$GENES_MANUAL" --out-bed "$GENES_BED" --out-tsv "$GENES_TSV"

log "references ready:"
ls -lh "$REF_DIR"
echo
echo "Check $GENES_TSV (gene coordinates + product) before continuing."
