#!/usr/bin/env bash
# Processes ONE fastq (one row of outputs/manifest.tsv).  Usage: 03_process_sample.sh <row-number>   (or SLURM array task)
#   raw QC -> NanoFilt -> human depletion (minimap2 vs GRCh38) -> map to H99 -> coverage / identity / variants
#   -> virulence-gene coverage + consensus -> (if enough data) Flye + Medaka
# Sample-level problems (empty file, no reads after filtering ...) are recorded in metrics.tsv and exit 0.
set -euo pipefail
HERE="$(dirname "$(readlink -f "$0")")"
source "$HERE/config.sh"; source "$HERE/lib.sh"

IDX="${1:-${SLURM_ARRAY_TASK_ID:-}}"
[[ -n "$IDX" ]] || die "usage: $0 <manifest-row-number>"
row="$(awk -F'\t' -v n="$IDX" 'NR==n+1' "$MANIFEST")"
[[ -n "$row" ]] || die "manifest row $IDX not found in $MANIFEST"
IFS=$'\t' read -r KEY RUN BARCODE FASTQ SAMPLE_ID BYTES <<<"$row"
THREADS="${SLURM_CPUS_PER_TASK:-$TASK_CPUS}"
ST=$(( THREADS > 4 ? 4 : THREADS ))

S="$SAMPLE_DIR/$KEY"
if [[ -e "$S/.done" && "${FORCE:-0}" != 1 ]]; then log "$KEY already done (FORCE=1 to redo)"; exit 0; fi
rm -rf "$S"; mkdir -p "$S/genes" "$S/assembly"
[[ -s "$HUMAN_MMI" && -s "$CRYPTO_FA.fai" && -s "$GENES_BED" && -s "$GENES_TSV" ]] \
  || die "references not ready (need $HUMAN_MMI, $CRYPTO_FA.fai, $GENES_BED). Run 01_prepare_references.sh and 01b_index_human.sh"
[[ -r "$FASTQ" ]] || die "cannot read $FASTQ"

M="$S/metrics.tsv"; : >"$M"
metric() { printf '%s\t%s\n' "$1" "$2" >>"$M"; }
finish() { metric status "$1"; touch "$S/.done"; log "$KEY finished: $1"; exit 0; }

log "=== $KEY  ($SAMPLE_ID)  $FASTQ ==="
metric sample_key "$KEY"; metric run "$RUN"; metric barcode "$BARCODE"; metric sample_id "$SAMPLE_ID"
metric min_q "$MIN_Q"; metric min_len "$MIN_LEN"; metric min_mapq "$MIN_MAPQ"; metric min_dp "$MIN_DP"

# ---- 1. raw counts ---------------------------------------------------------------------------------
read -r RAW_R RAW_B < <(fq_counts "$FASTQ")
metric raw_reads "$RAW_R"; metric raw_bases "$RAW_B"
(( RAW_R > 0 )) || finish empty_input

# ---- 2. NanoFilt ------------------------------------------------------------------------------------
log "NanoFilt -q $MIN_Q -l $MIN_LEN"
zcat -f -- "$FASTQ" | sx nanofilt NanoFilt -q "$MIN_Q" -l "$MIN_LEN" | gzip -c >"$S/filt.fastq.gz"
read -r FILT_R FILT_B < <(fq_counts "$S/filt.fastq.gz")
metric filtered_reads "$FILT_R"; metric filtered_bases "$FILT_B"
(( FILT_R > 0 )) || finish no_reads_after_filter

# ---- 3. human depletion -----------------------------------------------------------------------------
log "human depletion (minimap2 map-ont vs GRCh38; keeping reads with NO human alignment)"
sx minimap2 minimap2 -ax map-ont -t "$THREADS" --secondary=no "$HUMAN_MMI" "$S/filt.fastq.gz" \
  | sx samtools samtools view -u -f 4 - \
  | sx samtools samtools fastq -0 "$S/nohuman.fastq" - >/dev/null
if [[ -s "$S/nohuman.fastq" ]]; then gzip -f "$S/nohuman.fastq"; else rm -f "$S/nohuman.fastq"; : | gzip -c >"$S/nohuman.fastq.gz"; fi
read -r NH_R NH_B < <(fq_counts "$S/nohuman.fastq.gz")
metric nohuman_reads "$NH_R"; metric nohuman_bases "$NH_B"
metric human_reads_removed "$(( FILT_R - NH_R ))"
metric human_bases_removed "$(( FILT_B - NH_B ))"
rm -f "$S/filt.fastq.gz"
(( NH_R > 0 )) || finish no_nonhuman_reads

# ---- 4. map non-human reads to C. neoformans H99 ------------------------------------------------------------
BAM="$S/crypto.sorted.bam"
log "mapping non-human reads to H99"
sx minimap2 minimap2 -ax map-ont -t "$THREADS" --secondary=no "$CRYPTO_FA" "$S/nohuman.fastq.gz" \
  | sx samtools samtools sort -@ "$ST" -T "$S/sorttmp" -o "$BAM" -
sx samtools samtools index "$BAM"
sx samtools samtools stats -F 2304 "$BAM" >"$S/crypto.stats.txt"
metric crypto_mapped_reads "$(sx samtools samtools view -c -F 3844 "$BAM")"
metric crypto_mapped_reads_mapq "$(sx samtools samtools view -c -F 3844 -q "$MIN_MAPQ" "$BAM")"

# ---- 5. coverage (whole genome + virulence genes) -------------------------------------------------------
log "mosdepth"
sx mosdepth mosdepth -t "$ST" -Q "$MIN_MAPQ" -F 3844 --by "$GENES_BED" --thresholds 1,3,5,10 "$S/cov" "$BAM"
zcat "$S/cov.per-base.bed.gz" | awk -v OFS='\t' 'BEGIN{split("1 3 5 10 20",T," ")}
  {L=$3-$2; tot+=L; for(i=1;i<=5;i++) if($4>=T[i]) c[i]+=L}
  END{for(i=1;i<=5;i++) printf "%s\t%.6f\n", T[i], (tot?c[i]/tot:0)}' >"$S/breadth.tsv"
zcat "$S/cov.per-base.bed.gz" | awk -v m="$MIN_DP" -v OFS='\t' '$4<m{print $1,$2,$3}' >"$S/lowcov.bed"
COVERED="$(zcat "$S/cov.per-base.bed.gz" | awk -v m="$MIN_DP" '$4>=m{c+=$3-$2} END{print c+0}')"

# ---- 6. per-read identity to H99 (1 - minimap2 'de' divergence) ---------------------------------------
sx samtools samtools view -F 3844 -q "$MIN_MAPQ" "$BAM" \
  | awk '{for(i=12;i<=NF;i++) if($i ~ /^de:f:/){split($i,a,":"); print 1-a[3]; break}}' | sort -g >"$S/id.tmp"
awk -v OFS='\t' '{v[NR]=$1; s+=$1} END{ if(NR==0){print "NA","NA",0} else {m=(NR%2)?v[(NR+1)/2]:(v[NR/2]+v[NR/2+1])/2; printf "%.5f\t%.5f\t%d\n", s/NR, m, NR} }' \
  "$S/id.tmp" >"$S/read_identity.tsv"
rm -f "$S/id.tmp"

# ---- 7. variants vs H99 (haploid) ---------------------------------------------------------------------------
log "bcftools variant calling (haploid, depth >= $MIN_DP)"
sx bcftools bcftools mpileup -Ou -f "$CRYPTO_FA" -q "$MIN_MAPQ" -d 500 -a FORMAT/AD,FORMAT/DP --threads "$ST" "$BAM" \
  | sx bcftools bcftools call -mv --ploidy 1 -Ou --threads "$ST" \
  | sx bcftools bcftools filter -i "QUAL>=20 && INFO/DP>=${MIN_DP}" -Oz -o "$S/variants.vcf.gz"
sx bcftools bcftools index -f "$S/variants.vcf.gz"
NSNP="$(sx bcftools bcftools view -H -v snps "$S/variants.vcf.gz" | wc -l)"
NIND="$(sx bcftools bcftools view -H -v indels "$S/variants.vcf.gz" | wc -l)"
printf '%s\t%s\t%s\n' "$NSNP" "$NIND" "$COVERED" >"$S/variants_summary.tsv"

# ---- 8. virulence genes: reads, variants, consensus ------------------------------------------------------------
log "virulence genes"
GO="$S/genes"; : >"$GO/gene_reads.tsv"; : >"$GO/gene_variants.tsv"
while IFS=$'\t' read -r sym chrom gs ge strand _rest <&3; do
  [[ "$sym" == "symbol" || -z "$sym" ]] && continue
  reg="${chrom}:${gs}-${ge}"
  nr="$(sx samtools samtools view -c -F 3844 -q "$MIN_MAPQ" "$BAM" "$reg" </dev/null)"
  nv="$(sx bcftools bcftools view -H -r "$reg" "$S/variants.vcf.gz" </dev/null | wc -l)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sym" "$chrom" "$gs" "$ge" "$strand" "$nr" "$nv" >>"$GO/gene_reads.tsv"
  sx bcftools bcftools query -r "$reg" -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%TYPE\n' "$S/variants.vcf.gz" </dev/null \
    | awk -v s="$sym" -v OFS='\t' '{print s,$0}' >>"$GO/gene_variants.tsv"
  hdr=">${KEY}|${sym}|${chrom}:${gs}-${ge}(${strand})"
  sx samtools samtools faidx "$CRYPTO_FA" "$reg" </dev/null \
    | sx bcftools bcftools consensus -m "$S/lowcov.bed" "$S/variants.vcf.gz" >"$GO/$sym.fwd.fa"
  if [[ "$strand" == "-" ]]; then
    sx seqkit seqkit seq -t DNA -r -p "$GO/$sym.fwd.fa" </dev/null | awk -v h="$hdr" 'NR==1{print h; next}{print}' >"$GO/$sym.consensus.fa"
  else
    awk -v h="$hdr" 'NR==1{print h; next}{print}' "$GO/$sym.fwd.fa" >"$GO/$sym.consensus.fa"
  fi
  rm -f "$GO/$sym.fwd.fa"
done 3<"$GENES_TSV"

# ---- 9. assembly of the C. neoformans reads (only if there is enough data; failures are recorded, not fatal) -----------
MAPPED_BASES="$(awk -F'\t' '$1=="SN" && $2=="bases mapped (cigar):"{print $3}' "$S/crypto.stats.txt")"
MAPPED_BASES="${MAPPED_BASES:-0}"
ASM="$S/assembly"; ASM_STATUS="skipped_low_data"
if (( MAPPED_BASES >= ASM_MIN_MAPPED_BASES )); then
  log "assembly: extracting H99-mapped reads ($MAPPED_BASES bp)"
  sx samtools samtools fastq -F 3844 -0 "$S/crypto_reads.fastq" "$BAM" >/dev/null
  gzip -f "$S/crypto_reads.fastq"
  FINAL=""
  if sx flye flye --nano-hq "$S/crypto_reads.fastq.gz" --out-dir "$ASM/flye" --threads "$THREADS" >"$ASM/flye.log" 2>&1 \
     && [[ -s "$ASM/flye/assembly.fasta" ]]; then
    FINAL="$ASM/flye/assembly.fasta"; ASM_STATUS="flye_ok"
    if sx medaka medaka_consensus -i "$S/crypto_reads.fastq.gz" -d "$ASM/flye/assembly.fasta" -o "$ASM/medaka" \
         -t "$ST" -m "$MEDAKA_MODEL" >"$ASM/medaka.log" 2>&1 && [[ -s "$ASM/medaka/consensus.fasta" ]]; then
      FINAL="$ASM/medaka/consensus.fasta"; ASM_STATUS="flye_medaka_ok"
    else
      warn "$KEY: Medaka failed (see $ASM/medaka.log); using unpolished Flye assembly"; ASM_STATUS="flye_ok_medaka_failed"
    fi
  else
    warn "$KEY: Flye did not produce an assembly (expected at low coverage; see $ASM/flye.log)"; ASM_STATUS="flye_failed"
  fi
  if [[ -n "$FINAL" ]]; then
    sx seqkit seqkit stats -a -T "$FINAL" >"$ASM/assembly_stats.tsv" || warn "$KEY: seqkit stats failed"
    if sx minimap2 minimap2 -x asm20 -c --secondary=no -t "$ST" "$CRYPTO_FA" "$FINAL" >"$ASM/asm_vs_ref.paf" 2>"$ASM/asm_vs_ref.log"; then
      awk -F'\t' -v OFS='\t' '{m+=$10; b+=$11; q+=$4-$3; n++} END{ if(b>0) printf "%.5f\t%d\t%d\n", m/b, q, n; else print "NA",0,0 }' \
        "$ASM/asm_vs_ref.paf" >"$ASM/asm_vs_ref.tsv"
    fi
  fi
fi
metric assembly_status "$ASM_STATUS"
finish ok
