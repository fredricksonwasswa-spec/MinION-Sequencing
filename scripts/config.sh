#!/usr/bin/env bash
# =============================================================================
# config.sh  -  the ONLY file you normally need to edit. Every script sources it.
# =============================================================================

# ---- layout ------------------------------------------------------------------
export BASE="${BASE:-$HOME/Dr_Fred_Data}"
export SCRIPTS="${SCRIPTS:-$BASE/scripts}"
export CONT_DIR="$BASE/containers"
export REF_DIR="$BASE/references"
export OUT_DIR="$BASE/outputs"
export LOG_DIR="$BASE/logs"
export SAMPLE_DIR="$OUT_DIR/samples"
export SUMMARY_DIR="$OUT_DIR/summary"
export MANIFEST="$OUT_DIR/manifest.tsv"
export SAMPLE_MAP="${SAMPLE_MAP:-$BASE/sample_map.tsv}"   # optional: run<TAB>barcode<TAB>sample_id
export RUN_GLOB="demultiplexed_output*"                    # every matching folder in $BASE is analysed

# ---- references (auto-found in $REF_DIR, or downloaded, by 01_prepare_references.sh) ----
export HUMAN_FA="$REF_DIR/human_GRCh38.fa"
export HUMAN_MMI="$REF_DIR/human_GRCh38.map-ont.mmi"
export CRYPTO_FA="$REF_DIR/Cneoformans_H99.fa"
export CRYPTO_GFF="$REF_DIR/Cneoformans_H99.gff"
export GENES_TSV="$REF_DIR/virulence_genes.tsv"            # written by 01_prepare_references.sh
export GENES_BED="$REF_DIR/virulence_genes.bed"
export GENES_MANUAL="$REF_DIR/virulence_genes_manual.tsv"  # optional override: symbol chrom start end strand
# Optional: point at files you already have (they are COPIED into $REF_DIR)
export HUMAN_FA_SRC="${HUMAN_FA_SRC:-}"
export CRYPTO_FA_SRC="${CRYPTO_FA_SRC:-}"
export CRYPTO_GFF_SRC="${CRYPTO_GFF_SRC:-}"
# Download locations used only if nothing suitable is already in $REF_DIR
export HUMAN_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
export CRYPTO_FA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz"
export CRYPTO_GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/245/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.gff.gz"

# ---- analysis parameters ---------------------------------------------------------
export MIN_Q=10                 # NanoFilt mean read quality
export MIN_LEN=500              # NanoFilt minimum length (matches the 500 bp minimum in the original Table 2)
export MIN_MAPQ=10              # MAPQ used for coverage, read counts, variant calling
export MIN_DP=3                 # minimum depth to call a base / report a consensus base (lower = N)
export ASM_MIN_MAPPED_BASES=1000000   # attempt Flye only if >= this many bp of reads map to H99
export MEDAKA_MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_sup_v5.0.0}"   # must match the Dorado model used to basecall
export GENE_MIN_CALLED=0.25     # a sample's gene consensus enters the alignment if >= this fraction of bases are called

# ---- tool versions (quay.io/biocontainers; build tag is resolved automatically) --------
export V_NANOFILT=2.8.0
export V_SEQKIT=2.8.2
export V_MINIMAP2=2.28
export V_SAMTOOLS=1.20
export V_BCFTOOLS=1.20
export V_MOSDEPTH=0.3.8
export V_FLYE=2.9.5
export V_MEDAKA=2.0.1
export V_MAFFT=7.525
# To pin an exact image instead, e.g.:  export URI_FLYE="docker://quay.io/biocontainers/flye:2.9.5--py310h..."

# ---- SLURM (partition names from your cluster: shared / exclusive / gpu) ---------------
export PARTITION="${PARTITION:-shared}"
export ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-6}"
export TASK_CPUS="${TASK_CPUS:-8}"
export TASK_MEM="${TASK_MEM:-32G}"
export TASK_TIME="${TASK_TIME:-24:00:00}"

# ---- container runtime -------------------------------------------------------------
# If your cluster needs it, uncomment:  module load singularity
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$HOME/singularity_cache}"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
mkdir -p "$SINGULARITY_CACHEDIR" 2>/dev/null || true
