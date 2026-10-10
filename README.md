# Dr_Fred_Data re-analysis: direct MinION sequencing of C. neoformans from CSF (containerised)

Copy `scripts/` to `~/Dr_Fred_Data/scripts/`.

## Run order
    cd ~/Dr_Fred_Data/scripts
    bash 00_pull_containers.sh        # once, login node (internet): pulls 9 Singularity images
    bash 01_prepare_references.sh     # once, login node: finds/downloads GRCh38 + H99, resolves the 6 virulence genes
    cat ../references/virulence_genes.tsv     # CHECK coordinates / product
    bash 02_make_manifest.sh          # lists every barcode in every demultiplexed_output* folder
    #  fill ~/Dr_Fred_Data/sample_map_TEMPLATE.tsv (barcode -> MC81 ...) and save as sample_map.tsv
    bash run_all.sh                   # SLURM: human index -> per-sample array -> aggregation
    #  bash run_all.sh --local            sequential run
    #  bash run_all.sh --aggregate-only   after editing sample_map.tsv

## Per sample
raw counts -> NanoFilt (q>=10, >=500 bp) -> HUMAN REMOVAL (minimap2 map-ont vs GRCh38, unmapped reads kept)
-> map to H99 -> coverage (mosdepth) -> read identity -> haploid variants (bcftools) -> virulence-gene
coverage/variants/consensus -> Flye + Medaka if >=1 Mb of reads map to H99.

## Outputs (outputs/summary/)
- sample_summary.tsv: reads, % human, mean depth, breadth >=1/3/5/10/20x, identity to H99, assembly stats
- gene_summary.tsv: sample x gene reads, depth, breadth, variants, % consensus called
- gene_variants.tsv: every variant in the six genes
- alignments/<GENE>.aln.fa: MAFFT alignment of H99 gene + each sample consensus (supplementary data)
- gene_alignment_identity.tsv: mismatches/indels per gene vs H99
- outputs/samples/<run>__<barcode>/nohuman.fastq.gz is what to deposit in SRA (never raw reads: human DNA)

Settings (quality, min depth, Medaka model, SLURM) are in config.sh. MEDAKA_MODEL must match the Dorado model used.
