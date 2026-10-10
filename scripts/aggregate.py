#!/usr/bin/env python3
"""Aggregate per-sample results.
  --stage tables : sample_summary.tsv, gene_summary.tsv, gene_variants.tsv, alignments/<GENE>.unaligned.fa
  --stage align  : reads alignments/<GENE>.aln.fa (MAFFT) -> gene_alignment_identity.tsv
Pure standard library."""
import argparse, csv, glob, gzip, os, re, sys

def kv(path):
    d = {}
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    d[p[0]] = p[1]
    return d

def rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if l.strip()]

def stats_sn(path):
    d = {}
    for r in rows(path):
        if r[0] == "SN" and len(r) >= 3:
            d[r[1].rstrip(":")] = r[2]
    return d

def read_fasta(path):
    seqs, name = {}, None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]; seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(COMP)[::-1]

def f(x, nd=4):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "NA"

def load_manifest(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def stage_tables(a):
    man = load_manifest(a.manifest)
    fai = {r[0]: int(r[1]) for r in rows(a.crypto_fa + ".fai")}
    genome_len = sum(fai.values())
    genes = {r[0]: r for r in rows(a.genes_tsv)[1:]}          # symbol -> row
    genome = read_fasta(a.crypto_fa)
    os.makedirs(os.path.join(a.summary, "alignments"), exist_ok=True)

    samp_cols = ["sample_id", "sample_key", "run", "barcode", "status",
                 "raw_reads", "raw_bases", "filtered_reads", "filtered_bases",
                 "human_reads_removed", "human_pct_of_filtered_reads", "nohuman_reads", "nohuman_bases",
                 "crypto_mapped_reads", "crypto_mapped_pct_of_nohuman_reads", "crypto_mapped_bases",
                 "genome_size_bp", "mean_depth_x", "breadth_ge1x_pct", "breadth_ge3x_pct", "breadth_ge5x_pct",
                 "breadth_ge10x_pct", "breadth_ge20x_pct",
                 "read_identity_mean", "read_identity_median", "reads_used_for_identity",
                 "alignment_error_rate", "snps_vs_H99", "indels_vs_H99", "bases_called_at_min_depth",
                 "consensus_identity_snp_only_pct", "consensus_identity_snp_indel_pct",
                 "assembly_status", "asm_contigs", "asm_total_size", "asm_avg_len", "asm_largest", "asm_smallest",
                 "asm_N50", "asm_identity_to_H99", "asm_aligned_bp"]
    gene_cols = ["sample_id", "sample_key", "gene", "reads_overlapping", "mean_depth_x", "bases_ge1x_pct",
                 "bases_ge3x_pct", "bases_ge5x_pct", "bases_ge10x_pct", "variants_called", "consensus_called_pct",
                 "in_alignment"]
    samp_out, gene_out, var_out = [], [], []
    unaligned = {g: [] for g in genes}
    for g in genes:                                            # reference sequence first
        r = genes[g]; c, s, e, st = r[1], int(r[2]), int(r[3]), r[4]
        seq = genome[c][s - 1:e]
        unaligned[g].append((f"REF_H99|{g}|{c}:{s}-{e}({st})", revcomp(seq) if st == "-" else seq))

    for m in man:
        key, sid = m["sample_key"], m["sample_id"]
        d = os.path.join(a.samples, key)
        mt = kv(os.path.join(d, "metrics.tsv"))
        status = mt.get("status", "not_run")
        row = dict.fromkeys(samp_cols, "NA")
        row.update(sample_id=sid, sample_key=key, run=m["run"], barcode=m["barcode"], status=status)
        for k in ("raw_reads", "raw_bases", "filtered_reads", "filtered_bases", "human_reads_removed",
                  "nohuman_reads", "nohuman_bases", "crypto_mapped_reads", "assembly_status"):
            if k in mt: row[k] = mt[k]
        if "filtered_reads" in mt and int(mt["filtered_reads"]) > 0 and "human_reads_removed" in mt:
            row["human_pct_of_filtered_reads"] = f(100 * int(mt["human_reads_removed"]) / int(mt["filtered_reads"]), 2)
        if "nohuman_reads" in mt and int(mt["nohuman_reads"]) > 0 and "crypto_mapped_reads" in mt:
            row["crypto_mapped_pct_of_nohuman_reads"] = f(100 * int(mt["crypto_mapped_reads"]) / int(mt["nohuman_reads"]), 3)
        sn = stats_sn(os.path.join(d, "crypto.stats.txt"))
        if sn:
            row["crypto_mapped_bases"] = sn.get("bases mapped (cigar)", "NA")
            row["alignment_error_rate"] = f(sn.get("error rate"), 5)
        row["genome_size_bp"] = str(genome_len)
        summ = rows(os.path.join(d, "cov.mosdepth.summary.txt"))
        for r in summ:
            if r and r[0] == "total":
                row["mean_depth_x"] = f(r[3], 4)
        for r in rows(os.path.join(d, "breadth.tsv")):
            row[f"breadth_ge{r[0]}x_pct"] = f(100 * float(r[1]), 3)
        ri = rows(os.path.join(d, "read_identity.tsv"))
        if ri:
            row["read_identity_mean"], row["read_identity_median"], row["reads_used_for_identity"] = ri[0][0], ri[0][1], ri[0][2]
        vs = rows(os.path.join(d, "variants_summary.tsv"))
        if vs:
            snp, ind, cov = int(vs[0][0]), int(vs[0][1]), int(vs[0][2])
            row.update(snps_vs_H99=snp, indels_vs_H99=ind, bases_called_at_min_depth=cov)
            if cov > 0:
                row["consensus_identity_snp_only_pct"] = f(100 * (1 - snp / cov), 3)
                row["consensus_identity_snp_indel_pct"] = f(100 * (1 - (snp + ind) / cov), 3)
        ast = rows(os.path.join(d, "assembly", "assembly_stats.tsv"))
        if len(ast) >= 2:
            h = dict(zip(ast[0], ast[1]))
            row.update(asm_contigs=h.get("num_seqs", "NA"), asm_total_size=h.get("sum_len", "NA"),
                       asm_avg_len=h.get("avg_len", "NA"), asm_largest=h.get("max_len", "NA"),
                       asm_smallest=h.get("min_len", "NA"), asm_N50=h.get("N50", "NA"))
        av = rows(os.path.join(d, "assembly", "asm_vs_ref.tsv"))
        if av:
            row["asm_identity_to_H99"], row["asm_aligned_bp"] = av[0][0], av[0][1]
        samp_out.append(row)

        # ---- genes
        greads = {r[0]: r for r in rows(os.path.join(d, "genes", "gene_reads.tsv"))}
        thr, reg = {}, {}
        p = os.path.join(d, "cov.thresholds.bed.gz")
        if os.path.exists(p):
            with gzip.open(p, "rt") as fh:
                hdr = None
                for line in fh:
                    x = line.rstrip("\n").split("\t")
                    if line.startswith("#"):
                        hdr = x; continue
                    thr[x[3]] = (int(x[2]) - int(x[1]), dict(zip(hdr[4:], map(int, x[4:]))))
        p = os.path.join(d, "cov.regions.bed.gz")
        if os.path.exists(p):
            with gzip.open(p, "rt") as fh:
                for line in fh:
                    x = line.rstrip("\n").split("\t")
                    reg[x[3]] = x[4]
        for g in genes:
            grow = dict.fromkeys(gene_cols, "NA")
            grow.update(sample_id=sid, sample_key=key, gene=g)
            if g in greads:
                grow["reads_overlapping"], grow["variants_called"] = greads[g][5], greads[g][6]
            if g in reg: grow["mean_depth_x"] = f(reg[g], 3)
            if g in thr:
                L, t = thr[g]
                for T in (1, 3, 5, 10):
                    grow[f"bases_ge{T}x_pct"] = f(100 * t.get(f"{T}X", 0) / L, 2)
            cp = os.path.join(d, "genes", f"{g}.consensus.fa")
            called, included = None, "no"
            if os.path.exists(cp):
                seqs = read_fasta(cp)
                if seqs:
                    name, seq = next(iter(seqs.items()))
                    seq = seq.upper()
                    called = sum(seq.count(b) for b in "ACGT") / max(len(seq), 1)
                    grow["consensus_called_pct"] = f(100 * called, 2)
                    nreads = int(grow["reads_overlapping"]) if grow["reads_overlapping"] != "NA" else 0
                    if nreads > 0 and called >= a.min_called:
                        included = "yes"
                        unaligned[g].append((f"{sid}|{m['run']}|{m['barcode']}", seq))
            grow["in_alignment"] = included
            gene_out.append(grow)
        gv = os.path.join(d, "genes", "gene_variants.tsv")
        for r in rows(gv):
            var_out.append([sid, key] + r)

    def write(path, cols, data):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(cols)
            for r in data:
                w.writerow([r.get(c, "NA") for c in cols] if isinstance(r, dict) else r)
    write(os.path.join(a.summary, "sample_summary.tsv"), samp_cols, samp_out)
    write(os.path.join(a.summary, "gene_summary.tsv"), gene_cols, gene_out)
    write(os.path.join(a.summary, "gene_variants.tsv"),
          ["sample_id", "sample_key", "gene", "chrom", "pos", "ref", "alt", "qual", "depth", "type"], var_out)
    for g, recs in unaligned.items():
        with open(os.path.join(a.summary, "alignments", f"{g}.unaligned.fa"), "w") as fh:
            for n, s in recs:
                fh.write(f">{n}\n")
                for i in range(0, len(s), 80):
                    fh.write(s[i:i + 80] + "\n")
        print(f"{g}: {len(recs) - 1} sample consensus sequence(s) + reference")
    print(f"wrote tables to {a.summary}")

def stage_align(a):
    out = os.path.join(a.summary, "gene_alignment_identity.tsv")
    with open(out, "w") as o:
        o.write("gene\tsample\tcolumns_compared\tmismatches\tindel_columns\tidentity_pct\tn_columns_masked\n")
        for p in sorted(glob.glob(os.path.join(a.summary, "alignments", "*.aln.fa"))):
            g = os.path.basename(p).replace(".aln.fa", "")
            seqs = read_fasta(p)
            names = list(seqs)
            if not names: continue
            ref = seqs[names[0]].upper()
            for n in names[1:]:
                s = seqs[n].upper()
                comp = mm = ind = msk = 0
                for r, x in zip(ref, s):
                    if r == "-" and x == "-": continue
                    if x == "N" or r == "N":
                        msk += 1; continue
                    if r == "-" or x == "-":
                        ind += 1; continue
                    comp += 1
                    if r != x: mm += 1
                idp = f(100 * (1 - mm / comp), 3) if comp else "NA"
                o.write(f"{g}\t{n.split('|')[0]}\t{comp}\t{mm}\t{ind}\t{idp}\t{msk}\n")
    print(f"wrote {out}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["tables", "align"])
    ap.add_argument("--manifest"); ap.add_argument("--samples"); ap.add_argument("--summary", required=True)
    ap.add_argument("--crypto-fa"); ap.add_argument("--genes-tsv"); ap.add_argument("--min-called", type=float, default=0.25)
    a = ap.parse_args()
    (stage_tables if a.stage == "tables" else stage_align)(a)
