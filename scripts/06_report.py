#!/usr/bin/env python3
"""Report: which barcodes failed and why, and sequencing depth for those that worked.
Usage: python3 06_report.py [BASE_DIR]"""
import csv, os, sys
from collections import Counter, defaultdict

base = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Dr_Fred_Data")
S = os.path.join(base, "outputs", "summary")
rows = list(csv.DictReader(open(os.path.join(S, "sample_summary.tsv")), delimiter="\t"))
genes = list(csv.DictReader(open(os.path.join(S, "gene_summary.tsv")), delimiter="\t"))

def num(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d
def gi(x):
    v = num(x); return int(v) if v is not None else 0

REASON = {
 "empty_input": "Fastq file contained no reads (nothing sequenced/demultiplexed to this barcode)",
 "no_reads_after_filter": "Reads present but ALL failed NanoFilt (mean quality < Q or length < min length)",
 "no_nonhuman_reads": "All reads that passed QC were human; none remained after human-read removal",
 "not_run": "Sample did not run (check SLURM log)",
}
def tier(d):
    if d is None: return "no data"
    if d >= 10: return ">=10x"
    if d >= 1: return "1-10x"
    if d >= 0.1: return "0.1-1x"
    return "<0.1x"

failed, worked = [], []
for r in rows:
    st = r["status"]
    if st == "ok" and gi(r["crypto_mapped_reads"]) > 0:
        worked.append(r)
    else:
        if st == "ok":
            r = dict(r); r["status"] = "no_reads_mapped_to_H99"
            REASON["no_reads_mapped_to_H99"] = "Non-human reads present but none mapped to the C. neoformans H99 genome"
        failed.append(r)

by_reason = Counter(r["status"] for r in failed)
genes_by = defaultdict(list)
for g in genes: genes_by[g["sample_key"]].append(g)

with open(os.path.join(S, "report_failed.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["sample_id", "sample_key", "run", "barcode", "raw_reads", "raw_bases", "filtered_reads",
                "human_reads_removed", "nohuman_reads", "reason_code", "explanation"])
    for r in sorted(failed, key=lambda r: (r["status"], r["run"], r["barcode"])):
        w.writerow([r["sample_id"], r["sample_key"], r["run"], r["barcode"], r["raw_reads"], r["raw_bases"],
                    r["filtered_reads"], r["human_reads_removed"], r["nohuman_reads"], r["status"], REASON.get(r["status"], "")])

cols = ["sample_id", "sample_key", "run", "raw_reads", "filtered_reads", "human_pct_of_filtered_reads", "nohuman_reads",
        "crypto_mapped_reads", "crypto_mapped_bases", "mean_depth_x", "breadth_ge1x_pct", "breadth_ge10x_pct",
        "read_identity_mean", "consensus_identity_snp_only_pct", "depth_tier", "genes_ge3x_of_6", "assembly_status",
        "asm_contigs", "asm_total_size", "asm_N50"]
with open(os.path.join(S, "report_worked.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(cols)
    for r in sorted(worked, key=lambda r: -(num(r["mean_depth_x"], 0))):
        gg = genes_by.get(r["sample_key"], [])
        n3 = sum(1 for g in gg if (num(g["bases_ge3x_pct"], 0) >= 50))
        rr = dict(r); rr["depth_tier"] = tier(num(r["mean_depth_x"])); rr["genes_ge3x_of_6"] = f"{n3}/{len(gg)}"
        w.writerow([rr.get(c, "NA") for c in cols])

runs = defaultdict(Counter)
for r in rows:
    k = "worked" if r in worked else r["status"] if r["status"] != "ok" else "no_reads_mapped_to_H99"
    runs[r["run"]][k] += 1
cats = ["worked", "empty_input", "no_reads_after_filter", "no_nonhuman_reads", "no_reads_mapped_to_H99"]
with open(os.path.join(S, "report_by_run.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["run", "total_files"] + cats)
    for run in sorted(runs): w.writerow([run, sum(runs[run].values())] + [runs[run][c] for c in cats])

tot_raw = sum(gi(r["raw_reads"]) for r in rows)
L = []
L.append(f"Barcode files analysed: {len(rows)} from {len(runs)} sequencing runs; total raw reads {tot_raw:,}")
L.append(f"Worked (reads mapped to C. neoformans H99): {len(worked)}    Did not work: {len(failed)}")
L.append("\nWhy barcodes did not yield C. neoformans data:")
for k, n in by_reason.most_common(): L.append(f"  {n:4d}  {k}: {REASON.get(k, '')}")
hp = [num(r["human_pct_of_filtered_reads"]) for r in worked if num(r["human_pct_of_filtered_reads"]) is not None]
L.append("\nSamples that worked, by mean genome depth (human reads removed):")
tc = Counter(tier(num(r["mean_depth_x"])) for r in worked)
for t in (">=10x", "1-10x", "0.1-1x", "<0.1x"): L.append(f"  {t:7s} {tc.get(t, 0)}")
if hp: L.append(f"Human share of filtered reads in these samples: median {sorted(hp)[len(hp)//2]:.1f}% (range {min(hp):.1f}-{max(hp):.1f}%)")
L.append("\nTop samples (by mean depth of the ~19 Mb H99 genome):")
L.append("  %-22s %8s %9s %8s %9s %9s %8s" % ("sample", "rawreads", "%human", "H99reads", "depth(x)", "breadth1x", "identity"))
for r in sorted(worked, key=lambda r: -(num(r["mean_depth_x"], 0)))[:15]:
    L.append("  %-22s %8s %9s %8s %9s %9s %8s" % (r["sample_id"][:22], r["raw_reads"], r["human_pct_of_filtered_reads"],
             r["crypto_mapped_reads"], r["mean_depth_x"], r["breadth_ge1x_pct"], r["read_identity_mean"]))
L.append("\nFiles: report_failed.tsv, report_worked.tsv, report_by_run.tsv in " + S)
open(os.path.join(S, "report.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
