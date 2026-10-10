#!/usr/bin/env python3
"""Locate the virulence genes in the H99 GFF3 and write a BED (for mosdepth) + a TSV (symbol, coords, strand, provenance).

Order of evidence for each gene symbol:
  1. manual file (symbol chrom start end strand; 1-based inclusive)      -> always wins
  2. GFF gene/Name/gene_synonym/locus_tag/old_locus_tag equals the symbol
  3. candidate H99 locus tag, ACCEPTED ONLY IF the annotated product matches expected keywords
Anything not resolved makes the script exit with code 2 (unless --allow-partial)."""
import argparse, gzip, re, sys
from urllib.parse import unquote

SYMBOLS = ["CAP59", "CNA1", "URE1", "PLB1", "LAC1", "LAC2"]
# Recollected H99 locus tags: used ONLY as a last resort and ONLY if the product keywords also match.
CANDIDATE_TAGS = {"CAP59": "CNAG_00721", "CNA1": "CNAG_04796", "URE1": "CNAG_05540",
                  "PLB1": "CNAG_06085", "LAC1": "CNAG_03465", "LAC2": "CNAG_05607"}
PRODUCT_KEYS = {"CAP59": r"capsul|cap59", "CNA1": r"calcineurin", "URE1": r"urease",
                "PLB1": r"phospholipase|lysophospholipase", "LAC1": r"laccase|multicopper|oxidase",
                "LAC2": r"laccase|multicopper|oxidase"}

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def parse_attrs(s):
    d = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = unquote(v)
    return d

def load_gff(path):
    genes, parent_of, products = {}, {}, {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, typ, s, e, _, strand, _, at = f[:9]
            a = parse_attrs(at)
            if typ == "gene":
                genes[a.get("ID", f"{chrom}:{s}-{e}")] = dict(chrom=chrom, start=int(s), end=int(e),
                                                              strand=strand, attrs=a)
            else:
                if "ID" in a and "Parent" in a:
                    parent_of[a["ID"]] = a["Parent"].split(",")[0]
                if "product" in a and "Parent" in a:
                    products.setdefault(a["Parent"].split(",")[0], set()).add(a["product"])
    gene_products = {g: set() for g in genes}
    for pid, prods in products.items():
        cur, hops = pid, 0
        while cur not in genes and cur in parent_of and hops < 5:
            cur, hops = parent_of[cur], hops + 1
        if cur in genes:
            gene_products[cur] |= prods
    return genes, gene_products

def names_of(g):
    a, out = g["attrs"], set()
    for k in ("gene", "Name", "locus_tag", "old_locus_tag"):
        if a.get(k):
            out.add(a[k].upper())
    for k in ("gene_synonym", "gene_synonyms"):
        for x in re.split(r"[,;]", a.get(k, "")):
            if x.strip():
                out.add(x.strip().upper())
    return out

def read_fai(p):
    d = {}
    with open(p) as fh:
        for line in fh:
            f = line.split("\t")
            d[f[0]] = int(f[1])
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True); ap.add_argument("--fai", required=True)
    ap.add_argument("--manual"); ap.add_argument("--out-bed", required=True); ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--allow-partial", action="store_true")
    a = ap.parse_args()

    fai = read_fai(a.fai)
    genes, gprod = load_gff(a.gff)
    manual = {}
    if a.manual:
        try:
            with open(a.manual) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 5 and f[0].lower() != "symbol" and not line.startswith("#"):
                        manual[f[0].upper()] = (f[1], int(f[2]), int(f[3]), f[4])
        except FileNotFoundError:
            pass

    rows, missing = [], []
    for sym in SYMBOLS:
        if sym in manual:
            c, s, e, st = manual[sym]
            rows.append((sym, c, s, e, st, "manual", "", ""))
            continue
        hits = [gid for gid, g in genes.items() if sym in names_of(g)]
        how = "symbol"
        if not hits:
            tag = CANDIDATE_TAGS.get(sym)
            hits = [gid for gid, g in genes.items() if tag and tag.upper() in names_of(g)]
            hits = [gid for gid in hits
                    if re.search(PRODUCT_KEYS[sym], " ".join(gprod.get(gid, [])) + " " + genes[gid]["attrs"].get("description", ""), re.I)]
            how = "locus_tag(CHECK product)"
        if len(hits) > 1:
            exact = [gid for gid in hits if genes[gid]["attrs"].get("gene", "").upper() == sym]
            hits = exact if len(exact) == 1 else hits
        if len(hits) != 1:
            missing.append((sym, len(hits)))
            continue
        g = genes[hits[0]]
        rows.append((sym, g["chrom"], g["start"], g["end"], g["strand"], how,
                     g["attrs"].get("locus_tag", g["attrs"].get("Name", "")),
                     "; ".join(sorted(gprod.get(hits[0], [])))[:120]))

    for r in rows:
        if r[1] not in fai:
            sys.exit(f"ERROR: {r[0]} is on '{r[1]}', which is not in the reference FASTA index. "
                     f"GFF and FASTA must come from the same assembly.")
        if not (1 <= r[2] <= r[3] <= fai[r[1]]) or r[4] not in "+-":
            sys.exit(f"ERROR: bad coordinates/strand for {r[0]}: {r}")

    with open(a.out_bed, "w") as bed, open(a.out_tsv, "w") as tsv:
        tsv.write("symbol\tchrom\tstart\tend\tstrand\tmatched_by\tlocus_tag\tproduct\n")
        for sym, c, s, e, st, how, lt, prod in rows:
            bed.write(f"{c}\t{s-1}\t{e}\t{sym}\n")
            tsv.write(f"{sym}\t{c}\t{s}\t{e}\t{st}\t{how}\t{lt}\t{prod}\n")

    for r in rows:
        print("resolved  %-6s %s:%d-%d (%s)  via %s  %s | %s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))
    if any("CHECK" in r[5] for r in rows):
        print("NOTE: genes matched by locus tag must be confirmed against the 'product' shown above.")
    if missing:
        for sym, n in missing:
            print(f"NOT RESOLVED {sym} ({n} candidate(s))", file=sys.stderr)
        print("Find them with e.g.  grep -i -w 'CAP59' <gff> | head   and put the coordinates in the manual file:\n"
              "  symbol<TAB>chrom<TAB>start<TAB>end<TAB>strand   (1-based, inclusive)", file=sys.stderr)
        if not a.allow_partial:
            sys.exit(2)

if __name__ == "__main__":
    main()
