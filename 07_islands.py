#!/usr/bin/env python3
"""Step 07b — summarise IslandViewer 4 results (web server; not scriptable locally).

Upload (https://www.pathogenomics.sfu.ca/islandviewer/upload/) the GenBank files that
07_synteny_islands.sh stages in $OUT/islands/upload/ (same default parameters for both strains;
HD100/NC_005363.1 is also available pre-computed on the server). For each genome download
'All islands (tab-delimited)' and save as $OUT/islands/results/<genome_id>.tsv, then rerun.
"""
import argparse, glob, os, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()
wd = os.path.join(a.out, "islands")
files = sorted(glob.glob(os.path.join(wd, "results", "*.tsv")))
if not files:
    print(f"No IslandViewer results yet in {wd}/results — upload the files in {wd}/upload first.")
    raise SystemExit(0)

def genome_len(g):
    n = 0
    for line in open(os.path.join(a.out, "fna", f"{g}.fna")):
        if not line.startswith(">"):
            n += len(line.strip())
    return n

def merged_bp(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs + 1 if cs is not None else 0)

rows, genes = [], []
for f in files:
    g = os.path.basename(f)[:-4]
    d = pd.read_csv(f, sep="\t", dtype=str)
    d.columns = [re.sub(r"\s+", "_", c.strip().lower()) for c in d.columns]
    d["island_start"], d["island_end"] = d.island_start.astype(int), d.island_end.astype(int)
    L = genome_len(g)
    for meth, sub in list(d.groupby("method")) + [("All_methods(union)", d)]:
        isl = sub.drop_duplicates(["island_start", "island_end"])
        bp = merged_bp(list(zip(isl.island_start, isl.island_end)))
        rows.append(dict(genome=g, method=meth, n_islands=len(isl), island_bp=bp,
                         percent_genome=round(100 * bp / L, 2),
                         n_genes=sub.get("locus", sub.get("gene_id", pd.Series(dtype=str))).dropna().nunique()))
    d.insert(0, "genome", g); genes.append(d)
pd.DataFrame(rows).to_csv(os.path.join(wd, "island_summary.tsv"), sep="\t", index=False)
pd.concat(genes).to_csv(os.path.join(wd, "island_genes_all.tsv"), sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
