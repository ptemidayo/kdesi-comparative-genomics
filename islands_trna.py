#!/usr/bin/env python3
"""Are genomic islands (IslandViewer 4) next to tRNA / tmRNA genes?

tRNA genes are the classic insertion sites of integrative elements (prophages, ICEs): the integrase
recombines the element into the 3' end of a tRNA gene, leaving direct repeats at both island ends.
For each island this script reports the nearest tRNA/tmRNA (inside the island, or within --flank bp
of either edge on the same contig), the nearest integrase/recombinase/transposase, and whether the
island touches a contig end. A permutation test (islands of the same length placed at random on the
same contigs, --perm times) gives the number of islands expected near a tRNA by chance.

Island rows are matched to the annotation by locus tag (IslandViewer was run on the Prokka GenBank
files), so contig coordinates come from the Prokka GFF, not from IslandViewer's concatenated sequence.

Usage:
  python islands_trna.py --islands kdesi.xls --gff ~/bdello_analysis/gff/Bdellovibrio_bacteriovorus_kdesi.gff
  python islands_trna.py --islands HD100.xls --gff ~/bdello_analysis/gff/Bdellovibrio_bacteriovorus_HD100.gff
"""
import argparse, re, sys
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--islands", required=True, help="IslandViewer 4 download (.xls/.xlsx/.tsv/.csv)")
ap.add_argument("--gff", required=True, help="Prokka GFF of the same genome")
ap.add_argument("--flank", type=int, default=1000, help="bp beyond an island edge counted as adjacent")
ap.add_argument("--perm", type=int, default=10000)
ap.add_argument("--out", default=None, help="output TSV (default: <islands>_trna.tsv)")
a = ap.parse_args()

def read_table(p):
    for f in (lambda: pd.read_excel(p), lambda: pd.read_csv(p, sep="\t"), lambda: pd.read_csv(p)):
        try:
            d = f()
            if "Island start" in d.columns:
                return d
        except Exception:
            pass
    sys.exit(f"could not read {p} (expected IslandViewer columns 'Island start', 'Locus')")

isl = read_table(a.islands)

# ---- GFF: features and contig lengths
feat, clen = [], {}
for line in open(a.gff):
    if line.startswith("##FASTA"):
        break
    if line.startswith("##sequence-region"):
        _, c, s, e = line.split()[:4]; clen[c] = int(e); continue
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] in ("gene", "region"):
        continue
    attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    feat.append(dict(contig=f[0], type=f[2], start=int(f[3]), end=int(f[4]), strand=f[6],
                     locus=attr.get("locus_tag", attr.get("ID", "")), product=attr.get("product", "")))
feat = pd.DataFrame(feat)
for c, g in feat.groupby("contig"):
    clen.setdefault(c, int(g.end.max()))
rna = feat[feat.type.isin(["tRNA", "tmRNA"])].copy()
rna["label"] = rna.type + " " + rna["product"].str.replace("tRNA-", "", regex=False)
MOB = re.compile(r"integrase|recombinase|transposase|resolvase|invertase|phage|IS\d|insertion element|XerC|XerD|XerH", re.I)
loc2 = feat.set_index("locus")
print(f"{a.gff}: {len(clen)} contigs, {int((feat.type == 'CDS').sum())} CDS, "
      f"{int((rna.type == 'tRNA').sum())} tRNA, {int((rna.type == 'tmRNA').sum())} tmRNA")

# ---- islands in contig coordinates (via locus tags)
rows = []
for (s0, e0), g in isl.groupby(["Island start", "Island end"], sort=True):
    tags = [t for t in g["Locus"].astype(str) if t in loc2.index]
    if not tags:
        print(f"WARNING: island {s0}-{e0}: no locus tag found in the GFF; skipped"); continue
    sub = loc2.loc[tags]
    contigs = sub.contig.unique()
    if len(contigs) > 1:
        print(f"WARNING: island {s0}-{e0} spans contigs {list(contigs)}; using the one with most genes")
    c = sub.contig.value_counts().index[0]; sub = sub[sub.contig == c]
    rows.append(dict(island=f"{s0}-{e0}", contig=c, start=int(sub.start.min()), end=int(sub.end.max()),
                     length=int(sub.end.max() - sub.start.min() + 1), genes=len(g),
                     hypothetical=int(g["Product"].astype(str).str.contains("hypothetical", case=False).sum()),
                     mobility_genes="; ".join(sorted(set(p for p in g["Product"].astype(str) if MOB.search(p))))))
I = pd.DataFrame(rows).sort_values(["contig", "start"])
# merge overlapping predictions (IslandViewer lists overlapping islands from different methods separately)
merged = []
for r in I.to_dict("records"):
    if merged and merged[-1]["contig"] == r["contig"] and r["start"] <= merged[-1]["end"]:
        m = merged[-1]
        m["island"] += "+" + r["island"]; m["end"] = max(m["end"], r["end"])
        m["length"] = m["end"] - m["start"] + 1
        m["genes"] = max(m["genes"], r["genes"]); m["hypothetical"] = max(m["hypothetical"], r["hypothetical"])
        m["mobility_genes"] = "; ".join(sorted(set(filter(None, (m["mobility_genes"] + "; " + r["mobility_genes"]).split("; ")))))
    else:
        merged.append(dict(r))
I = pd.DataFrame(merged)
print(f"{len(I)} islands after merging overlapping predictions")

def near(c, s, e, fl):
    r = rna[(rna.contig == c) & (rna.end >= s - fl) & (rna.start <= e + fl)]
    out = []
    for x in r.itertuples():
        if x.start >= s and x.end <= e:
            pos = "inside"
        elif x.end < s:
            pos = f"{s - x.end} bp upstream"
        elif x.start > e:
            pos = f"{x.start - e} bp downstream"
        else:
            pos = "overlaps edge"
        out.append(f"{x.label} ({pos})")
    return out

res = []
for r in I.itertuples():
    hits = near(r.contig, r.start, r.end, a.flank)
    edge = min(r.start - 1, clen[r.contig] - r.end)
    # nearest tRNA on the same contig at any distance
    rc = rna[rna.contig == r.contig]
    d = min([max(0, max(x.start - r.end, r.start - x.end)) for x in rc.itertuples()], default=None)
    res.append(dict(island=r.island, contig=r.contig, contig_length=clen[r.contig], start=r.start, end=r.end,
                    length=r.length, genes=r.genes, hypothetical=r.hypothetical,
                    tRNA_within_flank="; ".join(hits) if hits else "none",
                    nearest_tRNA_bp=d if d is not None else "no tRNA on contig",
                    distance_to_contig_end_bp=edge, mobility_genes=r.mobility_genes or "none"))
R = pd.DataFrame(res)
obs = int((R.tRNA_within_flank != "none").sum())

# ---- permutation: same lengths, random positions on the same contigs
rng = np.random.default_rng(1)
by_c = {c: g[["start", "end"]].values for c, g in rna.groupby("contig")}
cnt = np.zeros(a.perm, int)
for k in range(a.perm):
    n = 0
    for r in I.itertuples():
        L = clen[r.contig]
        if r.length >= L:
            s = 1
        else:
            s = rng.integers(1, L - r.length + 2)
        e = s + r.length - 1
        t = by_c.get(r.contig)
        if t is not None and ((t[:, 1] >= s - a.flank) & (t[:, 0] <= e + a.flank)).any():
            n += 1
    cnt[k] = n
p = (np.sum(cnt >= obs) + 1) / (a.perm + 1)

out = a.out or re.sub(r"\.[^.]+$", "", a.islands) + "_trna.tsv"
R.to_csv(out, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(R[["island", "contig", "length", "genes", "tRNA_within_flank", "nearest_tRNA_bp",
         "distance_to_contig_end_bp", "mobility_genes"]].to_string(index=False))
print(f"\nIslands with a tRNA/tmRNA inside or within {a.flank} bp: {obs} of {len(R)}")
print(f"Expected by chance (same lengths, random positions on the same contigs): mean {cnt.mean():.2f}, "
      f"95% range {np.percentile(cnt, 2.5):.0f}-{np.percentile(cnt, 97.5):.0f}; permutation P = {p:.4f}")
print(f"Islands within 1 kb of a contig end: {int((R.distance_to_contig_end_bp <= 1000).sum())} of {len(R)}")
print(f"Results: {out}")
