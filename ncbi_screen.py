#!/usr/bin/env python3
"""Classify every NCBI Bdellovibrio genome against the curated 41-genome dataset, and re-estimate richness."""
import argparse, glob, json, os, re
from collections import defaultdict
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--ncbi-dir", required=True, help="ncbi_dataset/data")
ap.add_argument("--skani", required=True, help="skani triangle -E output")
ap.add_argument("--our-dir", required=True)
ap.add_argument("--metadata", required=True)
ap.add_argument("--prefix", default="ncbi_screen")
ap.add_argument("--min-compl", type=float, default=90)
ap.add_argument("--max-contam", type=float, default=5)
a = ap.parse_args()

base = lambda p: re.sub(r"\.(fna|fa|fasta|fas|fsa)(\.gz)?$", "", os.path.basename(p))

recs = {}
rep = os.path.join(a.ncbi_dir, "assembly_data_report.jsonl")
for line in open(rep):
    d = json.loads(line)
    acc = d.get("accession", "")
    ai = d.get("assemblyInfo", {}); bs = ai.get("biosample", {}) or {}
    attrs = {x.get("name", "").lower(): x.get("value", "") for x in bs.get("attributes", []) or []}
    org = d.get("organism", {}) or {}
    strain = (org.get("infraspecificNames", {}) or {}).get("strain", "") or attrs.get("strain", "") \
        or (org.get("infraspecificNames", {}) or {}).get("isolate", "") or attrs.get("isolate", "")
    warns = " ".join((ai.get("atypical", {}) or {}).get("warnings", []) or []).lower()
    pkg = str(bs.get("package", "")).lower()
    name = org.get("organismName", "")
    mag = ("metagenome" in warns or "metagenome_source" in attrs or "mimag" in pkg
           or "uncultured" in name.lower())
    ck = d.get("checkmInfo", {}) or {}
    st = d.get("assemblyStats", {}) or {}
    recs[acc] = dict(accession=acc, organism=name, strain=strain, assembly_level=ai.get("assemblyLevel", ""),
                     mag=bool(mag), atypical=warns, biosample=bs.get("accession", ""),
                     checkm_completeness=ck.get("completeness"), checkm_contamination=ck.get("contamination"),
                     size_mb=round(int(st.get("totalSequenceLength", 0) or 0) / 1e6, 3),
                     contigs=st.get("numberOfContigs"),
                     isolation_source=attrs.get("isolation_source", ""))
files = {}
for f in glob.glob(os.path.join(a.ncbi_dir, "*", "*.fna")):
    files[os.path.basename(os.path.dirname(f))] = f
missing = sorted(set(recs) - set(files))
print(f"NCBI: {len(recs)} assemblies in the report, {len(files)} FASTA files" +
      (f"; {len(missing)} without FASTA" if missing else ""))

meta = pd.read_csv(a.metadata, sep="\t", dtype=str).fillna("").set_index("genome")
gs_of = meta["genomospecies"].to_dict()
outg = set(meta.index[meta["outgroup"].astype(str).str.lower() == "yes"]) if "outgroup" in meta.columns else set()
ours = {base(f) for f in glob.glob(os.path.join(a.our_dir, "*")) if base(f) in gs_of}
print(f"Curated genomes found: {len(ours)} of {len(gs_of)}")

E = pd.read_csv(a.skani, sep="\t")
cols = {c.lower(): c for c in E.columns}
pick = lambda *ks: next(cols[k] for k in ks if k in cols)
rf, qf, an = pick("ref_file"), pick("query_file"), pick("ani")
afr, afq = pick("align_fraction_ref"), pick("align_fraction_query")
def gid(p):
    d = os.path.basename(os.path.dirname(p))
    return d if d in recs else base(p)
ani = defaultdict(dict)
for r in E.itertuples(index=False):
    x, y = gid(getattr(r, rf)), gid(getattr(r, qf))
    v, f1, f2 = float(getattr(r, an)), float(getattr(r, afr)), float(getattr(r, afq))
    for u, w in ((x, y), (y, x)):
        if v > ani[u].get(w, (0, 0, 0))[0]:
            ani[u][w] = (v, f1, f2)

def best(g, pool):
    hits = [(v[0], h) for h, v in ani[g].items() if h in pool]
    return max(hits) if hits else (0.0, "")

ingroup = {g for g in ours if g not in outg}
rows = []
for acc, r in recs.items():
    b_in, h_in = best(acc, ingroup); b_out, h_out = best(acc, outg & ours)
    full = ani[acc].get(h_in, (0, 0, 0))
    if b_in >= 99.9 and min(full[1], full[2]) >= 80:
        cls, gs = "same_strain_as_dataset", gs_of.get(h_in, "")
    elif b_out >= 95 and b_out >= b_in:
        cls, gs = "outgroup_like", "Out"
    elif b_in >= 95:
        cls, gs = "existing_GS", gs_of.get(h_in, "")
    else:
        cls, gs = "new_GS", ""
    r.update(cls=cls, genomospecies=gs, best_ingroup_hit=h_in, best_ingroup_ani=round(b_in, 2),
             best_outgroup_ani=round(b_out, 2))
    rows.append(r)
df = pd.DataFrame(rows)

def components(nodes, thr):
    nodes = list(nodes); parent = {n: n for n in nodes}
    def find(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]; n = parent[n]
        return n
    for n in nodes:
        for m, v in ani[n].items():
            if m in parent and v[0] >= thr:
                parent[find(n)] = find(m)
    comp = defaultdict(list)
    for n in nodes:
        comp[find(n)].append(n)
    return list(comp.values())

newmask = df.cls.isin(["existing_GS", "new_GS"])
df["strain_cluster"] = df.accession
df["representative"] = False
q = lambda r: (not r.mag, float(r.checkm_completeness or 0), -float(r.checkm_contamination or 100))
for comp in components(df.accession[newmask], 99.9):
    sub = df[df.accession.isin(comp)]
    repacc = max(sub.itertuples(), key=q).accession
    df.loc[df.accession.isin(comp), "strain_cluster"] = repacc
    df.loc[df.accession == repacc, "representative"] = True
df.loc[newmask & ~df.representative, "cls"] = "clone_of_new"
k = 0
for comp in sorted(components(df.accession[df.cls == "new_GS"], 95), key=len, reverse=True):
    k += 1
    df.loc[df.accession.isin(comp), "genomospecies"] = f"NEW{k:02d}"
for acc, rep in zip(df.accession, df.strain_cluster):
    if acc != rep and df.loc[df.accession == acc, "cls"].iloc[0] == "clone_of_new":
        df.loc[df.accession == acc, "genomospecies"] = df.loc[df.accession == rep, "genomospecies"].iloc[0]
df["high_quality"] = (pd.to_numeric(df.checkm_completeness, errors="coerce") >= a.min_compl) & \
                     (pd.to_numeric(df.checkm_contamination, errors="coerce") < a.max_contam)
gs_best = df[df.cls.isin(["new_GS", "clone_of_new"])].groupby("genomospecies").best_ingroup_ani.max()
df["needs_genus_check"] = df.genomospecies.str.startswith("NEW") & (df.genomospecies.map(gs_best).fillna(0) < 75)
order = ["accession", "organism", "strain", "cls", "genomospecies", "best_ingroup_hit", "best_ingroup_ani",
         "best_outgroup_ani", "mag", "high_quality", "checkm_completeness", "checkm_contamination",
         "assembly_level", "contigs", "size_mb", "strain_cluster", "representative", "needs_genus_check",
         "isolation_source", "biosample", "atypical"]
df = df[order].sort_values(["cls", "genomospecies", "accession"])
df.to_csv(f"{a.prefix}_genomes.tsv", sep="\t", index=False)

print("\nClassification of NCBI assemblies:")
print(pd.crosstab(df.cls, df.mag.map({True: "MAG", False: "isolate"})).to_string())
ex = df[df.cls == "existing_GS"]
if len(ex):
    print("\nNew strains in curated genomospecies (representatives):")
    print(ex[ex.representative].groupby("genomospecies").size().rename("new strains").to_string())
nw = df[(df.cls == "new_GS")]
print(f"\nNew genomospecies: {nw.genomospecies.nunique()} "
      f"({(~nw.mag).sum()} isolate and {nw.mag.sum()} MAG representatives; "
      f"{nw[nw.needs_genus_check].genomospecies.nunique()} genomospecies below 75% ANI to any curated "
      f"Bdellovibrio - check genus)")

def richness(counts):
    x = np.array([c for c in counts if c > 0]); N, S = int(x.sum()), len(x)
    f1, f2 = int((x == 1).sum()), int((x == 2).sum())
    f0 = (N - 1) / N * f1 * (f1 - 1) / (2 * (f2 + 1))
    var = (f1 * (f1 - 1) / (2 * (f2 + 1)) + f1 * (2 * f1 - 1) ** 2 / (4 * (f2 + 1) ** 2)
           + f1 ** 2 * f2 * (f1 - 1) ** 2 / (4 * (f2 + 1) ** 4))
    K = np.exp(1.96 * np.sqrt(np.log(1 + var / f0 ** 2))) if f0 > 0 else 1
    cov = 1 - f1 / N * ((N - 1) * f1 / ((N - 1) * f1 + 2 * f2)) if N > 1 else np.nan
    return dict(genomes=N, genomospecies=S, singletons=f1, doubletons=f2, chao1=round(S + f0, 1),
                chao1_lo=round(S + f0 / K, 1), chao1_hi=round(S + f0 * K, 1),
                percent_sampled=round(100 * S / (S + f0), 1), coverage=round(cov, 3))
base_counts = meta.loc[[g for g in meta.index if g not in outg], "genomospecies"].value_counts().to_dict()
def add(counts, sub):
    c = dict(counts)
    for gsx in sub.genomospecies:
        c[gsx] = c.get(gsx, 0) + 1
    return c
reps = df[df.representative & df.cls.isin(["existing_GS", "new_GS"])]
iso = reps[~reps.mag & reps.high_quality]
mags = reps[reps.mag & reps.high_quality]
scen = [("A: curated 38 genomes", base_counts),
        ("B: + high-quality isolates", add(base_counts, iso)),
        ("C: + high-quality isolates and MAGs", add(add(base_counts, iso), mags))]
R = pd.DataFrame([dict(scenario=n, **richness(c.values())) for n, c in scen])
R.to_csv(f"{a.prefix}_richness.tsv", sep="\t", index=False)
print("\nRichness (strain-level units; MAGs and isolates need CheckM completeness >= "
      f"{a.min_compl:g}% and contamination < {a.max_contam:g}%):")
print(R.to_string(index=False))
lowq = reps[~reps.high_quality]
if len(lowq):
    print(f"\nNot counted: {len(lowq)} new representatives below the quality threshold or without CheckM values")
print(f"\nTables: {a.prefix}_genomes.tsv, {a.prefix}_richness.tsv")
