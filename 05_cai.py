#!/usr/bin/env python3
"""Step 05 — Codon Adaptation Index (Sharp & Li 1987) and cross-strain translational profiles.

* reference set = each genome's own ribosomal-protein genes (Prokka product names,
  Pfam 'Ribosomal_*' domains as fallback)
* relative adaptiveness w = count / max count within each synonymous family
  (0.5 pseudocount for unobserved codons; Met, Trp and stop codons excluded; start codon skipped)
* CAI per CDS = geometric mean of w
* strain x category matrix of median CAI -> Euclidean distance -> Ward clustering
* concordance with ANI genomospecies: adjusted Rand index (+ permutation p-value)
"""
import argparse, os, re, sys, math
from collections import Counter, defaultdict
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from common import read_fasta, genome_ids, load_metadata, single_case_test

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--categories", default="",
                help="comma-separated categories for the matrix; default = all Pfam categories "
                     "with >= --min-genes genes in every genome")
ap.add_argument("--min-genes", type=int, default=3)
ap.add_argument("--min-ref", type=int, default=20)
ap.add_argument("--permutations", type=int, default=9999)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--foreground", default="")
ap.add_argument("--include-outgroups", action="store_true",
                help="use outgroup genomes in the matrix, clustering and ARI (default: excluded)")
a = ap.parse_args()
res = os.path.join(a.out, "cai"); os.makedirs(res, exist_ok=True)
ids = genome_ids(a.out)

# genetic code 11 (codon assignments identical to the standard code)
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a1 + a2 + a3: aas[16 * i + 4 * j + k]
        for i, a1 in enumerate(bases) for j, a2 in enumerate(bases) for k, a3 in enumerate(bases)}
FAMILIES = defaultdict(list)
for c, aa in CODE.items():
    if aa not in "*MW":
        FAMILIES[aa].append(c)

RIBO = re.compile(r"\b(30S|50S)\s+ribosomal\s+(sub[- ]?unit\s+)?protein\b|\bribosomal\s+protein\s+[LS]\d", re.I)
NOT_RIBO = re.compile(r"methyltransferase|acetyltransferase|modif|maturation|hibernation|silenc|"
                      r"pseudouridine|processing|assembly|kinase|synthase|deacylase|binding protein|"
                      r"glutamate|ligase|isomerase|inactivat", re.I)

def codons(seq):
    seq = seq.upper()
    return [seq[i:i + 3] for i in range(3, len(seq) - len(seq) % 3, 3)]   # skip start codon

def w_table(ref_seqs):
    cnt = Counter(c for s in ref_seqs for c in codons(s) if c in CODE)
    w = {}
    for aa, fam in FAMILIES.items():
        if len(fam) == 1:
            continue                                                        # non-synonymous
        x = {c: cnt.get(c, 0) or 0.5 for c in fam}
        m = max(x.values())
        for c in fam:
            w[c] = x[c] / m
    return w

def cai(seq, w):
    ls = [math.log(w[c]) for c in codons(seq) if c in w]
    return math.exp(sum(ls) / len(ls)) if ls else np.nan

pc_f = os.path.join(a.out, "pfam", "results", "protein_categories.tsv")
dom_f = os.path.join(a.out, "pfam", "results", "protein_domains_filtered.tsv")
pc = pd.read_csv(pc_f, sep="\t") if os.path.exists(pc_f) else None
dom = pd.read_csv(dom_f, sep="\t", usecols=["genome", "protein", "pfam_name"]) if os.path.exists(dom_f) else None

all_rows, ref_info = [], []
for g in ids:
    seqs, heads = read_fasta(os.path.join(a.out, "ffn", f"{g}.ffn"))
    ref = [k for k, h in heads.items() if RIBO.search(h) and not NOT_RIBO.search(h)]
    src = "annotation"
    if len(ref) < a.min_ref and dom is not None:
        ref = sorted(set(ref) | set(dom[(dom.genome == g) & dom.pfam_name.str.match(r"Ribosomal_[LS]")].protein))
        src = "annotation+Pfam"
    if len(ref) < a.min_ref:
        print(f"WARNING {g}: only {len(ref)} ribosomal reference genes")
    w = w_table([seqs[k] for k in ref])
    refset = set(ref)
    for k, s in seqs.items():
        if len(s) >= 300 and len(s) % 3 == 0:          # >= 100 codons, complete CDS
            all_rows.append((g, k, len(s) // 3, cai(s, w), k in refset))
    ref_info.append((g, len(ref), src))

cdf = pd.DataFrame(all_rows, columns=["genome", "gene", "codons", "CAI", "is_reference"])
cdf.to_csv(os.path.join(res, "cai_per_gene.tsv"), sep="\t", index=False)
pd.DataFrame(ref_info, columns=["genome", "n_ribosomal_reference", "source"])\
  .to_csv(os.path.join(res, "reference_sets.tsv"), sep="\t", index=False)
cdf.groupby("genome").CAI.describe().round(4).to_csv(os.path.join(res, "cai_genome_summary.tsv"), sep="\t")

# ---------- strain x category median CAI ----------
if pc is None:
    sys.exit("protein_categories.tsv missing — run step 03 first")
m = cdf.merge(pc[["genome", "protein", "category"]], left_on=["genome", "gene"], right_on=["genome", "protein"])
m = m[~m.is_reference]          # keep the reference genes out of the category medians
meta, ids = load_metadata(a.out, ids, exclude_outgroups=not a.include_outgroups)   # ids -> ingroup
if meta is None:
    sys.exit("metadata.tsv missing — run step 01 first")
meta = meta.loc[ids]
print(f"Translational-profile clustering on {len(ids)} genomes")
n = m.groupby(["genome", "category"]).size().unstack(fill_value=0).reindex(ids).fillna(0)
cats = [c.strip() for c in a.categories.split(",") if c.strip()] or \
       [c for c in n.columns if (n[c] >= a.min_genes).all()]
print(f"CAI matrix: {len(ids)} strains x {len(cats)} categories: {', '.join(cats)}")
M = m.groupby(["genome", "category"]).CAI.median().unstack().reindex(index=ids, columns=cats)
if M.isna().any().any():
    print("WARNING: missing cells filled with category median:", int(M.isna().sum().sum()))
    M = M.fillna(M.median())
M.round(4).to_csv(os.path.join(res, "median_cai_strain_x_category.tsv"), sep="\t")
if a.foreground in M.index:
    others = [g for g in ids if g != a.foreground]
    gm = cdf.groupby("genome").CAI.median().rename("genome_median_CAI")
    F = M.join(gm)
    sct = single_case_test(F.loc[a.foreground], F.loc[others])
    sct.to_csv(os.path.join(res, "foreground_vs_others_crawford_howell.tsv"), sep="\t", index=False)
    print(f"{a.foreground} median CAI vs {len(others)} others: "
          f"{(sct.p_bonferroni < 0.05).sum()} of {len(sct)} features differ (Bonferroni < 0.05)")
    d = pd.DataFrame(squareform(pdist(M.values)), index=ids, columns=ids)[a.foreground].drop(a.foreground)
    print(f"  nearest translational profile: {d.idxmin()} (distance {d.min():.4f})")

D = pdist(M.values, metric="euclidean")
pd.DataFrame(squareform(D), index=ids, columns=ids).round(5)\
  .to_csv(os.path.join(res, "euclidean_distance.tsv"), sep="\t")
Z = linkage(D, method="ward")

# ---------- ARI vs ANI genomospecies ----------
truth = meta.genomospecies.values
k_true = len(set(truth))
rng = np.random.default_rng(a.seed)
rows = []
for k in range(2, len(ids)):                       # every cut of the Ward tree
    lab = fcluster(Z, k, criterion="maxclust")
    r = dict(k=k, ARI=adjusted_rand_score(truth, lab), k_equals_n_genomospecies=(k == k_true), perm_p=np.nan)
    if k == k_true:                                 # permutation test at the genomospecies number
        perm = np.array([adjusted_rand_score(rng.permutation(truth), lab) for _ in range(a.permutations)])
        r["perm_p"] = (np.sum(perm >= r["ARI"]) + 1) / (a.permutations + 1)
    rows.append(r)
ari = pd.DataFrame(rows)
best = ari.loc[ari.ARI.idxmax()]
print(f"Best ARI over all cuts: {best.ARI:.3f} at k = {int(best.k)}")
ari.round(4).to_csv(os.path.join(res, "ARI_vs_genomospecies.tsv"), sep="\t", index=False)
pd.DataFrame({"genome": ids, "genomospecies": truth,
              f"cai_cluster_k{k_true}": fcluster(Z, k_true, criterion="maxclust")})\
  .to_csv(os.path.join(res, "cai_clusters.tsv"), sep="\t", index=False)
main = ari[ari.k == k_true]
if len(main):
    print(f"ARI (k = {k_true} genomospecies) = {main.ARI.iloc[0]:.3f}, permutation p = {main.perm_p.iloc[0]:.4f}")

# ---------- figure: dendrogram + heatmap ----------
fig = plt.figure(figsize=(12, max(6, len(ids) * 0.22)))
ax1 = fig.add_axes([0.02, 0.1, 0.18, 0.8])
dn = dendrogram(Z, orientation="left", labels=ids, ax=ax1, no_labels=True, color_threshold=0)
ax1.axis("off")
order = dn["ivl"]
ax2 = fig.add_axes([0.21, 0.1, 0.62, 0.8])
H = M.loc[order[::-1]]
im = ax2.imshow(H.values, aspect="auto", cmap="viridis")
lbl = meta["strain"] if "strain" in meta.columns else pd.Series(H.index, index=H.index)
ax2.set_yticks(range(len(H))); ax2.set_yticklabels([f"{lbl[g]} [{meta.loc[g,'genomospecies']}]" for g in H.index], fontsize=7)
ax2.yaxis.tick_right()
ax2.set_xticks(range(len(cats))); ax2.set_xticklabels(cats, rotation=60, ha="right", fontsize=7)
cax = fig.add_axes([0.93, 0.1, 0.015, 0.3]); fig.colorbar(im, cax=cax, label="median CAI")
fig.savefig(os.path.join(res, "cai_ward_heatmap.png"), dpi=300, bbox_inches="tight")
print(f"Results in {res}")
