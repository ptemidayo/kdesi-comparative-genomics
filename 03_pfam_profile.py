#!/usr/bin/env python3
"""Step 03 — Proteome-wide Pfam domain profiling.

* filters hmmscan domains (i-Evalue <= 1e-5, >= 30% HMM coverage)
* maps Pfam accessions to functional categories (pfam_categories.tsv)
* % of proteome per category (proteins with >= 1 category domain / total proteins x 100)
* MCP sub-profile (the Pfam entries flagged mcp=yes)
* z-scores across genomes, Mann–Whitney U between metadata groups (Bonferroni)
* PCA on standardized %-proteome values (scikit-learn)
"""
import argparse, glob, os, itertools, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from common import read_fasta, parse_domtbl, genome_ids, load_metadata, single_case_test, single_case_regression, bh_fdr
from math import comb

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True, help="pipeline OUT directory")
ap.add_argument("--map", required=True)
ap.add_argument("--evalue", type=float, default=1e-5)
ap.add_argument("--min-cov", type=float, default=0.30)
ap.add_argument("--foreground", default="kdesi")
ap.add_argument("--min-group-size", type=int, default=3)
ap.add_argument("--include-outgroups", action="store_true",
                help="use outgroup genomes in z-scores, Mann–Whitney and PCA (default: excluded)")
a = ap.parse_args()

res = os.path.join(a.out, "pfam", "results"); os.makedirs(res, exist_ok=True)
ids = genome_ids(a.out)
meta, stat_ids = load_metadata(a.out, ids, exclude_outgroups=not a.include_outgroups)
print(f"{len(ids)} genomes; {len(stat_ids)} used for z-scores / Mann–Whitney / PCA"
      + ("" if len(stat_ids) == len(ids) else f" (excluded: {', '.join(g for g in ids if g not in stat_ids)})"))
pmap = pd.read_csv(a.map, sep="\t", dtype=str).fillna("")
pmap["pfam_acc"] = pmap["pfam_acc"].str.split(".").str[0]
acc2cat = dict(zip(pmap.pfam_acc, pmap.category))
categories = list(dict.fromkeys(pmap.category))
mcp_accs = list(pmap.loc[pmap.mcp.str.lower() == "yes", "pfam_acc"])
print(f"Map: {len(pmap)} Pfam accessions, {len(categories)} categories, {len(mcp_accs)} MCP entries")

# ---------------- parse domains ----------------
dom_rows, nprot = [], {}
for g in ids:
    seqs, _ = read_fasta(os.path.join(a.out, "faa", f"{g}.faa"))
    nprot[g] = len(seqs)
    for prot, acc, name, ev, cov in parse_domtbl(os.path.join(a.out, "pfam", "domtbl", f"{g}.domtbl"),
                                                   a.evalue, a.min_cov):
        dom_rows.append((g, prot, acc, name, ev, round(cov, 3)))
dom = pd.DataFrame(dom_rows, columns=["genome", "protein", "pfam_acc", "pfam_name", "i_evalue", "hmm_cov"])
dom.to_csv(os.path.join(res, "protein_domains_filtered.tsv"), sep="\t", index=False)
pd.Series(nprot, name="n_proteins").rename_axis("genome").to_csv(os.path.join(res, "proteome_sizes.tsv"), sep="\t")

# Sanity check the mapping table against what hmmscan actually reports
obs_name = dom.drop_duplicates("pfam_acc").set_index("pfam_acc")["pfam_name"]
chk = pmap[["pfam_acc", "pfam_name", "category"]].copy()
chk["observed_name"] = chk.pfam_acc.map(obs_name).fillna("")
chk["total_domain_hits"] = chk.pfam_acc.map(dom.pfam_acc.value_counts()).fillna(0).astype(int)
chk["flag"] = np.where(chk.total_domain_hits == 0, "NO_HITS",
               np.where(chk.observed_name.str.lower() != chk.pfam_name.str.lower(), "NAME_DIFFERS", "ok"))
chk.to_csv(os.path.join(res, "pfam_map_check.tsv"), sep="\t", index=False)
bad = chk[chk.flag != "ok"]
if len(bad):
    print(f"WARNING: {len(bad)} map entries flagged (see pfam_map_check.tsv):")
    print(bad.to_string(index=False))

# ---------------- all-Pfam count matrix ----------------
allpf = (dom.drop_duplicates(["genome", "protein", "pfam_acc"])
            .groupby(["genome", "pfam_acc"]).size().unstack(fill_value=0).reindex(ids).fillna(0).astype(int))
allpf.to_csv(os.path.join(res, "pfam_all_protein_counts.tsv"), sep="\t")

# ---------------- category assignment ----------------
d = dom[dom.pfam_acc.isin(acc2cat)].copy()
d["category"] = d.pfam_acc.map(acc2cat)
pc = (d.groupby(["genome", "protein", "category"])["pfam_acc"]
        .apply(lambda s: ",".join(sorted(set(s)))).reset_index())
selmap = dict(zip(pmap.category, pmap.selection_class))
pc["selection_class"] = pc.category.map(selmap)
pc.to_csv(os.path.join(res, "protein_categories.tsv"), sep="\t", index=False)

counts = (pc.groupby(["genome", "category"]).size().unstack(fill_value=0)
            .reindex(index=ids, columns=categories).fillna(0).astype(int))
pct = counts.div(pd.Series(nprot), axis=0) * 100
counts.to_csv(os.path.join(res, "category_protein_counts.tsv"), sep="\t")
pct.round(4).to_csv(os.path.join(res, "category_percent_proteome.tsv"), sep="\t")

# ---------------- MCP profile ----------------
mcp = allpf.reindex(columns=mcp_accs, fill_value=0)
mcp.columns = [f"{c}_{pmap.set_index('pfam_acc').loc[c, 'pfam_name']}" for c in mcp.columns]
mcp_prot = dom[dom.pfam_acc == "PF00015"].groupby("genome").protein.nunique().reindex(ids).fillna(0)
mcp["MCP_proteins(PF00015)"] = mcp_prot.astype(int)
mcp["MCP_percent_proteome"] = (mcp_prot / pd.Series(nprot) * 100).round(4)
mcp.to_csv(os.path.join(res, "mcp_domain_profile.tsv"), sep="\t")

# ---------------- z-scores ----------------
ref = pct.loc[stat_ids]                      # reference = ingroup genomes
z = (pct - ref.mean()) / ref.std(ddof=1)
z.round(3).to_csv(os.path.join(res, "category_zscores.tsv"), sep="\t")

# ---------------- foreground vs all other ingroup genomes (single-case test) ----------------
if a.foreground in pct.index:
    others = [g for g in stat_ids if g != a.foreground]
    feats_fg = pd.concat([pct, mcp[["MCP_percent_proteome"]]], axis=1)
    sc_t = single_case_test(feats_fg.loc[a.foreground], feats_fg.loc[others])
    sc_t.to_csv(os.path.join(res, "foreground_vs_others_crawford_howell.tsv"), sep="\t", index=False)
    print(f"\n{a.foreground} vs {len(others)} other ingroup genomes (Crawford–Howell; Bonferroni / BH):")
    print(f"  {'feature':28s} {'value':>7s} {'ref mean':>8s} {'z':>6s} {'pctl':>5s} {'p':>8s} {'p_bonf':>7s}")
    for r in sc_t.itertuples():
        flag = "**" if r.p_bonferroni < 0.05 else ("*" if r.q_BH < 0.05 else "")
        print(f"  {r.feature:28s} {r.foreground:7.3f} {r.ref_mean:8.3f} {r.z:6.2f} {r.percentile:5.0f} "
              f"{r.p:8.4f} {r.p_bonferroni:7.4f} {flag}")
    print("  ** Bonferroni < 0.05   * BH q < 0.05 only")

    # ---------- size-adjusted test: gene counts vs proteome size ----------
    size = pd.Series(nprot).reindex(ids).astype(float)
    cnt = counts.copy()
    cnt["MCP_proteins"] = mcp_prot.reindex(ids)
    rows = []
    for c in cnt.columns:
        r = single_case_regression(size[a.foreground], float(cnt.loc[a.foreground, c]),
                                   size[others].values, cnt.loc[others, c].values)
        rows.append(dict(feature=c, **r))
    adj = pd.DataFrame(rows)
    adj["p_bonferroni"] = np.minimum(adj.p * adj.p.notna().sum(), 1.0)
    adj["q_BH"] = bh_fdr(adj.p)
    adj.to_csv(os.path.join(res, "foreground_vs_others_size_adjusted.tsv"), sep="\t", index=False)
    print(f"\nSize-adjusted ({a.foreground} gene count vs count expected for its proteome of "
          f"{int(size[a.foreground])} proteins; regression on {len(others)} genomes, Crawford & Garthwaite):")
    print(f"  {'feature':28s} {'obs':>5s} {'exp':>7s} {'resid':>7s} {'t':>6s} {'p':>8s} {'p_bonf':>7s}")
    for r in adj.itertuples():
        flag = "**" if r.p_bonferroni < 0.05 else ("*" if r.q_BH < 0.05 else "")
        print(f"  {r.feature:28s} {r.observed:5.0f} {r.expected:7.1f} {r.residual:+7.1f} {r.t:6.2f} "
              f"{r.p:8.4f} {r.p_bonferroni:7.4f} {flag}")

    # family level (exploratory): every mapped Pfam family present in >= 5 genomes, BH-corrected
    fam = allpf.reindex(columns=[x for x in pmap.pfam_acc if x in allpf.columns]).fillna(0)
    fam = fam.loc[:, (fam.loc[stat_ids] > 0).sum() >= 5]
    frows = []
    for acc in fam.columns:
        r = single_case_regression(size[a.foreground], float(fam.loc[a.foreground, acc]),
                                   size[others].values, fam.loc[others, acc].values)
        frows.append(dict(pfam_acc=acc, pfam_name=pmap.set_index("pfam_acc").loc[acc, "pfam_name"],
                          category=acc2cat.get(acc, ""), ref_median=fam.loc[others, acc].median(),
                          ref_max=fam.loc[others, acc].max(), **r))
    fadj = pd.DataFrame(frows)
    fadj["q_BH"] = bh_fdr(fadj.p)
    fadj.sort_values("p").to_csv(os.path.join(res, "foreground_pfam_families_size_adjusted.tsv"), sep="\t", index=False)
    top = fadj.sort_values("p").head(8)
    print(f"\nPfam families, size-adjusted (exploratory; {len(fadj)} families, BH):")
    for r in top.itertuples():
        print(f"  {r.pfam_acc} {r.pfam_name:22s} obs {r.observed:3.0f}  exp {r.expected:5.1f}  "
              f"median {r.ref_median:4.1f}  max {r.ref_max:3.0f}  p={r.p:.4f}  q={r.q_BH:.3f}")

# ---------------- Mann–Whitney U between groups (Bonferroni) ----------------
if meta is not None:
    sm = meta.loc[stat_ids]
    feats = pd.concat([pct, mcp[["MCP_percent_proteome"]]], axis=1)
    gcol = sm["group"].fillna("").astype(str).str.strip() if "group" in sm.columns else pd.Series("", index=sm.index)
    sm = sm.assign(group=gcol)
    groups = [g for g, n in sm.group[sm.group != ""].value_counts().items() if n >= a.min_group_size]
    rows = []
    for g1, g2 in itertools.combinations(sorted(groups), 2):
        s1, s2 = sm.index[sm.group == g1], sm.index[sm.group == g2]
        for f in feats.columns:
            x, y = feats.loc[s1, f], feats.loc[s2, f]
            if pd.concat([x, y]).nunique() == 1:          # identical values -> no test
                U, p = np.nan, 1.0
            else:
                U, p = mannwhitneyu(x, y, alternative="two-sided")
            rbc = 1 - 2 * U / (len(x) * len(y)) if not np.isnan(U) else np.nan
            rows.append(dict(feature=f, group1=g1, group2=g2, n1=len(x), n2=len(y),
                             median1=x.median(), median2=y.median(), U=U, rank_biserial=rbc, p=p,
                             min_attainable_p=min(1.0, 2 / comb(len(x) + len(y), len(x)))))
    mw = pd.DataFrame(rows)
    if len(mw):
        mw["p_bonferroni"] = np.minimum(mw.p * len(mw), 1.0)
        mw["significant"] = mw.p_bonferroni < 0.05
        mw["can_reach_significance"] = mw.min_attainable_p * len(mw) < 0.05
        mw.to_csv(os.path.join(res, "mannwhitney_bonferroni.tsv"), sep="\t", index=False)
        print(f"\nMann–Whitney: {len(mw)} tests across {len(groups)} groups "
              f"({', '.join(f'{g} n={(sm.group == g).sum()}' for g in sorted(groups))}); "
              f"{mw.significant.sum()} significant after Bonferroni")
        if not mw.can_reach_significance.any():
            print("  WARNING: with these group sizes and this many tests, NO comparison can reach "
                  "Bonferroni p < 0.05 (smallest possible p is shown in min_attainable_p). "
                  "Use larger groups in the 'group' column of metadata.tsv.")
    else:
        print(f"\nBetween-group Mann–Whitney skipped (optional): fill the 'group' column of the metadata "
              f"table (e.g. habitat) with >= 2 groups of >= {a.min_group_size} genomes")
else:
    print("metadata.tsv not found — Mann–Whitney skipped")

# ---------------- PCA ----------------
scaler = StandardScaler().fit(pct.loc[stat_ids].fillna(0))          # fit on ingroup only
X = scaler.transform(pct.loc[stat_ids].fillna(0))
pca = PCA(n_components=min(5, X.shape[1], X.shape[0])).fit(X)
sc = pd.DataFrame(pca.transform(X), index=stat_ids, columns=[f"PC{i+1}" for i in range(pca.n_components_)])
sc.round(4).to_csv(os.path.join(res, "pca_scores.tsv"), sep="\t")
pd.DataFrame(pca.components_.T, index=pct.columns, columns=sc.columns).round(4)\
  .to_csv(os.path.join(res, "pca_loadings.tsv"), sep="\t")
pd.Series(pca.explained_variance_ratio_, index=sc.columns, name="explained_variance_ratio")\
  .to_csv(os.path.join(res, "pca_variance.tsv"), sep="\t")

fig, ax = plt.subplots(figsize=(9, 6.5))
lab = meta.loc[stat_ids, "genomospecies"] if meta is not None else pd.Series("all", index=stat_ids)
short = meta.loc[stat_ids, "strain"] if meta is not None and "strain" in meta.columns else pd.Series(stat_ids, index=stat_ids)
colors = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors)
markers = "osD^vP*X"
for i, gs in enumerate(sorted(lab.unique())):
    m = lab == gs
    ax.scatter(sc.PC1[m], sc.PC2[m], s=45, label=gs, alpha=.9, color=colors[i % len(colors)],
               marker=markers[i % len(markers)], edgecolor="k", linewidth=.3)
for g in stat_ids:
    ax.annotate(short[g], (sc.loc[g, "PC1"], sc.loc[g, "PC2"]), fontsize=6, xytext=(3, 2),
                textcoords="offset points", weight="bold" if g == a.foreground else "normal")
ev = pca.explained_variance_ratio_ * 100
ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
ax.legend(fontsize=6, title="Genomospecies", frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1))
ax.set_title("PCA of standardized % proteome (Pfam categories)")
fig.tight_layout(); fig.savefig(os.path.join(res, "pca_pfam_categories.png"), dpi=300)
print(f"\nPCA: PC1 {ev[0]:.1f}%, PC2 {ev[1]:.1f}%. Results in {res}")
