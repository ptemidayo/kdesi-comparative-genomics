#!/usr/bin/env python3
"""Which functional traits are species-level, and how many Bdellovibrio species remain unsampled?
Part A: share of variance between genomospecies (eta^2) per trait, across ingroup genomes of genomospecies with
>= 2 members (W and outgroups excluded), permutation test of genomospecies labels, Benjamini-Hochberg.
Part B: Chao1 (bias-corrected, log-normal 95% CI), sample coverage, extrapolation (Chao et al. 2014),
accumulation curve from 1,000 random genome orders over all 38 ingroup genomes."""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import bh_fdr

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--metadata", required=True)
ap.add_argument("--perm", type=int, default=9999)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
wd = os.path.join(a.out, "species_traits"); os.makedirs(wd, exist_ok=True)
rng = np.random.default_rng(a.seed)

meta = pd.read_csv(a.metadata, sep="\t", dtype=str).fillna("").set_index("genome")
yes = lambda c: meta[c].str.lower().eq("yes") if c in meta.columns else pd.Series(False, index=meta.index)
out_mask, exc_mask = yes("outgroup"), yes("exclude_from_stats")
ingroup = meta[~out_mask]
stat = meta[~out_mask & ~exc_mask]

R = os.path.join(a.out, "pfam", "results")
rd = lambda f: pd.read_csv(f, sep="\t", index_col=0)
counts = rd(os.path.join(R, "category_protein_counts.tsv"))
pct = rd(os.path.join(R, "category_percent_proteome.tsv"))
size = rd(os.path.join(R, "proteome_sizes.tsv")).iloc[:, 0]
mcp = rd(os.path.join(R, "mcp_domain_profile.tsv"))["MCP_proteins(PF00015)"]
genomes = [g for g in stat.index if g in counts.index]

def resid(y):
    x = size.loc[genomes].astype(float); yy = y.loc[genomes].astype(float)
    b, c0 = np.polyfit(x, yy, 1)
    return yy - (c0 + b * x)

traits = {}
traits[("genome", "proteome_size")] = size.loc[genomes].astype(float)
for c in counts.columns:
    traits[("pfam_percent", c)] = pct.loc[genomes, c]
    traits[("pfam_size_adjusted", c)] = resid(counts[c])
traits[("chemoreceptors", "MCP_proteins")] = mcp.loc[genomes].astype(float)
traits[("chemoreceptors", "MCP_proteins_size_adjusted")] = resid(mcp)
sp_f = os.path.join(a.out, "signalp", "results", "signalp_summary.tsv")
if os.path.exists(sp_f):
    sp = rd(sp_f)
    traits[("secretome", "signal_peptide_percent")] = sp.loc[genomes, "signal_peptide_percent"].astype(float)
    for c in [c for c in sp.columns if c not in ("OTHER", "total_proteins", "signal_peptide_total",
                                                  "signal_peptide_percent")]:
        traits[("secretome", f"{c}_percent")] = (sp.loc[genomes, c] / sp.loc[genomes, "total_proteins"] * 100)
cai_f = os.path.join(a.out, "cai", "median_cai_strain_x_category.tsv")
if os.path.exists(cai_f):
    cai = rd(cai_f)
    for c in cai.columns:
        traits[("codon_usage_CAI", c)] = cai.loc[[g for g in genomes if g in cai.index], c]

gs = stat.loc[genomes, "genomospecies"]
multi = gs[gs.map(gs.value_counts()) >= 2]
n, k = len(multi), multi.nunique()
codes = pd.factorize(multi.values)[0]

def eta2(y, codes):
    y = np.asarray(y, float); tot = ((y - y.mean()) ** 2).sum()
    if tot == 0:
        return np.nan
    w = sum(((y[codes == j] - y[codes == j].mean()) ** 2).sum() for j in np.unique(codes))
    return 1 - w / tot

perms = [rng.permutation(codes) for _ in range(a.perm)]
rows = []
for (fam, name), y in traits.items():
    y = y.reindex(multi.index)
    if y.isna().any():
        continue
    yv = y.values; obs = eta2(yv, codes)
    null = np.array([eta2(yv, p) for p in perms])
    p = (1 + np.sum(null >= obs)) / (1 + len(null))
    within = max(np.ptp(yv[codes == j]) for j in np.unique(codes))
    allv = traits[(fam, name)].astype(float)
    rows.append(dict(family=fam, trait=name, eta2_between_genomospecies=obs, perm_p=p,
                     max_within_genomospecies_range=within, genus_range=float(np.ptp(allv)),
                     within_over_genus_range=within / np.ptp(allv) if np.ptp(allv) else np.nan))
res = pd.DataFrame(rows)
res["q_BH"] = bh_fdr(res.perm_p)
res = res.sort_values(["q_BH", "eta2_between_genomospecies"], ascending=[True, False])
res.to_csv(os.path.join(wd, "species_level_traits.tsv"), sep="\t", index=False)
print(f"Part A: {len(res)} traits; {n} genomes in {k} genomospecies with >=2 members "
      f"(eta2 expected by chance = {(k - 1) / (n - 1):.2f}); {a.perm} permutations")
print(f"Species-level (q < 0.05): {(res.q_BH < 0.05).sum()} of {len(res)}")
for fam, sub in res.groupby("family", sort=False):
    print(f"  {fam:20s} {(sub.q_BH < 0.05).sum():2d} of {len(sub):2d} significant; median eta2 "
          f"{sub.eta2_between_genomospecies.median():.2f}")
show = res[["family", "trait", "eta2_between_genomospecies", "q_BH", "within_over_genus_range"]].copy()
print(show.head(25).round(3).to_string(index=False))
print("...\nLeast species-level:")
print(show.tail(6).round(3).to_string(index=False))

x = ingroup["genomospecies"].value_counts().values
N, S = int(x.sum()), len(x)
f1, f2 = int((x == 1).sum()), int((x == 2).sum())
f0 = (N - 1) / N * f1 * (f1 - 1) / (2 * (f2 + 1))
chao = S + f0
var = (f1 * (f1 - 1) / (2 * (f2 + 1)) + f1 * (2 * f1 - 1) ** 2 / (4 * (f2 + 1) ** 2)
       + f1 ** 2 * f2 * (f1 - 1) ** 2 / (4 * (f2 + 1) ** 4))
K = np.exp(1.96 * np.sqrt(np.log(1 + var / f0 ** 2))) if f0 > 0 else 1
cov = 1 - f1 / N * ((N - 1) * f1 / ((N - 1) * f1 + 2 * f2))
ext = lambda m: S + f0 * (1 - (1 - f1 / (N * f0 + f1)) ** m) if f0 > 0 else S
rich = pd.DataFrame([dict(genomes=N, genomospecies_observed=S, singletons=f1, doubletons=f2,
                          chao1=chao, chao1_lo=S + f0 / K, chao1_hi=S + f0 * K,
                          percent_observed=100 * S / chao, sample_coverage=cov,
                          p_next_genome_new=1 - cov, expected_at_2N=ext(N), expected_at_5N=ext(4 * N))])
rich.round(3).to_csv(os.path.join(wd, "genomospecies_richness.tsv"), sep="\t", index=False)
print(f"\nPart B: {S} genomospecies from {N} genomes ({f1} singletons, {f2} doubletons)")
print(f"  Chao1 = {chao:.1f} (95% CI {S + f0 / K:.1f}-{S + f0 * K:.1f}); observed = {100 * S / chao:.0f}% of estimate")
print(f"  sample coverage = {cov:.3f}; probability that the next genome is a new genomospecies = {1 - cov:.2f}")
print(f"  expected genomospecies with {2 * N} genomes: {ext(N):.1f}; with {5 * N}: {ext(4 * N):.1f}")

labels = np.repeat(np.arange(S), x)
curves = np.array([[len(set(o[:i])) for i in range(1, N + 1)]
                   for o in (rng.permutation(labels) for _ in range(1000))])
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(5.5, 4))
xs = np.arange(1, N + 1)
ax.fill_between(xs, np.percentile(curves, 2.5, 0), np.percentile(curves, 97.5, 0), color="#9ecae1", alpha=0.6,
                lw=0, label="Observed (95% of random orders)")
ax.plot(xs, curves.mean(0), color="#08519c", lw=2)
xe = np.arange(N, 5 * N + 1)
ax.plot(xe, [ext(m - N) for m in xe], color="#08519c", lw=2, ls="--", label="Extrapolated (Chao et al. 2014)")
ax.axhline(chao, color="#a50f15", lw=1, ls=":", label=f"Chao1 = {chao:.0f}")
ax.set_xlabel("Genomes sampled"); ax.set_ylabel("Genomospecies")
ax.legend(frameon=False, fontsize=8, loc="lower right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(wd, "genomospecies_accumulation.png"), dpi=300)
print(f"\nResults and figure in {wd}")
