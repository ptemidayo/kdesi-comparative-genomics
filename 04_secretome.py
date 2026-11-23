#!/usr/bin/env python3
"""Step 04b — Summarise SignalP-6 predictions and benchmark the foreground secretome.

Outputs ($OUT/signalp/results):
  signalp_all_proteins.tsv       one row per protein with SP type
  signalp_summary.tsv            per genome: counts and % proteome per SP type
  secreted_by_category.tsv       per genome x Pfam category: SP-bearing proteins
  secretome_benchmark.tsv        foreground vs benchmark strains (median, range, z, rank)
"""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import genome_ids

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--foreground", default="kdesi")
ap.add_argument("--benchmark", default="", help="space-separated genome IDs; empty = all others")
a = ap.parse_args()
res = os.path.join(a.out, "signalp", "results"); os.makedirs(res, exist_ok=True)
ids = genome_ids(a.out)

rows = []
for g in ids:
    f = os.path.join(a.out, "signalp", g, "prediction_results.txt")
    for line in open(f):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        rows.append((g, p[0].split()[0], p[1]))
sp = pd.DataFrame(rows, columns=["genome", "protein", "prediction"])
sp.to_csv(os.path.join(res, "signalp_all_proteins.tsv"), sep="\t", index=False)

summ = sp.groupby(["genome", "prediction"]).size().unstack(fill_value=0).reindex(ids).fillna(0).astype(int)
summ["total_proteins"] = sp.groupby("genome").size().reindex(ids)
types = [c for c in summ.columns if c not in ("OTHER", "total_proteins")]
summ["signal_peptide_total"] = summ[types].sum(axis=1)
summ["signal_peptide_percent"] = (summ.signal_peptide_total / summ.total_proteins * 100).round(3)
summ.to_csv(os.path.join(res, "signalp_summary.tsv"), sep="\t")

feat = summ[types + ["signal_peptide_total", "signal_peptide_percent"]].copy()
cat_f = os.path.join(a.out, "pfam", "results", "protein_categories.tsv")
if os.path.exists(cat_f):
    pc = pd.read_csv(cat_f, sep="\t")
    sec = pc.merge(sp[sp.prediction != "OTHER"], on=["genome", "protein"])
    bycat = sec.groupby(["genome", "category"]).size().unstack(fill_value=0).reindex(ids).fillna(0).astype(int)
    bycat.to_csv(os.path.join(res, "secreted_by_category.tsv"), sep="\t")
    feat = feat.join(bycat.add_prefix("secreted_"))

fg = a.foreground
bench = a.benchmark.split() if a.benchmark.strip() else [g for g in ids if g != fg]
bench = [b for b in bench if b in feat.index and b != fg]
B = feat.loc[bench]
out = pd.DataFrame({
    fg: feat.loc[fg],
    "benchmark_median": B.median(), "benchmark_min": B.min(), "benchmark_max": B.max(),
    "z_vs_benchmark": (feat.loc[fg] - B.mean()) / B.std(ddof=1).replace(0, np.nan),
    "rank_of_foreground(1=highest)": [int((feat.loc[bench + [fg], c] > feat.loc[fg, c]).sum() + 1) for c in feat.columns],
    "n_benchmark": len(bench)})
out.round(3).to_csv(os.path.join(res, "secretome_benchmark.tsv"), sep="\t")
print(f"Secretome benchmark: {fg} vs {len(bench)} strains ({', '.join(bench)})")
print(out.round(2).to_string())

# ---------- single-case tests against the reference genomes (as in the Pfam analysis) ----------
from common import load_metadata, single_case_test, single_case_regression, bh_fdr
meta, stat = load_metadata(a.out, ids)            # drops exclude_from_stats (W) and outgroups
ref = [g for g in stat if g != fg]
cnt_cols = types + ["signal_peptide_total"]
pct = summ[cnt_cols].div(summ.total_proteins, axis=0) * 100
ch = single_case_test(pct.loc[fg], pct.loc[ref])
ch.insert(0, "measure", "percent_of_proteome")
reg = []
for c in cnt_cols:
    r = single_case_regression(summ.loc[fg, "total_proteins"], summ.loc[fg, c],
                               summ.loc[ref, "total_proteins"], summ.loc[ref, c])
    reg.append(dict(feature=c, **r))
reg = pd.DataFrame(reg)
reg["p_bonferroni"] = (reg.p * reg.p.notna().sum()).clip(upper=1)
tests = ch.merge(reg.add_prefix("sizeadj_").rename(columns={"sizeadj_feature": "feature"}), on="feature")
tests.to_csv(os.path.join(res, "foreground_secretome_tests.tsv"), sep="\t", index=False)
print(f"\n{fg} vs {len(ref)} reference ingroup genomes (excluded: outgroups, exclude_from_stats):")
show = tests[["feature", "foreground", "ref_mean", "ref_min", "ref_max", "z", "percentile", "p_bonferroni",
              "sizeadj_observed", "sizeadj_expected", "sizeadj_t", "sizeadj_p_bonferroni"]]
print(show.round(3).to_string(index=False))
if os.path.exists(cat_f):
    rows = []
    for c in bycat.columns:
        y = bycat[c]
        if (y.loc[ref] > 0).sum() < 5:
            continue
        r = single_case_regression(summ.loc[fg, "total_proteins"], y.loc[fg],
                                   summ.loc[ref, "total_proteins"], y.loc[ref])
        rows.append(dict(category=c, **r))
    cr = pd.DataFrame(rows)
    if len(cr):
        cr["q_BH"] = bh_fdr(cr.p)
        cr.sort_values("p").to_csv(os.path.join(res, "foreground_secreted_by_category_sizeadj.tsv"), sep="\t", index=False)
        print("\nSecreted proteins per Pfam category, size-adjusted (exploratory, BH):")
        print(cr.sort_values("p").head(6)[["category", "observed", "expected", "t", "p", "q_BH"]].round(3).to_string(index=False))
