#!/usr/bin/env python3
"""Figure 5: the kdesi-specific features at a glance.

(A) Chemoreceptor (MCP) number in every ingroup genome, grouped by genomospecies (kdesi highlighted,
    its closest relatives Bd5 and GS10 marked) - shows MCP number is fixed within species and reduced in kdesi.
(B) Map of the kdesi chromosome: genomic islands, tRNA/tmRNA genes, MCP genes, the two RM pairs unique to
    kdesi, the MetNIQ methionine transporter operon and the divergent CheR-CheD pair.

Inputs (defaults are the pipeline paths):
  --pfam      pfam/results/protein_domains_filtered.tsv   (MCP = proteins with PF00015)
  --metadata  metadata_genomospecies.tsv                    (genome, strain, genomospecies, outgroup)
  --gff       Prokka GFF of kdesi                           (tRNA, MCP, Che and Met gene positions)
  --islands   kdesi_trna.tsv from islands_trna.py           (island coordinates)
Usage:  python fig5_kdesi.py --islands /path/to/kdesi_trna.tsv
"""
import argparse, glob, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.expanduser
ap = argparse.ArgumentParser()
ap.add_argument("--pfam", default=H("~/bdello_analysis/pfam/results/protein_domains_filtered.tsv"))
ap.add_argument("--metadata", default=H("~/bdello_pipeline/metadata_genomospecies.tsv"))
ap.add_argument("--gff", default=H("~/bdello_analysis/gff/Bdellovibrio_bacteriovorus_kdesi.gff"))
ap.add_argument("--islands", default=None, help="kdesi_trna.tsv (searched under ~ if not given)")
ap.add_argument("--che-pair", default="kdesi_02253,kdesi_02254", help="locus tags of the divergent CheR-CheD pair")
ap.add_argument("--met", default="kdesi_01482,kdesi_01483,kdesi_01484", help="metN, metP, metQ locus tags")
ap.add_argument("--rm", default="1039648-1041791,2657766-2660355", help="RM pair coordinates (start-end)")
ap.add_argument("--outdir", default=".")
a = ap.parse_args()
isl_path = a.islands or next(iter(glob.glob(H("~/**/kdesi_trna.tsv"), recursive=True)), None)
for lab, p in (("--pfam", a.pfam), ("--metadata", a.metadata), ("--gff", a.gff), ("--islands", isl_path)):
    if not p or not os.path.exists(H(p)):
        raise SystemExit(f"{lab} file not found: {p}")
print("islands:", isl_path)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": "#6b6b66",
                     "axes.labelcolor": "#2b2b29", "xtick.color": "#4a4a46", "ytick.color": "#4a4a46"})
BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#9c9b95", "#2b2b29"

# ---------------- Panel A data: MCP number per genome
dom = pd.read_csv(H(a.pfam), sep="\t")
dom["acc"] = dom.pfam_acc.astype(str).str.split(".").str[0]
mcp = dom[dom.acc == "PF00015"].groupby("genome").protein.nunique()
meta = pd.read_csv(H(a.metadata), sep="\t")
meta = meta[meta.outgroup.astype(str).str.lower() != "yes"].copy()
if "exclude_from_stats" in meta:          # strain W: excluded from all trait statistics, so also from the figure
    drop = meta.exclude_from_stats.astype(str).str.lower() == "yes"
    print("excluded as in the statistics:", ", ".join(meta.loc[drop, "strain"].astype(str)))
    meta = meta[~drop].copy()
meta["mcp"] = meta.genome.map(mcp).fillna(0).astype(int)
kd = meta[meta.genome.str.contains("kdesi")].iloc[0]
order = (meta.groupby("genomospecies").mcp.median().sort_values(ascending=False).index.tolist())
xpos = {g: i for i, g in enumerate(order)}

# ---------------- Panel B data: kdesi chromosome features
feat, clen = [], {}
for line in open(H(a.gff)):
    if line.startswith("##FASTA"):
        break
    if line.startswith("##sequence-region"):
        x = line.split(); clen[x[1]] = int(x[3]); continue
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] not in ("CDS", "tRNA", "tmRNA"):
        continue
    at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    feat.append(dict(type=f[2], start=int(f[3]), end=int(f[4]), locus=at.get("locus_tag", at.get("ID", ""))))
feat = pd.DataFrame(feat)
L = max(clen.values()) if clen else int(feat.end.max())
kd_mcp = dom[(dom.acc == "PF00015") & dom.genome.str.contains("kdesi")].protein.unique()
mcp_pos = feat[feat.locus.isin(kd_mcp)]
rna = feat[feat.type.isin(["tRNA", "tmRNA"])]
isl = pd.read_csv(H(isl_path), sep="\t")
def span(tags):
    s = feat[feat.locus.isin(tags.split(","))]
    return (int(s.start.min()), int(s.end.max())) if len(s) else None
che, met = span(a.che_pair), span(a.met)
rm = [tuple(map(int, r.split("-"))) for r in a.rm.split(",")]
print(f"MCP genes located: {len(mcp_pos)} of {len(kd_mcp)}; tRNA/tmRNA {len(rna)}; islands {len(isl)}; "
      f"CheR-CheD {che}; MetNIQ {met}")

# ---------------- figure
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.6)
ax = fig.add_subplot(gs[0])
rng = np.random.default_rng(3)
rel = {"GS09", "GS10"}
for r in meta.itertuples():
    x = xpos[r.genomospecies] + rng.uniform(-0.18, 0.18)
    is_k = "kdesi" in r.genome
    col = BLUE if is_k else (ORANGE if r.genomospecies in rel else GREY)
    ax.plot(x, r.mcp, "o", ms=7 if is_k else 5, color=col, mec="white", mew=1, zorder=3 if is_k else 2)
for g, sub in meta.groupby("genomospecies"):
    if len(sub) > 1:
        ax.plot([xpos[g] - 0.3, xpos[g] + 0.3], [sub.mcp.median()] * 2, color=INK, lw=1, zorder=1)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=6.5)
ax.set_ylabel("Chemoreceptor (MCP) genes per genome"); ax.set_xlim(-0.7, len(order) - 0.3)
ax.grid(axis="y", color="#e4e3dd", lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.annotate(f"kdesi ({kd.mcp})", (xpos[kd.genomospecies], kd.mcp), xytext=(10, -14), textcoords="offset points",
            fontsize=7.5, color=INK, fontweight="bold", arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
ax.plot([], [], "o", color=BLUE, label="kdesi (GS05)"); ax.plot([], [], "o", color=ORANGE, label="closest relatives (GS09, GS10)")
ax.plot([], [], "o", color=GREY, label="other genomes"); ax.plot([], [], color=INK, lw=1, label="median (genomospecies with >1 genome)")
ax.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=4, handletextpad=0.3, columnspacing=1.0)
ax.yaxis.get_major_locator().set_params(integer=True)
ax.set_xlabel("Genomospecies (ordered by median MCP number)")
ax.text(-0.07, 1.12, "A", transform=ax.transAxes, fontweight="bold", fontsize=10)

bx = fig.add_subplot(gs[1])
bx.set_xlim(0, L / 1e6); bx.set_ylim(-1.6, 2.2); bx.set_yticks([])
for s in ("top", "right", "left"):
    bx.spines[s].set_visible(False)
bx.set_xlabel("Position on the kdesi chromosome (Mb)")
bx.add_patch(plt.Rectangle((0, -0.08), L / 1e6, 0.16, color="#d9d8d2", lw=0))
for r in isl.itertuples():
    bx.add_patch(plt.Rectangle((r.start / 1e6, -0.32), (r.end - r.start) / 1e6, 0.64, color=ORANGE, lw=0, zorder=3))
for r in rna.itertuples():
    bx.plot([r.start / 1e6] * 2, [0.38, 0.62], color="#6b6b66", lw=0.7)
for r in mcp_pos.itertuples():
    bx.plot([r.start / 1e6] * 2, [-0.38, -0.62], color=BLUE, lw=0.9)
labels = [((s + e) / 2e6, "RM pair (NspV-like)" if i == 0 else "RM pair (AccI-like)") for i, (s, e) in enumerate(rm)]
if met: labels.append(((met[0] + met[1]) / 2e6, "MetNIQ transporter"))
if che: labels.append(((che[0] + che[1]) / 2e6, "divergent CheR–CheD"))
labels.sort()
ylev = [1.05, 1.55, 2.05]
for i, (x, txt) in enumerate(labels):
    y = ylev[i % len(ylev)]
    bx.plot([x, x], [0.15, y - 0.12], color=INK, lw=0.6)
    bx.plot(x, 0, marker="v", color=INK, ms=4, zorder=4)
    bx.text(x, y, txt, ha="center", va="bottom", fontsize=7, color=INK)
bx.text(0, -1.05, "islands", color=ORANGE, fontsize=7, fontweight="bold")
bx.text(0.45, -1.05, "| tRNA / tmRNA genes (above)", color="#6b6b66", fontsize=7)
bx.text(1.75, -1.05, "| MCP genes (below)", color=BLUE, fontsize=7)
bx.set_title("B", loc="left", fontweight="bold", fontsize=10)
fig.subplots_adjust(left=0.09, right=0.98, top=0.90, bottom=0.08)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(a.outdir, f"Figure5_kdesi_features.{ext}"), dpi=600)
print(meta.groupby("genomospecies").mcp.agg(["count", "min", "max"]).loc[order].to_string())
print(f"Saved Figure5_kdesi_features.png/.pdf in {os.path.abspath(a.outdir)}")
