#!/usr/bin/env python3
"""skani evidence for genomospecies membership of the curated and newly added Bdellovibrio genomes.

Genomes: the curated ingroup (metadata) + NCBI assemblies retained by the screen (new strains of curated
genomospecies or new genomospecies; GTDB genus Bdellovibrio or >=95% ANI to a curated genome; CheckM2
completeness >= 90% and contamination < 5%).
Outputs (--outdir): skani_ani_matrix.tsv, assignment_evidence.tsv (one row per added genome: the rule applied
and the numbers behind it), group_membership.tsv (members, within-group minimum ANI, closest other group),
skani_groups_heatmap.png/.pdf.
"""
import argparse, os, re
from collections import defaultdict
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--skani", required=True)
ap.add_argument("--screen", required=True, help="ncbi_screen_genomes_gtdb.tsv")
ap.add_argument("--checkm2", required=True, help="CheckM2 quality_report.tsv for the candidates")
ap.add_argument("--metadata", required=True)
ap.add_argument("--outdir", default="skani_groups")
ap.add_argument("--threshold", type=float, default=95.0)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
base = lambda p: re.sub(r"\.(fna|fa|fasta|fas|fsa)(\.gz)?$", "", os.path.basename(p))

meta = pd.read_csv(a.metadata, sep="\t", dtype=str).fillna("").set_index("genome")
ing = meta[meta.get("outgroup", pd.Series("", index=meta.index)).str.lower() != "yes"]
scr = pd.read_csv(a.screen, sep="\t")
cm = pd.read_csv(a.checkm2, sep="\t").set_index("Name")
scr["cm_compl"] = scr.accession.map(cm["Completeness"]); scr["cm_cont"] = scr.accession.map(cm["Contamination"])
keep = scr[scr.representative.astype(str).str.lower().eq("true") & scr.cls.isin(["existing_GS", "new_GS"])
           & ((scr.gtdb_genus == "Bdellovibrio") | (scr.cls == "existing_GS"))
           & (scr.cm_compl >= 90) & (scr.cm_cont < 5)].copy()

def short(r):
    s = str(r.strain) if isinstance(r.strain, str) and r.strain.strip() not in ("", "missing", "nan") else ""
    if not s:                                    # e.g. 'Bdellovibrio sp. ArHS' with no strain field
        m = re.match(r"^(?:uncultured )?Bdellovibrio(?: bacteriovorus| sp\.)?\s*(\S+)?$", str(r.organism))
        s = m.group(1) if m and m.group(1) else ""
    if not s:
        return r.accession.split(".")[0]
    s = re.sub(r"(_CONCOCT_v[\d.]+_MAG|_metabat|_concoct|_Refined|_sub)$", "", s)
    return s if len(s) <= 22 else s[:20] + "…"
keep["label"] = keep.apply(short, axis=1)
names = {g: meta.loc[g, "strain"] or g for g in ing.index}
names.update(dict(zip(keep.accession, keep.label)))
kind = {g: "curated" for g in ing.index}
kind.update({r.accession: ("MAG" if r.mag else "isolate") for r in keep.itertuples()})
G = list(ing.index) + list(keep.accession)
print(f"{len(ing)} curated ingroup genomes + {len(keep)} added NCBI genomes "
      f"({(~keep.mag).sum()} isolates, {keep.mag.sum()} MAGs)")

# ---- ANI matrix from the skani edge list (pairs below the screen threshold are absent -> 'below 70')
E = pd.read_csv(a.skani, sep="\t"); cols = {c.lower(): c for c in E.columns}
def gid(p):
    d = os.path.basename(os.path.dirname(p))
    return d if d.startswith(("GCA_", "GCF_")) else base(p)
idx = {g: i for i, g in enumerate(G)}
M = np.full((len(G), len(G)), np.nan); AF = {}
for r in E.itertuples(index=False):
    x, y = gid(getattr(r, cols["ref_file"])), gid(getattr(r, cols["query_file"]))
    if x in idx and y in idx:
        v = float(getattr(r, cols["ani"]))
        i, j = idx[x], idx[y]
        if np.isnan(M[i, j]) or v > M[i, j]:
            M[i, j] = M[j, i] = v
            AF[(x, y)] = AF[(y, x)] = (float(getattr(r, cols["align_fraction_ref"])),
                                       float(getattr(r, cols["align_fraction_query"])))
np.fill_diagonal(M, 100.0)
pd.DataFrame(M, index=[names[g] for g in G], columns=[names[g] for g in G]).round(2) \
  .to_csv(os.path.join(a.outdir, "skani_ani_matrix.tsv"), sep="\t", na_rep="<70")

# ---- final group labels: curated genomospecies; new genomospecies renumbered nGS01.. (single linkage at threshold)
group = {g: meta.loc[g, "genomospecies"] for g in ing.index}
for r in keep[keep.cls == "existing_GS"].itertuples():
    group[r.accession] = r.genomospecies
newg = list(keep[keep.cls == "new_GS"].accession)
parent = {g: g for g in newg}
def find(g):
    while parent[g] != g:
        parent[g] = parent[parent[g]]; g = parent[g]
    return g
for x in newg:
    for y in newg:
        if x < y and M[idx[x], idx[y]] >= a.threshold:
            parent[find(x)] = find(y)
comps = defaultdict(list)
for g in newg:
    comps[find(g)].append(g)
for k, members in enumerate(sorted(comps.values(), key=lambda c: names[c[0]]), 1):
    for g in members:
        group[g] = f"nGS{k:02d}"

# ---- assignment evidence for each added genome
curated_ing = list(ing.index)
rows = []
for r in keep.itertuples():
    i = idx[r.accession]
    vals = [(M[i, idx[c]], c) for c in curated_ing if not np.isnan(M[i, idx[c]])]
    vals.sort(reverse=True)
    best_v, best_c = vals[0] if vals else (np.nan, "")
    by_gs = defaultdict(float)
    for v, c in vals:
        by_gs[group[c]] = max(by_gs[group[c]], v)
    ranked = sorted(by_gs.items(), key=lambda t: -t[1])
    second = ranked[1] if len(ranked) > 1 else ("", np.nan)
    n95 = sum(v >= a.threshold for v, _ in vals)
    af = AF.get((r.accession, best_c), (np.nan, np.nan))
    if r.cls == "existing_GS":
        rule = (f"ANI {best_v:.2f}% (>= {a.threshold:g}%) to {names[best_c]} ({group[best_c]}); "
                f"next genomospecies {second[0]} at {second[1]:.2f}%")
    else:
        rule = (f"below {a.threshold:g}% to all {len(curated_ing)} curated genomes (highest "
                + (f"{best_v:.2f}% to {names[best_c]}, {group[best_c]}" if vals else "<70%") + ")")
    allv = [(M[i, idx[o]], o) for o in G if o != r.accession and not np.isnan(M[i, idx[o]])]
    ba_v, ba_g = max(allv) if allv else (np.nan, "")
    if r.cls == "existing_GS":
        rel = f"same genomospecies (>= {a.threshold:g}%)"
    elif not vals:
        rel = "new species, distant (<70% to all curated genomes; below skani range)"
    elif best_v >= 85:
        rel = f"new species, close relative (85-{a.threshold:g}%)"
    else:
        rel = "new species (70-85%)"
    rows.append(dict(genome=r.label, accession=r.accession, type="MAG" if r.mag else "isolate",
                     assigned_group=group[r.accession], relationship=rel, closest_curated=names.get(best_c, ""),
                     closest_curated_group=group.get(best_c, ""),
                     ani_to_closest_curated=f"{best_v:.2f}" if vals else "<70",
                     closest_any_genome=names.get(ba_g, "") if allv else "",
                     ani_to_closest_any=f"{ba_v:.2f}" if allv else "<70",
                     af_ref=af[0] if vals else "", af_query=af[1] if vals else "",
                     curated_genomes_at_or_above_threshold=n95,
                     second_group=second[0], ani_to_second_group=f"{second[1]:.2f}" if second[0] else "",
                     gtdb_species=str(r.gtdb_taxonomy).split("s__")[-1] if "s__" in str(r.gtdb_taxonomy) else "",
                     checkm2_completeness=r.cm_compl, checkm2_contamination=r.cm_cont,
                     isolation_source=r.isolation_source, rule=rule))
ev = pd.DataFrame(rows).sort_values(["assigned_group", "genome"])
ev.to_csv(os.path.join(a.outdir, "assignment_evidence.tsv"), sep="\t", index=False)

# ---- group membership summary (+ concordance with an independent single-linkage clustering of all genomes)
par = {g: g for g in G}
def f2(g):
    while par[g] != g:
        par[g] = par[par[g]]; g = par[g]
    return g
for x in G:
    for y in G:
        if x < y and M[idx[x], idx[y]] >= a.threshold:
            par[f2(x)] = f2(y)
skc = {g: f2(g) for g in G}
mem = []
for gname in sorted(set(group.values()), key=lambda s: (s.startswith("n"), s)):
    gm = [g for g in G if group[g] == gname]
    ii = [idx[g] for g in gm]
    within = np.nanmin(M[np.ix_(ii, ii)]) if len(gm) > 1 else np.nan
    others = [g for g in G if group[g] != gname]
    sub = M[np.ix_(ii, [idx[o] for o in others])] if others else np.full((1, 1), np.nan)
    out, outg = np.nan, ""
    if not np.all(np.isnan(sub)):
        out = np.nanmax(sub)
        outg = group[others[np.nanargmax(sub) % len(others)]]
    joined = sorted({group[o] for o in G if skc[o] in {skc[g] for g in gm} and group[o] != gname})
    mem.append(dict(group=gname, n=len(gm), curated=sum(kind[g] == "curated" for g in gm),
                    added_isolates=sum(kind[g] == "isolate" for g in gm), added_MAGs=sum(kind[g] == "MAG" for g in gm),
                    members=", ".join(names[g] for g in gm),
                    min_within_ani=round(within, 2) if not np.isnan(within) else None,
                    max_ani_to_other_group=f"{out:.2f}" if not np.isnan(out) else "<70", closest_other_group=outg,
                    merged_by_skani_single_linkage_with=", ".join(joined)))
ms = pd.DataFrame(mem)
ms.to_csv(os.path.join(a.outdir, "group_membership.tsv"), sep="\t", index=False)
print(f"{ms.shape[0]} groups: {sum(~ms.group.str.startswith('n'))} curated genomospecies, "
      f"{sum(ms.group.str.startswith('n'))} new; groups merged by skani single linkage: "
      f"{', '.join(sorted(set(ms[ms.merged_by_skani_single_linkage_with != ''].group))) or 'none'}")

# ---- heatmap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
D = 100 - np.nan_to_num(M, nan=70.0); D = (D + D.T) / 2; np.fill_diagonal(D, 0)
order = leaves_list(linkage(squareform(D, checks=False), "average"))
Mo = M[np.ix_(order, order)]; Go = [G[i] for i in order]
t = (a.threshold - 75) / 25                     # below the threshold: light-to-mid blue; at/above: darkest blue
cmap = LinearSegmentedColormap.from_list("blue", [(0, "#cde2fb"), (t * 0.5, "#86b6ef"), (t - 1e-6, "#2a78d6"),
                                                  (t, "#104281"), (1, "#0d366b")])
cmap.set_bad("#f4f4f2")
n = len(G); s = max(7.0, 0.16 * n + 2.5)
fig, ax = plt.subplots(figsize=(s + 1.6, s))
im = ax.imshow(np.where(np.isnan(Mo), np.nan, np.clip(Mo, 75, 100)), cmap=cmap, vmin=75, vmax=100,
               interpolation="nearest")
labels = []
for g in Go:
    tag = {"curated": "", "isolate": " *", "MAG": " †"}[kind[g]]
    labels.append(f"{names[g]}{tag}")
ax.set_xticks(range(n)); ax.set_yticks(range(n))
ax.set_xticklabels(labels, rotation=90, fontsize=6); ax.set_yticklabels(labels, fontsize=6)
for tl, g in zip(ax.get_yticklabels(), Go):
    if names[g] == "kdesi":
        tl.set_fontweight("bold")
for tl, g in zip(ax.get_xticklabels(), Go):
    if names[g] == "kdesi":
        tl.set_fontweight("bold")
k = 0
while k < n:                                           # box contiguous runs of the same group with >= 2 members
    j = k
    while j + 1 < n and group[Go[j + 1]] == group[Go[k]]:
        j += 1
    if j > k:
        ax.add_patch(Rectangle((k - 0.5, k - 0.5), j - k + 1, j - k + 1, fill=False, ec="#1a1a1a", lw=1.0))
    ax.text(n - 0.2 + 0.6, (k + j) / 2, group[Go[k]], va="center", ha="left", fontsize=6, color="#333333",
            transform=ax.transData, clip_on=False)
    k = j + 1
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.08)
cb.set_label("skani ANI (%)", fontsize=8); cb.ax.tick_params(labelsize=7)
cb.ax.axhline(a.threshold, color="#ffffff", lw=1.0)
ax.set_title("Genomospecies of curated and added Bdellovibrio genomes (skani ANI)", fontsize=9, loc="left")
fig.text(0.01, 0.01, f"* added isolate genome;  † added MAG;  boxes: genomospecies with ≥2 genomes;  darkest blue: "
         f"ANI ≥ {a.threshold:g}% (species threshold);  grey: ANI < 70% (not reported by skani)",
         fontsize=6, color="#555555")
fig.tight_layout(rect=(0, 0.02, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(a.outdir, f"skani_groups_heatmap.{ext}"), dpi=300)
print(f"\nAssignment of added genomes:")
print(ev[["genome", "type", "assigned_group", "relationship", "closest_curated", "ani_to_closest_curated",
          "closest_any_genome", "ani_to_closest_any"]].to_string(index=False))
print(f"\nOutputs in {a.outdir}/")
