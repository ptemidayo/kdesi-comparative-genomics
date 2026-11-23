#!/usr/bin/env python3
"""Supplementary Figures S1-S3 as PNG (600 dpi), SVG and PDF.

  Figure S1  REALPHY core-SNP phylogeny of the 41 genomes, rooted on the outgroups
             (needs the REALPHY tree file; found automatically or given with --tree)
  Figure S2  skani ANI heatmap: converts the existing skani_groups_heatmap.pdf to SVG and PNG
             (pdftocairo/pdftoppm from poppler), so the figure is exactly the one already made
  Figure S3  PCA of the functional profiles of the 37 ingroup genomes (pfam/results/pca_scores.tsv)

Usage:  python make_suppl_figs.py [--tree /path/to/realphy_tree] [--outdir ~/bdello_analysis/supplementary_figures]
"""
import argparse, glob, os, re, shutil, subprocess
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.path.expanduser
ap = argparse.ArgumentParser()
ap.add_argument("--tree", default=None, help="REALPHY (RAxML/PhyML) tree in Newick format")
ap.add_argument("--metadata", default=H("~/bdello_pipeline/metadata_genomospecies.tsv"))
ap.add_argument("--pca-dir", default=H("~/bdello_analysis/pfam/results"))
ap.add_argument("--skani-pdf", default=None)
ap.add_argument("--outgroups", default="JSS,qaytius,HCB154")
ap.add_argument("--outdir", default=H("~/bdello_analysis/supplementary_figures"))
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "svg.fonttype": "none",
                     "axes.edgecolor": "#6b6b66", "axes.labelcolor": "#2b2b29",
                     "xtick.color": "#4a4a46", "ytick.color": "#4a4a46"})
BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#9c9b95", "#2b2b29"
meta = pd.read_csv(H(a.metadata), sep="\t")
done = []

def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(a.outdir, f"{name}.{ext}"), dpi=600, bbox_inches="tight")
    plt.close(fig); done.append(name)

def strain_of(label):
    lab = label.replace(".fasta", "").replace(".fna", "").replace(".fa", "").strip("'\"")
    for r in meta.itertuples():
        if lab == r.genome or lab.endswith("_" + str(r.strain)) or lab == str(r.strain):
            return r.strain, r.genomospecies
    return lab, ""

# ======================= Figure S1: REALPHY tree =======================
def find_tree():
    if a.tree:
        return H(a.tree)
    pats = ["**/*realphy*/**/*RAxML_bestTree*", "**/*RealPhy*/**/*RAxML_bestTree*", "**/*realphy*/**/*.tree",
            "**/*RealPhy*/**/*.tree", "**/*realphy*/**/*.nwk", "**/*realphy*/**/*phyml_tree*"]
    for p in pats:
        hits = glob.glob(os.path.join(H("~"), p), recursive=True)
        if hits:
            return sorted(hits, key=os.path.getmtime)[-1]
    return None

def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = 0; nid = [0]; nodes = {}
    def new(name="", length=0.0):
        nid[0] += 1; nodes[nid[0]] = dict(name=name, len=length, kids=[]); return nid[0]
    def label():
        nonlocal pos
        m = re.match(r"([^:,()]*)(?::([0-9.eE+-]+))?", s[pos:]); pos += m.end()
        return m.group(1).strip("'\""), float(m.group(2)) if m.group(2) else 0.0
    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1; kids = [node()]
            while s[pos] == ",":
                pos += 1; kids.append(node())
            pos += 1  # ')'
            name, ln = label(); n = new(name, ln); nodes[n]["kids"] = kids; return n
        name, ln = label(); return new(name, ln)
    root = node(); return nodes, root

def reroot(nodes, root, outs):
    # undirected edges with lengths; internal labels are support values attached to the child side of an edge
    adj = {n: [] for n in nodes}
    for p, d in nodes.items():
        for c in d["kids"]:
            sup = nodes[c]["name"] if nodes[c]["kids"] else ""
            adj[p].append((c, nodes[c]["len"], sup)); adj[c].append((p, nodes[c]["len"], sup))
    leaves = {n for n, d in nodes.items() if not d["kids"]}
    tipname = {n: nodes[n]["name"] for n in leaves}
    def side(u, v):  # leaves reachable from v without passing u
        st, seen, out = [v], {u, v}, set()
        while st:
            x = st.pop()
            if x in leaves: out.add(x)
            for y, _, _ in adj[x]:
                if y not in seen: seen.add(y); st.append(y)
        return out
    target = {n for n in leaves if any(o in tipname[n] for o in outs)}
    best = None
    for u in adj:
        for v, ln, sup in adj[u]:
            if side(u, v) == target:
                best = (u, v, ln); break
        if best: break
    if not best:
        print("  outgroups are not monophyletic in this tree; drawing it as rooted in the file")
        return nodes, root
    u, v, ln = best
    new = {}
    k = [0]
    def build(x, parent, length, sup):
        k[0] += 1; me = k[0]
        kids = [build(y, x, l, s_) for y, l, s_ in adj[x] if y != parent]
        new[me] = dict(name=tipname.get(x, sup), len=length, kids=kids); return me
    r = 10**6
    new[r] = dict(name="", len=0.0, kids=[build(v, u, ln / 2, ""), build(u, v, ln / 2, "")])
    return new, r

tree_file = find_tree()
if tree_file:
    print("Figure S1 tree:", tree_file)
    nodes, root = parse_newick(open(tree_file).read())
    nodes, root = reroot(nodes, root, a.outgroups.split(","))
    xs, ys, order = {}, {}, []
    def setx(n, x):
        xs[n] = x
        for c in nodes[n]["kids"]: setx(c, x + nodes[c]["len"])
    setx(root, 0.0)
    def sety(n):
        if not nodes[n]["kids"]:
            ys[n] = len(order); order.append(n); return ys[n]
        ys[n] = sum(sety(c) for c in nodes[n]["kids"]) / len(nodes[n]["kids"]); return ys[n]
    sety(root)
    ntip = len(order); maxx = max(xs.values())
    fig, ax = plt.subplots(figsize=(7.2, 0.19 * ntip + 0.8))
    for n, d in nodes.items():
        for c in d["kids"]:
            ax.plot([xs[n], xs[c]], [ys[c], ys[c]], color=INK, lw=0.8)
        if d["kids"]:
            ax.plot([xs[n], xs[n]], [min(ys[c] for c in d["kids"]), max(ys[c] for c in d["kids"])], color=INK, lw=0.8)
            sup = d["name"]
            if sup and re.match(r"^[0-9.]+$", sup) and n != root:
                v = float(sup); v = v * 100 if v <= 1 else v
                if v >= 70:
                    ax.text(xs[n] - maxx * 0.005, ys[n] - 0.25, f"{v:.0f}", ha="right", va="bottom", fontsize=5.5, color="#6b6b66")
    for n in order:
        st, gs = strain_of(nodes[n]["name"])
        is_k = "kdesi" in nodes[n]["name"]
        ax.text(xs[n] + maxx * 0.01, ys[n], f"{st}" + (f"  ({gs})" if gs else ""), va="center", fontsize=7,
                color=BLUE if is_k else INK, fontweight="bold" if is_k else "normal")
    sb = 10 ** (len(str(int(maxx * 100))) - 3) if maxx >= 0.01 else 0.001
    sb = max(sb, round(maxx / 10, 3)) or 0.01
    ax.plot([0, sb], [ntip + 0.6, ntip + 0.6], color=INK, lw=1)
    ax.text(sb / 2, ntip + 1.1, f"{sb:g} substitutions per site", ha="center", va="top", fontsize=6.5)
    ax.set_xlim(-maxx * 0.02, maxx * 1.35); ax.set_ylim(ntip + 1.6, -1); ax.axis("off")
    save(fig, "FigureS1_REALPHY_tree")
else:
    print("Figure S1: REALPHY tree not found - give it with --tree /path/to/tree")

# ======================= Figure S2: skani heatmap (convert) =======================
pdf = H(a.skani_pdf) if a.skani_pdf else next(iter(sorted(glob.glob(H("~/**/skani_groups_heatmap.pdf"), recursive=True))), None)
if pdf and shutil.which("pdftocairo"):
    base = os.path.join(a.outdir, "FigureS2_skani_heatmap")
    shutil.copy(pdf, base + ".pdf")
    subprocess.run(["pdftocairo", "-svg", pdf, base + ".svg"], check=True)
    subprocess.run(["pdftoppm", "-png", "-r", "600", "-singlefile", pdf, base], check=True)
    done.append("FigureS2_skani_heatmap"); print("Figure S2 from:", pdf)
else:
    print("Figure S2: skani_groups_heatmap.pdf or pdftocairo (conda install -c conda-forge poppler) not found")

# ======================= Figure S3: PCA =======================
sc_f, var_f = os.path.join(H(a.pca_dir), "pca_scores.tsv"), os.path.join(H(a.pca_dir), "pca_variance.tsv")
if os.path.exists(sc_f):
    sc = pd.read_csv(sc_f, sep="\t", index_col=0)
    var = pd.read_csv(var_f, sep="\t", index_col=0).iloc[:, 0] if os.path.exists(var_f) else None
    size_f = os.path.join(H(a.pca_dir), "proteome_sizes.tsv")
    size = pd.read_csv(size_f, sep="\t", index_col=0).iloc[:, 0] if os.path.exists(size_f) else None
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    rel = {"GS09", "GS10"}
    for g in sc.index:
        st, gs = strain_of(g)
        is_k = "kdesi" in g
        col = BLUE if is_k else (ORANGE if gs in rel else GREY)
        ms = 8 if is_k else 6
        ax.plot(sc.loc[g, "PC1"], sc.loc[g, "PC2"], "o", ms=ms, color=col, mec="white", mew=0.8, zorder=3 if is_k else 2)
        if is_k or gs in rel:
            ax.annotate(st, (sc.loc[g, "PC1"], sc.loc[g, "PC2"]), xytext=(5, 3), textcoords="offset points",
                        fontsize=6.5, color=INK, fontweight="bold" if is_k else "normal")
    lab = lambda i: f"PC{i} ({100 * var.iloc[i - 1]:.1f}% of variance)" if var is not None else f"PC{i}"
    ax.set_xlabel(lab(1)); ax.set_ylabel(lab(2))
    ax.axhline(0, color="#e4e3dd", lw=0.6, zorder=0); ax.axvline(0, color="#e4e3dd", lw=0.6, zorder=0)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.plot([], [], "o", color=BLUE, label="kdesi (GS05)"); ax.plot([], [], "o", color=ORANGE, label="closest relatives (GS09, GS10)")
    ax.plot([], [], "o", color=GREY, label="other ingroup genomes")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3,
              handletextpad=0.2, columnspacing=0.8)
    if size is not None:
        rho = sc.PC1.corr(size.reindex(sc.index), method="spearman")
        print(f"Figure S3: PC1 vs proteome size, Spearman rho = {rho:.2f} (manuscript: -0.73)")
    save(fig, "FigureS3_functional_PCA")
else:
    print("Figure S3: pca_scores.tsv not found in", a.pca_dir)

print("\nSaved in", a.outdir, ":", ", ".join(done) if done else "nothing")
