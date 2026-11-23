#!/usr/bin/env python3
"""Step 06c — PAML codeml (v4.10.x) branch-site test (foreground = kdesi) + M0 purifying test.

Per selected orthogroup:
  M0        model=0 NSsites=0, free omega             -> overall omega
  M0_w1     model=0 NSsites=0, omega fixed at 1       -> LRT for purifying selection (omega < 1)
  A_alt     model=2 NSsites=2, fix_omega=0 (starts omega=0.5 and 1.5; best lnL kept)
  A_null    model=2 NSsites=2, fix_omega=1 omega=1
Branch-site LRT: 2*dlnL vs chi2(df=1) (conservative; 50:50 mixture p also reported).
All p-values Benjamini–Hochberg adjusted (alpha = 0.05). BEB sites (P > 0.95) listed.

Two-stage branch-site test (v14, for laptop-scale runs): every orthogroup is first screened with
Model A alternative (starts 0.5 and 1.5) and null with branch lengths fixed at that orthogroup's M0
estimates (fix_blength = 2, about 20x faster). Orthogroups with screening BH q < --confirm-q (default
0.2) are re-run with the full models (branch lengths re-estimated, fix_blength = 1 from the M0 values,
both starts), and the full result replaces the screen. Only full-model results can be called
significant. In simulations (25 taxa) screening LRTs were within ~1 unit of, and mostly below, the
full-model LRTs, so genes far from the screening threshold cannot reach q < 0.05 with full models.
M0_w1 also starts from the M0 branch lengths (re-estimated).

Two phases: (1) M0 / M0_w1 on every selected orthogroup; orthogroups with synonymous
saturation (any branch dS > --ds-max under M0) are excluded and listed; (2) branch-site
models on the remaining orthogroups only. Finally, Fisher's exact test (one-sided, BH) asks
whether positively selected genes are enriched in any functional class / Pfam category.
"""
import argparse, os, re, subprocess, sys
from multiprocessing import Pool
import pandas as pd
from scipy.stats import chi2, fisher_exact
sys.path.insert(0, os.path.dirname(__file__))
from common import read_fasta, bh_fdr

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--strains", required=True, help="space-separated: foreground first")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--alpha", type=float, default=0.05)
ap.add_argument("--codeml", default="codeml")
ap.add_argument("--confirm-q", type=float, default=0.2,
                help="screening BH q below which the full branch-site models are run")
ap.add_argument("--ds-max", type=float, default=3.0,
                help="exclude orthogroups with any M0 branch dS above this (synonymous saturation)")
a = ap.parse_args()
strains = a.strains.split()
wd = os.path.join(a.out, "selection"); run_dir = os.path.join(wd, "codeml_runs")
os.makedirs(run_dir, exist_ok=True)

cand = pd.read_csv(os.path.join(wd, "orthogroup_candidates.tsv"), sep="\t")
sel = cand[cand.selected.astype(str).str.upper() == "TRUE"]
print(f"{len(sel)} orthogroups selected for codeml")

CTL = """seqfile = {seq}
treefile = {tree}
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = {model}
NSsites = {ns}
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = {fixw}
omega = {w}
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
fix_blength = {fb}
method = 0
"""
MODELS = {  # name: (tree, model, NSsites, fix_omega, omega)
    "M0": ("tree_M0.nwk", 0, 0, 0, 0.5),
    "M0_w1": ("tree_M0.nwk", 0, 0, 1, 1),
    "A_alt_w0.5": ("tree_branchsite.nwk", 2, 2, 0, 0.5),
    "A_alt_w1.5": ("tree_branchsite.nwk", 2, 2, 0, 1.5),
    "A_null": ("tree_branchsite.nwk", 2, 2, 1, 1),
    "S_alt_w0.5": ("tree_branchsite.nwk", 2, 2, 0, 0.5),   # screen: branch lengths fixed at M0
    "S_alt_w1.5": ("tree_branchsite.nwk", 2, 2, 0, 1.5),
    "S_null": ("tree_branchsite.nwk", 2, 2, 1, 1),
}
LEN_TREE = re.compile(r"^\(+S\d+: [\d.]+")

def write_len_trees(og):
    """Trees carrying this orthogroup's M0 branch lengths (fixed in the screen, starting values otherwise)."""
    f = os.path.join(run_dir, og, "M0", "out.txt")
    txt = open(f).read() if os.path.exists(f) else ""
    lines = [l.strip() for l in txt.split("\n") if LEN_TREE.match(l.strip())]
    if "lnL(" not in txt or not lines:
        return False
    open(os.path.join(run_dir, og, "tree_m0_len.nwk"), "w").write(lines[0] + "\n")
    open(os.path.join(run_dir, og, "tree_bs_len.nwk"), "w").write(re.sub(r"\bS1:", "S1 #1:", lines[0]) + "\n")
    return True

def done(og, name):
    f = os.path.join(run_dir, og, name, "out.txt")
    return os.path.exists(f) and "lnL(" in open(f).read()

def run(job):
    og, name = job
    d = os.path.join(run_dir, og, name); os.makedirs(d, exist_ok=True)
    if done(og, name):
        return
    tree, model, ns, fixw, w = MODELS[name]
    tree, fb = os.path.join(wd, tree), 0
    if name != "M0":
        lt = os.path.join(run_dir, og, "tree_m0_len.nwk" if name == "M0_w1" else "tree_bs_len.nwk")
        if os.path.exists(lt):
            tree, fb = lt, (2 if name.startswith("S_") else 1)
        elif name.startswith("S_"):
            return                                   # no M0 lengths: screen impossible (M0 failed)
    seqs, _ = read_fasta(os.path.join(wd, "alignments", f"{og}.codon.fasta"))
    with open(os.path.join(d, "aln.phy"), "w") as fh:
        fh.write(f" {len(strains)} {len(seqs[strains[0]])}\n")
        for i, s in enumerate(strains):
            fh.write(f"S{i + 1}  {seqs[s]}\n")       # same aliases as 06b_tree.py
    with open(os.path.join(d, "codeml.ctl"), "w") as fh:
        fh.write(CTL.format(seq="aln.phy", tree=tree, model=model, ns=ns, fixw=fixw, w=w, fb=fb))
    subprocess.run([a.codeml, "codeml.ctl"], cwd=d, input="\n", text=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def parse(path):
    txt = open(path).read() if os.path.exists(path) else ""
    r = {}
    m = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s+(-?[\d.]+)", txt)
    if m:
        r["np"], r["lnL"] = int(m.group(1)), float(m.group(2))
    m = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", txt)
    if m:
        r["omega"] = float(m.group(1))
    m = re.search(r"proportion\s+([\d.\s]+)\n\s*background w\s+([\d.\s]+)\n\s*foreground w\s+([\d.\s]+)", txt)
    if m:
        p, bw, fw = (list(map(float, x.split())) for x in m.groups())
        r.update(p0=p[0], p1=p[1], p2a=p[2], p2b=p[3], w0=bw[0], w2_foreground=fw[2])
    beb = txt.split("Bayes Empirical Bayes (BEB)")
    sites = []
    if len(beb) > 1:
        sec = beb[1].split("Positive sites for foreground lineages")[-1]
        for line in sec.split("\n")[1:]:
            mm = re.match(r"^\s*(\d+)\s+(\S)\s+([\d.]+)(\**)", line)
            if mm:
                sites.append((int(mm.group(1)), mm.group(2), float(mm.group(3))))
            elif sites and not line.strip():
                break
    r["BEB_sites_P>0.95"] = ";".join(f"{n}{aa}({pp:.3f})" for n, aa, pp in sites if pp > 0.95)
    r["n_BEB_P>0.95"] = sum(pp > 0.95 for _, _, pp in sites)
    return r

def m0_branches(path):
    """Per-branch dS and the dS tree length from the M0 output."""
    txt = open(path).read() if os.path.exists(path) else ""
    # branch table rows: "26..1  t  N  S  dN/dS  dN  dS  N*dN  S*dS" (plain decimals, one line each)
    br = re.findall(r"^[ \t]*(\d+)\.\.(\d+)[ \t]+(?:\d+\.\d+[ \t]+){5}(\d+\.\d+)", txt, re.M)
    tl = re.search(r"tree length for dS:\s*([\d.]+)", txt)
    return [(int(a_), int(b_), float(ds)) for a_, b_, ds in br], (float(tl.group(1)) if tl else float("nan"))

alias_f = os.path.join(wd, "codeml_aliases.tsv")
alias = dict(pd.read_csv(alias_f, sep="\t").values) if os.path.exists(alias_f) else {}
ntax = len(strains)
def branch_name(a_, b_):
    return alias.get(f"S{b_}", f"S{b_}") if b_ <= ntax else f"internal {a_}..{b_}"

# ---------- phase 1: M0 and saturation filter ----------
with Pool(a.threads) as p:
    p.map(run, [(og, "M0") for og in sel.orthogroup], chunksize=1)
    for og in sel.orthogroup:
        write_len_trees(og)
    p.map(run, [(og, "M0_w1") for og in sel.orthogroup], chunksize=1)
sat_rows, keep = [], []
for og in sel.orthogroup:
    br, tl = m0_branches(os.path.join(run_dir, og, "M0", "out.txt"))
    if not br:
        sat_rows.append(dict(orthogroup=og, max_branch_dS=None, worst_branch="", tree_length_dS=tl,
                             status="M0_failed")); continue
    a_, b_, mx = max(br, key=lambda x: x[2])
    ok = mx <= a.ds_max
    sat_rows.append(dict(orthogroup=og, max_branch_dS=mx, worst_branch=branch_name(a_, b_),
                         tree_length_dS=tl, status="kept" if ok else "saturated"))
    if ok:
        keep.append(og)
sat = pd.DataFrame(sat_rows)
sat.to_csv(os.path.join(wd, "saturation_filter.tsv"), sep="\t", index=False)
print(f"Saturation filter (max branch dS <= {a.ds_max}): kept {len(keep)} of {len(sel)}; "
      f"saturated {int((sat.status == 'saturated').sum())}; M0 failed {int((sat.status == 'M0_failed').sum())}")
if (sat.status == "saturated").any():
    print("  branches most often responsible: " + ", ".join(
        f"{k} ({v})" for k, v in sat[sat.status == "saturated"].worst_branch.value_counts().head(5).items()))
sel = sel[sel.orthogroup.isin(keep)]

# ---------- phase 2: branch-site models ----------
def lrt(og, alts, null):
    R = [parse(os.path.join(run_dir, og, m, "out.txt")) for m in alts]
    n0 = parse(os.path.join(run_dir, og, null, "out.txt"))
    best = max((r.get("lnL", -1e18) for r in R))
    return 2 * (best - n0["lnL"]) if "lnL" in n0 and best > -1e18 else None

jobs = [(og, m) for og in sel.orthogroup for m in ("S_alt_w0.5", "S_alt_w1.5", "S_null")]
print(f"Screening {len(sel)} orthogroups with branch lengths fixed at M0 ({len(jobs)} codeml runs)")
with Pool(a.threads) as p:
    p.map(run, jobs, chunksize=1)
scr = pd.DataFrame(dict(orthogroup=list(sel.orthogroup),
                        LRT=[lrt(og, ("S_alt_w0.5", "S_alt_w1.5"), "S_null") for og in sel.orthogroup]))
scr["p"] = [chi2.sf(max(0.0, v), 1) if v is not None else float("nan") for v in scr.LRT]
scr["q_BH"] = bh_fdr(scr.p)
scr.to_csv(os.path.join(wd, "branchsite_screen.tsv"), sep="\t", index=False)
confirm = list(scr.orthogroup[(scr.q_BH < a.confirm_q) | scr.p.isna()])
jobs = [(og, m) for og in confirm for m in ("A_alt_w0.5", "A_alt_w1.5", "A_null")]
print(f"Screen: nominal P < 0.05 in {int((scr.p < 0.05).sum())}, BH q < {a.confirm_q} in "
      f"{int((scr.q_BH < a.confirm_q).sum())}, failed {int(scr.p.isna().sum())} of {len(scr)}")
print(f"Full branch-site models for {len(confirm)} orthogroups ({len(jobs)} codeml runs)")
with Pool(a.threads) as p:
    p.map(run, jobs, chunksize=1)

rows = []
for _, c in sel.iterrows():
    og = c.orthogroup
    R = {m: parse(os.path.join(run_dir, og, m, "out.txt")) for m in MODELS}
    full = "lnL" in R["A_null"] and any("lnL" in R[m] for m in ("A_alt_w0.5", "A_alt_w1.5"))
    pre = "A_" if full else "S_"
    alt = max((R[pre + "alt_w0.5"], R[pre + "alt_w1.5"]), key=lambda x: x.get("lnL", -1e18))
    R["A_null"] = R[pre + "null"]
    if not all("lnL" in x for x in (R["M0"], R["M0_w1"], alt, R["A_null"])):
        print(f"WARNING: codeml failed for {og} (see {run_dir}/{og})"); continue
    lrt_bs = max(0.0, 2 * (alt["lnL"] - R["A_null"]["lnL"]))
    lrt_pur = max(0.0, 2 * (R["M0"]["lnL"] - R["M0_w1"]["lnL"]))
    srow = sat[sat.orthogroup == og].iloc[0]
    rows.append(dict(orthogroup=og, selection_class=c.selection_class, pfam_category=c.pfam_category,
                     max_branch_dS=srow.max_branch_dS, tree_length_dS=srow.tree_length_dS,
                     foreground_protein=c.foreground_protein, codons=c.codons,
                     omega_M0=R["M0"]["omega"], lnL_M0=R["M0"]["lnL"], lnL_M0_w1=R["M0_w1"]["lnL"],
                     LRT_purifying=lrt_pur, p_purifying=chi2.sf(lrt_pur, 1),
                     lnL_A_alt=alt["lnL"], lnL_A_null=R["A_null"]["lnL"],
                     branchsite_test="full" if full else "screen_fixed_M0_lengths",
                     LRT_screen=lrt(og, ("S_alt_w0.5", "S_alt_w1.5"), "S_null"), LRT_branchsite=lrt_bs, p_branchsite=chi2.sf(lrt_bs, 1),
                     p_branchsite_mixture=0.5 * chi2.sf(lrt_bs, 1),
                     foreground_w2=alt.get("w2_foreground"), p2a_plus_p2b=alt.get("p2a", 0) + alt.get("p2b", 0),
                     **{k: alt[k] for k in ("n_BEB_P>0.95", "BEB_sites_P>0.95")}))
df = pd.DataFrame(rows)
df["q_branchsite_BH"] = bh_fdr(df.p_branchsite)
df["q_purifying_BH"] = bh_fdr(df.p_purifying)
df["positive_selection_foreground"] = (df.q_branchsite_BH < a.alpha) & (df.foreground_w2 > 1) \
    & (df.branchsite_test == "full")
df["purifying_selection"] = (df.omega_M0 < 1) & (df.q_purifying_BH < a.alpha)
df.to_csv(os.path.join(wd, "codeml_branchsite_results.tsv"), sep="\t", index=False)

n = len(df)
print(f"\nOrthogroups analysed: {n} ({(df.branchsite_test == 'full').sum()} with full branch-site models)")
print(f"Purifying selection (M0 omega < 1, BH q < {a.alpha}): {df.purifying_selection.sum()}/{n}")
print(f"Branch-site positive selection on {strains[0]} (BH q < {a.alpha}, w2 > 1): "
      f"{df.positive_selection_foreground.sum()}/{n}")
print(df.groupby("selection_class")[["purifying_selection", "positive_selection_foreground"]].sum().to_string())

# ---------- functional enrichment of positively selected genes ----------
pos = df.positive_selection_foreground
enr = []
if pos.sum() > 0:
    sets = {}
    for col, label in (("selection_class", "class"), ("pfam_category", "pfam_category")):
        for i_, v in df[col].astype(object).fillna("").astype(str).items():
            for term in [t for t in v.split(",") if t and t != "unclassified" and not t.startswith("multiple:")]:
                sets.setdefault((label, term), set()).add(i_)
    for (label, term), members in sorted(sets.items()):
        inn = df.index.isin(list(members))
        a11, a12 = int((inn & pos).sum()), int((inn & ~pos).sum())
        a21, a22 = int((~inn & pos).sum()), int((~inn & ~pos).sum())
        odds, p = fisher_exact([[a11, a12], [a21, a22]], alternative="greater")
        enr.append(dict(level=label, term=term, n_genes=a11 + a12, n_positive=a11,
                        expected=round((a11 + a12) * pos.mean(), 2), odds_ratio=odds, p=p))
    enr = pd.DataFrame(enr)
if len(enr) == 0:
    print("\nNo functional terms to test for enrichment.")
else:
    enr["q_BH"] = bh_fdr(enr.p)
    enr.sort_values("p").to_csv(os.path.join(wd, "positive_selection_enrichment.tsv"), sep="\t", index=False)
    top = enr.sort_values("p").head(8)
    print("\nEnrichment of positively selected genes (Fisher one-sided, BH):")
    for r in top.itertuples():
        print(f"  {r.level:13s} {r.term:28s} {r.n_positive}/{r.n_genes} positive "
              f"(expected {r.expected})  p={r.p:.4f}  q={r.q_BH:.4f}")
