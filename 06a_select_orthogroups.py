#!/usr/bin/env python3
"""Step 06a — pick single-copy orthogroups for the branch-site test.

For every OrthoFinder single-copy orthogroup across the six strains:
  * functional class from the foreground protein's Pfam category (pfam_categories.tsv: selection_class)
  * MAFFT protein alignment -> nucleotide-guided back-translation to a codon alignment
  * filters: >= 50 codons, < 30% gaps, >= 10 segregating sites
Then samples N orthogroups per class (fixed seed) and writes PHYLIP codon alignments for codeml.
All passing candidates are listed so you can hand-pick instead.
"""
import argparse, glob, os, subprocess, sys, random
from multiprocessing import Pool
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import read_fasta

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--strains", required=True, help="space-separated: foreground first")
ap.add_argument("--orthofinder-dir", required=True)
ap.add_argument("--per-class", type=int, default=6)
ap.add_argument("--mode", choices=["all", "per_class"], default="all",
                help="all = every single-copy orthogroup passing the filters; per_class = N per class")
ap.add_argument("--min-codons", type=int, default=50)
ap.add_argument("--max-gap", type=float, default=0.30)
ap.add_argument("--min-seg", type=int, default=10)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--seed", type=int, default=42)
a = ap.parse_args()

strains = a.strains.split(); fg = strains[0]
wd = os.path.join(a.out, "selection"); aln_dir = os.path.join(wd, "alignments")
os.makedirs(aln_dir, exist_ok=True)

runs = sorted(glob.glob(os.path.join(a.orthofinder_dir, "OrthoFinder", "Results_*", "Orthogroups", "Orthogroups.tsv")))
if not runs:
    sys.exit(f"No OrthoFinder Orthogroups.tsv under {a.orthofinder_dir}/OrthoFinder — check logs/orthofinder.log")
res = os.path.dirname(os.path.dirname(runs[-1]))
og = pd.read_csv(os.path.join(res, "Orthogroups", "Orthogroups.tsv"), sep="\t", dtype=str).set_index("Orthogroup")
missing_cols = [x for x in strains if x not in og.columns]
if missing_cols:
    sys.exit(f"Orthogroups.tsv lacks columns for: {missing_cols}")
og = og[strains]
sco_f = os.path.join(res, "Orthogroups", "Orthogroups_SingleCopyOrthologues.txt")
if os.path.exists(sco_f):
    sco = [l.strip() for l in open(sco_f) if l.strip()]
else:   # OrthoFinder 3 no longer writes this file: derive single-copy orthogroups ourselves
    one = og.apply(lambda col: col.fillna("").astype(str).map(lambda v: v.strip() != "" and "," not in v))
    sco = list(og.index[one.all(axis=1)])
    print(f"(single-copy orthogroups derived from Orthogroups.tsv: exactly one gene in each of {len(strains)} genomes)")
og = og.loc[sco]
og = og.apply(lambda col: col.astype(str).str.strip().str.split().str[0])   # gene ID only (drop product text)
print(f"{len(og)} single-copy orthogroups across {len(strains)} strains ({res})")

pc = pd.read_csv(os.path.join(a.out, "pfam", "results", "protein_categories.tsv"), sep="\t").fillna("")
pc = pc[pc.genome == fg]
cat = pc.groupby("protein").category.apply(lambda s: ",".join(sorted(set(s))))   # all Pfam categories
pcs = pc[pc.selection_class != ""]
cls = pcs.groupby("protein").selection_class.apply(lambda s: sorted(set(s)))

prot = {s: read_fasta(os.path.join(a.out, "faa", f"{s}.faa"))[0] for s in strains}
nuc = {s: read_fasta(os.path.join(a.out, "ffn", f"{s}.ffn"))[0] for s in strains}
STOPS = {"TAA", "TAG", "TGA"}

def process(item):
    ogid, members = item
    fgp = members[fg]
    c = cls.get(fgp)
    if c is None:
        if a.mode == "per_class":
            return None
        c = ["unclassified"]
    elif len(c) != 1:
        if a.mode == "per_class":
            return None
        c = ["multiple:" + "+".join(c)]
    pfa = os.path.join(aln_dir, f"{ogid}.prot.faa")
    with open(pfa, "w") as fh:
        for s in strains:
            fh.write(f">{s}\n{prot[s][members[s]].rstrip('*')}\n")
    aln = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", pfa],
                         capture_output=True, text=True, check=True).stdout
    paln = {}
    for blk in aln.split(">")[1:]:
        h, *sq = blk.split("\n"); paln[h.strip()] = "".join(sq).upper()
    caln = {}
    for s in strains:
        nt = nuc[s][members[s]].upper()
        if len(nt) >= 3 and nt[-3:] in STOPS:
            nt = nt[:-3]
        if len(nt) != 3 * len(paln[s].replace("-", "")):
            return dict(orthogroup=ogid, status=f"length_mismatch_{s}")
        i, out = 0, []
        for aa in paln[s]:
            if aa == "-":
                out.append("---")
            else:
                cod = nt[i:i + 3]; i += 3
                out.append("---" if cod in STOPS else cod)   # internal stops masked
        caln[s] = "".join(out)
    L = len(caln[fg]); ncod = L // 3
    gaps = sum(v.count("-") for v in caln.values()) / (L * len(strains))
    seg = sum(1 for j in range(L)
              if all(caln[s][j] != "-" for s in strains) and len({caln[s][j] for s in strains}) > 1)
    ok = ncod >= a.min_codons and gaps < a.max_gap and seg >= a.min_seg
    with open(os.path.join(aln_dir, f"{ogid}.codon.fasta"), "w") as fh:
        for s in strains:
            fh.write(f">{s}\n{caln[s]}\n")
    with open(os.path.join(aln_dir, f"{ogid}.prot.aln.faa"), "w") as fh:
        for s in strains:
            fh.write(f">{s}\n{paln[s]}\n")
    return dict(orthogroup=ogid, selection_class=c[0], pfam_category=cat.get(fgp, ""),
                foreground_protein=fgp, codons=ncod, gap_fraction=round(gaps, 4),
                segregating_sites=seg, status="pass" if ok else "fail_filters",
                **{f"protein_{s}": members[s] for s in strains})

with Pool(a.threads) as p:
    rows = [r for r in p.map(process, [(i, og.loc[i].to_dict()) for i in og.index]) if r]
cand = pd.DataFrame(rows)
passing = cand[cand.status == "pass"]
print(passing.groupby("selection_class").size().rename("passing_orthogroups").to_string())

random.seed(a.seed)
chosen = []
if a.mode == "all":
    chosen = sorted(passing.orthogroup)
for c, sub in ([] if a.mode == "all" else passing.groupby("selection_class")):
    ids = sorted(sub.orthogroup)
    pick = ids if len(ids) <= a.per_class else sorted(random.sample(ids, a.per_class))
    if len(ids) < a.per_class:
        print(f"WARNING: class {c} has only {len(ids)} passing orthogroups")
    chosen += pick
cand["selected"] = cand.orthogroup.isin(chosen)
cand.to_csv(os.path.join(wd, "orthogroup_candidates.tsv"), sep="\t", index=False)

print(f"Selected {len(chosen)} orthogroups (column 'selected' in orthogroup_candidates.tsv).")
print("To hand-pick (e.g. exactly 40), edit TRUE/FALSE in that column and rerun step 06c only.")
