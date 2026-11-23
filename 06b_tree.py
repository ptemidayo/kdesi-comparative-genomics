#!/usr/bin/env python3
"""Step 06b — six-taxon tree for codeml.

Either prunes your existing IQ-TREE tree (--species-tree) to the six strains, or builds one
with IQ-TREE from the concatenated protein alignments of all passing single-copy orthogroups.
Writes an unrooted, topology-only tree (codeml re-estimates branch lengths) in two versions:
  tree_M0.nwk         no labels
  tree_branchsite.nwk foreground strain labelled '#1'
"""
import argparse, os, shutil, subprocess, sys
from io import StringIO
import pandas as pd
from Bio import Phylo
sys.path.insert(0, os.path.dirname(__file__))
from common import read_fasta

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--strains", required=True, help="space-separated: foreground first")
ap.add_argument("--species-tree", default="")
ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()
strains = a.strains.split(); fg = strains[0]
wd = os.path.join(a.out, "selection")

tree_file = a.species_tree
if not tree_file:
    cand = pd.read_csv(os.path.join(wd, "orthogroup_candidates.tsv"), sep="\t")
    ogs = cand.loc[cand.status == "pass", "orthogroup"]
    cat = {s: [] for s in strains}
    for og in ogs:
        seqs, _ = read_fasta(os.path.join(wd, "alignments", f"{og}.prot.aln.faa"))
        for s in strains:
            cat[s].append(seqs[s])
    cfa = os.path.join(wd, "concat_sco.faa")
    with open(cfa, "w") as fh:
        for s in strains:
            fh.write(f">{s}\n{''.join(cat[s])}\n")
    print(f"Concatenated {len(ogs)} orthogroups; running IQ-TREE (ModelFinder, 1000 UFBoot)")
    iq = next((b for b in ("iqtree3", "iqtree2", "iqtree") if shutil.which(b)), None)
    if iq is None:
        sys.exit("IQ-TREE not found (looked for iqtree3, iqtree2, iqtree)")
    subprocess.run([iq, "-s", cfa, "-m", "MFP", "-B", "1000", "-T", str(a.threads),
                    "--prefix", os.path.join(wd, "concat_sco"), "-redo", "-quiet"], check=True)
    tree_file = os.path.join(wd, "concat_sco.treefile")

t = Phylo.read(tree_file, "newick")
leaves = [c.name or "" for c in t.get_terminals()]

# ---- map tree tips to genome IDs (exact ID, else by strain name from metadata) ----
import re
meta_f = os.path.join(a.out, "metadata.tsv")
strain_of = {}
if os.path.exists(meta_f):
    md = pd.read_csv(meta_f, sep="\t", dtype=str).fillna("")
    if "strain" in md.columns:
        strain_of = dict(zip(md.genome, md.strain))
from common import tip_matcher
match = tip_matcher(strain_of)

tipmap = {tip: match(tip) for tip in leaves}
pd.DataFrame({"tree_tip": list(tipmap), "genome": [v or "" for v in tipmap.values()]})\
  .to_csv(os.path.join(wd, "tree_tip_mapping.tsv"), sep="\t", index=False)
amb = {k: v for k, v in tipmap.items() if v and v.startswith("AMBIGUOUS")}
mapped = {}
for tip, g in tipmap.items():
    if g and not g.startswith("AMBIGUOUS"):
        mapped.setdefault(g, []).append(tip)
dup = {g: tips for g, tips in mapped.items() if len(tips) > 1 and g in strains}
missing = [s_ for s_ in strains if s_ not in mapped]
if amb or dup or missing:
    if amb:
        print("Ambiguous tips:", amb)
    if dup:
        print("Several tips map to the same genome:", dup)
    if missing:
        print("No tree tip found for:", ", ".join(missing))
    sys.exit(f"Could not match the tree to the codeml strains — see {wd}/tree_tip_mapping.tsv")
renamed = 0
for leaf in t.get_terminals():
    g = tipmap.get(leaf.name or "")
    if g and not g.startswith("AMBIGUOUS") and g != leaf.name:
        leaf.name = g; renamed += 1
print(f"Tree: {len(leaves)} tips; all {len(strains)} codeml strains found"
      + (f" ({renamed} tip names matched by strain name)" if renamed else ""))
for leaf in list(t.get_terminals()):
    if leaf.name not in strains:
        t.prune(leaf)

# unroot: a bifurcating root becomes a trifurcation
root = t.root
if len(root.clades) == 2:
    inner = next((c for c in root.clades if not c.is_terminal()), None)
    if inner is not None:
        root.clades.remove(inner); root.clades.extend(inner.clades)

alias = {s: f"S{i + 1}" for i, s in enumerate(strains)}   # short names for codeml
with open(os.path.join(wd, "codeml_aliases.tsv"), "w") as fh:
    fh.write("alias\tgenome\n" + "".join(f"{v}\t{k}\n" for k, v in alias.items()))

def nwk(c, label):
    if c.is_terminal():
        return alias[c.name] + (" #1" if label and c.name == fg else "")
    return "(" + ",".join(nwk(x, label) for x in c.clades) + ")"

for fname, lab in (("tree_M0.nwk", False), ("tree_branchsite.nwk", True)):
    with open(os.path.join(wd, fname), "w") as fh:
        fh.write(nwk(t.root, lab) + ";\n")
print(open(os.path.join(wd, "tree_branchsite.nwk")).read().strip())
print("Aliases: " + ", ".join(f"{v}={k}" for k, v in alias.items()))
