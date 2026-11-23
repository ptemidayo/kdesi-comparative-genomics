#!/usr/bin/env python3
"""Check whether an existing tree can be used for codeml (step 06).
Usage: python3 scripts/check_tree.py <tree file>      (run from the bdello_pipeline folder)
Reports: number of tips, how each tip maps to a genome, missing genomes, and whether the
codeml representatives (one per genomospecies) are all present."""
import os, subprocess, sys
import pandas as pd
from Bio import Phylo
sys.path.insert(0, os.path.dirname(__file__))
from common import tip_matcher

tree = sys.argv[1]
OUT = subprocess.run(["bash", "-c", f"source {os.path.dirname(__file__)}/../config.sh && echo $OUT"],
                     capture_output=True, text=True).stdout.strip()
md = pd.read_csv(os.path.join(OUT, "metadata.tsv"), sep="\t", dtype=str).fillna("")
match = tip_matcher(dict(zip(md.genome, md.strain)))
t = Phylo.read(tree, "newick")
tips = [c.name or "" for c in t.get_terminals()]
rows = [(tip, match(tip) or "") for tip in tips]
print(f"{tree}\n{len(tips)} tips; rooted as given: {t.rooted}")
print(f"{'tree tip':45s} -> genome")
for tip, g in rows:
    print(f"  {tip[:43]:43s} -> {g or 'NO MATCH'}")
found = {g for _, g in rows if g and not g.startswith("AMBIGUOUS")}
miss = md[~md.genome.isin(found)]
print(f"\nGenomes not in the tree: {', '.join(miss.strain) if len(miss) else 'none'}")
yes = lambda c: md[c].str.lower().eq("yes") if c in md.columns else False
pool = md[~(yes("outgroup") | yes("exclude_from_stats"))]
gs_missing = sorted(set(pool.genomospecies) - set(pool[pool.genome.isin(found)].genomospecies))
print(f"Genomospecies with no genome in the tree: {', '.join(gs_missing) if gs_missing else 'none'}")
problems = [r for r in rows if not r[1] or r[1].startswith("AMBIGUOUS")]
print("\nRESULT:", "usable for codeml" if not gs_missing and not
      [r for r in problems if r[1].startswith("AMBIGUOUS")] else "NOT usable as is — send this output to Claude")
