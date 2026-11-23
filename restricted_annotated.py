#!/usr/bin/env python3
"""Which annotated kdesi-restricted proteins survived the tblastn check?

Joins restricted_tblastn_summary.tsv (from check_restricted_tblastn.py) with the protein names in the
query FASTA and lists, for every protein with a functional annotation (not 'hypothetical protein'),
whether it is still restricted to kdesi or had a homologue in another genome's DNA.

Usage:
  python restricted_annotated.py --summary ~/bdello_analysis/restricted_tblastn/restricted_tblastn_summary.tsv \
      --query restricted_kdesi.faa [--gff ~/bdello_analysis/gff/Bdellovibrio_bacteriovorus_kdesi.gff]
If the FASTA headers carry no product names, give the kdesi GFF/GenBank-derived GFF with --gff.
"""
import argparse, os, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--summary", required=True)
ap.add_argument("--query", required=True, help="the FASTA used as tblastn query")
ap.add_argument("--gff", default=None, help="optional GFF to take product names from")
a = ap.parse_args()

s = pd.read_csv(os.path.expanduser(a.summary), sep="\t")
prod = {}
for l in open(os.path.expanduser(a.query)):
    if l.startswith(">"):
        f = l[1:].strip().split(None, 1)
        prod[f[0]] = f[1] if len(f) > 1 else ""
if a.gff:
    for l in open(os.path.expanduser(a.gff)):
        if l.startswith("##FASTA"):
            break
        if l.startswith("#") or "\t" not in l:
            continue
        x = l.rstrip("\n").split("\t")
        if len(x) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in x[8].split(";") if "=" in kv)
        name = at.get("product", "")
        for key in ("ID", "locus_tag", "protein_id", "Name"):
            k = at.get(key, "").replace("cds-", "")
            if k in s.protein.values and not prod.get(k):
                prod[k] = name
s["product"] = s.protein.map(prod).fillna("")
s["product"] = s["product"].str.replace(r"\s*\[.*?\]$", "", regex=True)
annot = s[~s["product"].str.contains("hypothetical", case=False) & (s["product"].str.strip() != "")]
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 60)
print(f"{len(s)} proteins checked; restricted after tblastn: {(s.status == 'restricted_confirmed').sum()}")
print(f"\nAnnotated proteins ({len(annot)}):")
print(annot[["protein", "product", "status", "n_genomes_with_homologue", "best_genome",
             "best_identity_pct", "best_query_cov_pct"]].to_string(index=False))
lost = s[s.status == "homologue_found"]
print(f"\nAll {len(lost)} proteins with a homologue at the DNA level:")
print(lost[["protein", "product", "n_genomes_with_homologue", "best_genome", "best_identity_pct",
            "best_query_cov_pct"]].to_string(index=False))
