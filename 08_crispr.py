#!/usr/bin/env python3
"""Step 08 — CRISPR-Cas survey across all genomes.

1. cas gene detection: curated regular expressions on Prokka gene/product names
   + Pfam Cas-family domains (Cas_*, CRISPR_*, Csm*, Cmr*, ...) from step 03.
2. CRISPR direct-repeat arrays: MinCED (CRT algorithm; >= 3 repeats by default).
Genome classified: 'complete' (array + >= 2 cas genes incl. cas1/cas2 or effector),
'orphan_array', 'cas_only', or 'none'.
"""
import argparse, os, re, subprocess, sys
from collections import Counter
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import genome_ids

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()
wd = os.path.join(a.out, "crispr"); os.makedirs(wd, exist_ok=True)
ids = genome_ids(a.out)

CAS_RE = re.compile(
    r"\b(cas\d{1,2}[a-z]?\d?|cas3'{0,2}|csn\d|csm\d|cmr\d|cse\d|csy\d|csa\d|csb\d|csc\d|csf\d|"
    r"csx\d+|csh\d|cpf1|c2c\d|crispr[- ]associated)\b|CRISPR[- ]associated|CRISPR system", re.I)
PFAM_CAS_RE = re.compile(r"^(Cas_|CRISPR_|Csm\d|Cmr\d|Cse\d|Csy\d|Csa\d|Csn\d|Cas\d|RAMPs|DevR)", re.I)
CORE_RE = re.compile(r"\bcas(1|2|9|12|13)[a-z]?\b|Cas_Cas1|CRISPR_Cas2|Cas9|Cas12|Cas13|cpf1", re.I)

dom_f = os.path.join(a.out, "pfam", "results", "protein_domains_filtered.tsv")
dom = pd.read_csv(dom_f, sep="\t") if os.path.exists(dom_f) else pd.DataFrame(columns=["genome", "protein", "pfam_name"])

cas_rows, arr_rows, summ = [], [], []
for g in ids:
    # --- cas genes from annotation
    hits = {}
    for line in open(os.path.join(a.out, "gff", f"{g}.gff")):
        if line.startswith("##FASTA"):
            break
        p = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(p) < 9 or p[2] != "CDS":
            continue
        at = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        txt = f"{at.get('gene', '')} {at.get('product', '')}"
        if CAS_RE.search(txt):
            hits[at["ID"]] = dict(genome=g, protein=at["ID"], contig=p[0], start=int(p[3]), end=int(p[4]),
                                  annotation=txt.strip(), pfam="")
    for prot, sub in dom[(dom.genome == g) & dom.pfam_name.str.match(PFAM_CAS_RE)].groupby("protein"):
        h = hits.setdefault(prot, dict(genome=g, protein=prot, contig="", start=None, end=None, annotation="", pfam=""))
        h["pfam"] = ",".join(sorted(set(sub.pfam_name)))
    cas_rows += hits.values()

    # --- CRISPR arrays with MinCED
    mout = os.path.join(wd, f"{g}.minced.txt")
    if not os.path.exists(mout):
        subprocess.run(["minced", os.path.join(a.out, "fna", f"{g}.fna"), mout,
                        os.path.join(wd, f"{g}.minced.gff")], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    seq, arr = None, None
    for line in open(mout):
        m = re.match(r"Sequence '(\S+)'", line)
        if m:
            seq = m.group(1)
        m = re.match(r"CRISPR (\d+)\s+Range: (\d+) - (\d+)", line)
        if m:
            arr = dict(genome=g, contig=seq, array=int(m.group(1)), start=int(m.group(2)),
                       end=int(m.group(3)), repeats=[], n_repeats=0)
        if arr is not None and re.match(r"^\d+\s+[ACGTN]+", line):
            arr["repeats"].append(line.split()[1])
        m = re.match(r"Repeats: (\d+)", line)
        if m and arr is not None:
            arr["n_repeats"] = int(m.group(1))
            arr["consensus_repeat"] = Counter(arr.pop("repeats")).most_common(1)[0][0]
            arr["repeat_length"] = len(arr["consensus_repeat"])
            arr_rows.append(arr); arr = None

    ga = [r for r in arr_rows if r["genome"] == g]
    gc = list(hits.values())
    core = any(CORE_RE.search(f"{h['annotation']} {h['pfam']}") for h in gc)
    status = ("complete" if ga and len(gc) >= 2 and core else "orphan_array" if ga
              else "cas_only" if gc else "none")
    summ.append(dict(genome=g, n_arrays=len(ga), total_repeats=sum(r["n_repeats"] for r in ga),
                     n_cas_genes=len(gc), core_cas_present=core, status=status))

pd.DataFrame(cas_rows).to_csv(os.path.join(wd, "cas_genes.tsv"), sep="\t", index=False)
pd.DataFrame(arr_rows).to_csv(os.path.join(wd, "crispr_arrays.tsv"), sep="\t", index=False)
S = pd.DataFrame(summ); S.to_csv(os.path.join(wd, "crispr_cas_summary.tsv"), sep="\t", index=False)
print(S.to_string(index=False))
print("\n" + S.status.value_counts().to_string())
