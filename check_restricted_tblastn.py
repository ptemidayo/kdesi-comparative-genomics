#!/usr/bin/env python3
"""Confirm that kdesi-restricted proteins have no homologue in the other genomes at the DNA level.

BLASTp against predicted proteins can miss homologues that were not annotated (missed ORFs,
frameshifts, pseudogenes). This script searches each restricted protein with tblastn against the
GENOME sequences of all other genomes and applies the same rule as the protein search:
a protein stays "restricted" only if no genome has a hit reaching >= MIN_ID % identity over
>= MIN_COV of the protein length. HSPs from the same genome region are merged on the query, so a
homologue broken by a frameshift into several pieces is still counted over its full length.

Usage:
  python3 scripts/check_restricted_tblastn.py --query restricted_kdesi.faa \
      --genomes /home/sarah_elufisan/bdello_genome --exclude kdesi \
      --out /home/sarah_elufisan/bdello_analysis/restricted_tblastn --threads 3
Outputs (in --out): tblastn_hits.tsv (all HSPs), restricted_tblastn_summary.tsv (one row per protein).
"""
import argparse, glob, os, subprocess, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--query", required=True, help="FASTA of kdesi-restricted proteins")
ap.add_argument("--genomes", required=True, help="folder with genome FASTA files (.fna/.fa/.fasta/.fas)")
ap.add_argument("--exclude", default="kdesi", help="skip genome files whose name contains this text")
ap.add_argument("--out", required=True)
ap.add_argument("--evalue", type=float, default=1e-5)
ap.add_argument("--min-id", type=float, default=30.0)
ap.add_argument("--min-cov", type=float, default=0.50)
ap.add_argument("--window", type=int, default=3000, help="max distance (bp) between HSPs merged as one locus")
ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

# ---- query lengths ----
qlen, cur = {}, None
for l in open(a.query):
    if l.startswith(">"):
        cur = l[1:].split()[0]; qlen[cur] = 0
    elif cur:
        qlen[cur] += len(l.strip().rstrip("*"))
print(f"{len(qlen)} restricted proteins in {a.query}")

# ---- one nucleotide database of all other genomes, contigs tagged by genome ----
files = sorted(f for ext in ("fna", "fa", "fasta", "fas") for f in glob.glob(os.path.join(a.genomes, f"*.{ext}")))
files = [f for f in files if a.exclude.lower() not in os.path.basename(f).lower()]
db_fa = os.path.join(a.out, "other_genomes.fna")
with open(db_fa, "w") as out:
    for f in files:
        g = os.path.basename(f).rsplit(".", 1)[0]
        n = 0
        for l in open(f):
            if l.startswith(">"):
                n += 1; out.write(f">{g}|c{n}\n")
            else:
                out.write(l.strip() + "\n")
print(f"Database: {len(files)} genomes (excluded files containing '{a.exclude}')")
subprocess.run(["makeblastdb", "-in", db_fa, "-dbtype", "nucl", "-out", os.path.join(a.out, "other_genomes")],
               check=True, stdout=subprocess.DEVNULL)

hits_f = os.path.join(a.out, "tblastn_hits.tsv")
cols = "qseqid sseqid pident nident length qstart qend sstart send evalue bitscore"
subprocess.run(["tblastn", "-query", a.query, "-db", os.path.join(a.out, "other_genomes"),
                "-evalue", str(a.evalue), "-seg", "no", "-max_target_seqs", "5000",
                "-num_threads", str(a.threads), "-outfmt", "6 " + cols, "-out", hits_f], check=True)

# ---- merge HSPs per (query, genome, contig, nearby region) ----
loci = defaultdict(list)
for l in open(hits_f):
    q, s, pid, nid, ln, qs, qe, ss, se, ev, bs = l.rstrip("\n").split("\t")
    genome, contig = s.split("|", 1)
    mid = (int(ss) + int(se)) // 2
    key = (q, genome, contig)
    for locus in loci[key]:
        if abs(locus["mid"] - mid) <= a.window:
            break
    else:
        locus = {"mid": mid, "hsps": []}; loci[key].append(locus)
    locus["hsps"].append((int(qs), int(qe), int(nid), int(ln)))

def merged_cov(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if cs is None or s > ce + 1:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs + 1 if cs is not None else 0)

best = {}   # query -> (passes, cov, ident, genome)
pass_genomes = defaultdict(set)
for (q, genome, contig), lst in loci.items():
    for locus in lst:
        h = locus["hsps"]
        cov = merged_cov([(s, e) for s, e, _, _ in h]) / qlen[q]
        ident = 100 * sum(n for *_, n, _ in h) / sum(l for *_, l in h)
        ok = cov >= a.min_cov and ident >= a.min_id
        if ok:
            pass_genomes[q].add(genome)
        score = (ok, cov * ident)
        if q not in best or score > best[q][0]:
            best[q] = (score, cov, ident, genome)

rows = []
for q in qlen:
    if q in best:
        (_, _), cov, ident, g = best[q]
    else:
        cov, ident, g = 0.0, 0.0, ""
    status = "homologue_found" if pass_genomes[q] else "restricted_confirmed"
    rows.append((q, qlen[q], status, len(pass_genomes[q]), g, round(ident, 1), round(100 * cov, 1),
                 ",".join(sorted(pass_genomes[q]))))
summ = os.path.join(a.out, "restricted_tblastn_summary.tsv")
with open(summ, "w") as fh:
    fh.write("protein\tlength_aa\tstatus\tn_genomes_with_homologue\tbest_genome\tbest_identity_pct\t"
             "best_query_cov_pct\tgenomes_with_homologue\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")

conf = sum(r[2] == "restricted_confirmed" for r in rows)
print(f"\nRestricted after tblastn (no hit >= {a.min_id}% identity over >= {100*a.min_cov:.0f}% length): "
      f"{conf} of {len(rows)}")
for r in rows:
    if r[2] == "homologue_found":
        print(f"  homologue found: {r[0]} ({r[1]} aa) in {r[3]} genome(s); best {r[4]} "
              f"{r[5]}% identity over {r[6]}%")
print(f"Summary: {summ}")
