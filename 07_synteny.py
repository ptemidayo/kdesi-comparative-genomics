#!/usr/bin/env python3
"""Step 07a — synteny of the foreground (single closed chromosome) vs the type strain HD100.

Bidirectional best hits (BBH) with DIAMOND blastp — the same reciprocal-best-hit logic as the
RAST sequence-based comparison — then gene-order collinearity and a dot plot.
(If you prefer to reproduce the Methods literally, run RAST 'Sequence Based Comparison'
on the web and keep this as a scriptable cross-check.)
"""
import argparse, os, subprocess, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--query", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--min-block", type=int, default=5)
ap.add_argument("--max-step", type=int, default=3, help="max gene-rank gap allowed inside a block")
a = ap.parse_args()
wd = os.path.join(a.out, "synteny"); os.makedirs(wd, exist_ok=True)
faa = lambda g: os.path.join(a.out, "faa", f"{g}.faa")

def blast(q, s):
    o = os.path.join(wd, f"{q}_vs_{s}.tsv")
    if not os.path.exists(o):
        db = os.path.join(wd, s)
        subprocess.run(["diamond", "makedb", "--in", faa(s), "-d", db, "--quiet"], check=True)
        subprocess.run(["diamond", "blastp", "-q", faa(q), "-d", db, "-o", o, "--more-sensitive",
                        "-e", "1e-5", "--max-target-seqs", "5", "-p", str(a.threads), "--quiet",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
                        "qcovhsp"], check=True)
    d = pd.read_csv(o, sep="\t", header=None,
                    names=["q", "s", "pident", "len", "evalue", "bits", "qcov"])
    return d.sort_values("bits", ascending=False).drop_duplicates("q")

f, r = blast(a.query, a.ref), blast(a.ref, a.query)
best_r = dict(zip(r.q, r.s))
bbh = f[f.apply(lambda x: best_r.get(x.s) == x.q, axis=1)].rename(columns={"q": "query_gene", "s": "ref_gene"})

def gene_order(g):
    rows, contigs = [], []
    for line in open(os.path.join(a.out, "gff", f"{g}.gff")):
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        if p[0] not in contigs:
            contigs.append(p[0])
        gid = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv).get("ID")
        rows.append((gid, contigs.index(p[0]), int(p[3]), p[6]))
    d = pd.DataFrame(rows, columns=["gene", "contig_idx", "start", "strand"]).sort_values(["contig_idx", "start"])
    d["rank"] = np.arange(len(d))
    return d.set_index("gene"), len(contigs)

qo, nq = gene_order(a.query); ro, nr = gene_order(a.ref)
bbh = bbh[bbh.query_gene.isin(qo.index) & bbh.ref_gene.isin(ro.index)].copy()
bbh["q_rank"] = bbh.query_gene.map(qo["rank"]); bbh["r_rank"] = bbh.ref_gene.map(ro["rank"])
bbh["q_start"] = bbh.query_gene.map(qo.start); bbh["r_start"] = bbh.ref_gene.map(ro.start)
bbh = bbh.sort_values("q_rank").reset_index(drop=True)

# collinear blocks along query order (either orientation; circular wrap ignored)
blocks, cur = [], [0]
def same_block(i, j, direction):
    step = bbh.r_rank[j] - bbh.r_rank[i]
    return 0 < abs(step) <= a.max_step and (direction == 0 or np.sign(step) == direction) \
        and bbh.q_rank[j] - bbh.q_rank[i] <= a.max_step
direction = 0
for j in range(1, len(bbh)):
    if same_block(cur[-1], j, direction):
        direction = int(np.sign(bbh.r_rank[j] - bbh.r_rank[cur[-1]])); cur.append(j)
    else:
        blocks.append((cur, direction)); cur, direction = [j], 0
blocks.append((cur, direction))
brows = [dict(block=n + 1, n_genes=len(c), orientation="forward" if d >= 0 else "inverted",
              query_start_gene=bbh.query_gene[c[0]], query_end_gene=bbh.query_gene[c[-1]],
              ref_start_gene=bbh.ref_gene[c[0]], ref_end_gene=bbh.ref_gene[c[-1]])
         for n, (c, d) in enumerate([b for b in blocks if len(b[0]) >= a.min_block])]
bdf = pd.DataFrame(brows)
bbh.to_csv(os.path.join(wd, f"BBH_{a.query}_vs_{a.ref}.tsv"), sep="\t", index=False)
bdf.to_csv(os.path.join(wd, f"collinear_blocks_{a.query}_vs_{a.ref}.tsv"), sep="\t", index=False)

in_blocks = bdf.n_genes.sum() if len(bdf) else 0
summary = {
    "query_contigs": nq, "ref_contigs": nr, "query_CDS": len(qo), "ref_CDS": len(ro),
    "BBH_pairs": len(bbh), "BBH_percent_query_CDS": round(100 * len(bbh) / len(qo), 2),
    "median_BBH_identity": round(bbh.pident.median(), 2),
    f"collinear_blocks(>={a.min_block}_genes)": len(bdf),
    "inverted_blocks": int((bdf.orientation == "inverted").sum()) if len(bdf) else 0,
    "percent_BBH_in_blocks": round(100 * in_blocks / max(len(bbh), 1), 2),
    "spearman_rank_corr": round(bbh[["q_rank", "r_rank"]].corr("spearman").iloc[0, 1], 4)}
pd.Series(summary).to_csv(os.path.join(wd, "synteny_summary.tsv"), sep="\t", header=False)
print(pd.Series(summary).to_string())

fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.scatter(bbh.q_start / 1e6, bbh.r_start / 1e6, s=1.5, c=bbh.pident, cmap="viridis")
ax.set_xlabel(f"{a.query} position (Mb)"); ax.set_ylabel(f"{a.ref} position (Mb)")
ax.set_title(f"BBH synteny: {a.query} vs {a.ref} ({len(bbh)} pairs)")
fig.tight_layout(); fig.savefig(os.path.join(wd, f"synteny_dotplot_{a.query}_vs_{a.ref}.png"), dpi=300)
