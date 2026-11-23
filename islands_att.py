#!/usr/bin/env python3
"""Look for att-site direct repeats at tRNA-associated genomic islands.

When an integrase inserts an element into a tRNA (or tmRNA) gene, the core attachment sequence (attB) at
the 3' end of the tRNA is duplicated: one copy stays with the tRNA, the other marks the far end of the
inserted element. A direct repeat of the tRNA 3' end near the opposite island boundary is therefore the
signature of site-specific integration.

For every tRNA/tmRNA inside an island or within --flank bp of it (from islands_trna.py output), the
script takes the 3' end of the gene (last --tail nt plus --down nt beyond it, in the gene's orientation),
and finds the longest exact direct repeat of that sequence in the island +/- --pad bp, ignoring
matches that overlap any tRNA/tmRNA gene (tRNA genes share conserved motifs). It reports where the
repeat lies (near the opposite island edge or not) and the chance of a match that long by luck.

Control: the same search around every tRNA/tmRNA NOT associated with an island, using a window of the
median island search size, shows how often repeats of each length appear without an island.

Usage:
  python islands_att.py --islands-tsv kdesi_trna.tsv --gff ~/bdello_analysis/gff/Bdellovibrio_bacteriovorus_kdesi.gff
"""
import argparse, math
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--islands-tsv", required=True, help="output of islands_trna.py (<islands>_trna.tsv)")
ap.add_argument("--gff", required=True, help="Prokka GFF with the ##FASTA section")
ap.add_argument("--flank", type=int, default=1000, help="tRNA within this distance of an island is tested")
ap.add_argument("--tail", type=int, default=50, help="nt from the 3' end of the tRNA used as query")
ap.add_argument("--down", type=int, default=20, help="nt beyond the 3' end added to the query")
ap.add_argument("--pad", type=int, default=2000, help="search the island extended by this many bp each side")
ap.add_argument("--min-len", type=int, default=15, help="shortest repeat reported as a candidate att site")
ap.add_argument("--edge", type=int, default=2000, help="'near the opposite edge' = within this many bp")
ap.add_argument("--out", default=None)
a = ap.parse_args()

# ---- read GFF features and sequences
feat, seq, cur, fasta = [], {}, None, False
for line in open(a.gff):
    line = line.rstrip("\n")
    if line.startswith("##FASTA"):
        fasta = True; continue
    if fasta:
        if line.startswith(">"):
            cur = line[1:].split()[0]; seq[cur] = []
        elif cur:
            seq[cur].append(line.strip().upper())
        continue
    if line.startswith("#") or not line.strip():
        continue
    f = line.split("\t")
    if len(f) >= 9 and f[2] in ("tRNA", "tmRNA"):
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        feat.append(dict(contig=f[0], type=f[2], start=int(f[3]), end=int(f[4]), strand=f[6],
                         product=attr.get("product", "").replace("%2C", ",")))
seq = {k: "".join(v) for k, v in seq.items()}
if not seq:
    raise SystemExit("no ##FASTA section in the GFF; Prokka GFFs normally include one")
rna = pd.DataFrame(feat)
rna["label"] = rna.type + " " + rna["product"].str.replace("tRNA-", "", regex=False)
isl = pd.read_csv(a.islands_tsv, sep="\t")
print(f"{len(seq)} contigs; {len(rna)} tRNA/tmRNA genes; {len(isl)} islands")

def query(r):
    """3' end of the gene (+ downstream), in forward genome coordinates (1-based inclusive)."""
    if r.strand == "+":
        s, e = r.end - a.tail + 1, r.end + a.down
    else:
        s, e = r.start - a.down, r.start + a.tail - 1
    return max(1, s), e

def masked(c):
    m = np.zeros(len(seq[c]) + 2, bool)
    for r in rna[rna.contig == c].itertuples():
        m[r.start:r.end + 1] = True
    return m
MASK = {c: masked(c) for c in seq}

def longest_repeat(c, qs, qe, ws, we, k=10):
    """Longest exact match of genome[qs..qe] inside genome[ws..we], same orientation, outside tRNA genes
    and not overlapping the query itself. Returns (length, start, end, sequence)."""
    g, m = seq[c], MASK[c]
    q = g[qs - 1:qe]; best = (0, None, None, "")
    idx = {}
    for i in range(len(q) - k + 1):
        idx.setdefault(q[i:i + k], []).append(i)
    ws, we = max(1, ws), min(len(g), we)
    for j in range(ws - 1, we - k + 1):
        for i in idx.get(g[j:j + k], ()):
            # extend
            l, r_ = 0, k
            while i - l - 1 >= 0 and j - l - 1 >= ws - 1 and q[i - l - 1] == g[j - l - 1]:
                l += 1
            while i + r_ < len(q) and j + r_ < we and q[i + r_] == g[j + r_]:
                r_ += 1
            s0, e0 = j - l + 1, j + r_          # 1-based genome coordinates of the match
            L = e0 - s0 + 1
            if L <= best[0]:
                continue
            if not (e0 < qs or s0 > qe):        # overlaps the query copy itself
                continue
            if m[s0:e0 + 1].any():              # inside another tRNA/tmRNA gene
                continue
            best = (L, s0, e0, g[s0 - 1:e0])
    return best

def p_chance(L, qlen, wlen):
    return min(1.0, qlen * wlen * 0.25 ** L)

rows, used = [], set()
for I in isl.itertuples():
    c = I.contig
    near = rna[(rna.contig == c) & (rna.end >= I.start - a.flank) & (rna.start <= I.end + a.flank)]
    for r in near.itertuples():
        used.add(r.Index)
        qs, qe = query(r)
        ws, we = I.start - a.pad, I.end + a.pad
        L, s0, e0, rep = longest_repeat(c, qs, qe, ws, we)
        # which island edge is the tRNA at, and is the repeat near the other one?
        mid = (r.start + r.end) / 2
        t_edge = "start" if abs(mid - I.start) <= abs(mid - I.end) else "end"
        if L:
            opp = I.end if t_edge == "start" else I.start
            d_opp = min(abs(s0 - opp), abs(e0 - opp))
            where = (f"near the opposite island edge ({d_opp} bp)" if d_opp <= a.edge
                     else f"inside/around island, {d_opp} bp from the opposite edge")
        else:
            d_opp, where = None, ""
        rows.append(dict(island=I.island, contig=c, island_start=I.start, island_end=I.end, tRNA=r.label,
                         tRNA_start=r.start, tRNA_end=r.end, strand=r.strand,
                         tRNA_position=("inside" if r.start >= I.start and r.end <= I.end else "flank"),
                         repeat_len=L, repeat_start=s0, repeat_end=e0, repeat_seq=rep,
                         repeat_location=where, P_chance=p_chance(L, qe - qs + 1, we - ws + 1) if L else 1.0,
                         candidate_att=L >= a.min_len))
R = pd.DataFrame(rows)

# ---- control: tRNAs not associated with islands, same window size
win = int(np.median([I.end - I.start + 1 + 2 * a.pad for I in isl.itertuples()]))
ctrl = []
for r in rna.drop(index=list(used)).itertuples():
    qs, qe = query(r)
    mid = (r.start + r.end) // 2
    L = longest_repeat(r.contig, qs, qe, mid - win // 2, mid + win // 2)[0]
    ctrl.append(L)
ctrl = np.array(ctrl)

out = a.out or a.islands_tsv.replace("_trna.tsv", "") + "_att.tsv"
R.to_csv(out, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 45)
print(R[["island", "tRNA", "tRNA_position", "repeat_len", "repeat_seq", "repeat_location", "P_chance"]]
      .to_string(index=False))
n_isl = R.groupby("island").candidate_att.any()
print(f"\nCandidate att repeats (>= {a.min_len} bp exact direct repeat of a tRNA 3' end): "
      f"{int(R.candidate_att.sum())} of {len(R)} tRNA/tmRNA genes; {int(n_isl.sum())} of {len(n_isl)} tRNA-associated islands")
if len(ctrl):
    print(f"Control: {len(ctrl)} tRNA/tmRNA genes away from islands, window {win} bp: longest repeat median "
          f"{np.median(ctrl):.0f} bp (range {ctrl.min()}-{ctrl.max()}); >= {a.min_len} bp in "
          f"{int((ctrl >= a.min_len).sum())} of {len(ctrl)}")
print(f"Results: {out}")
