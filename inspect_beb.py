#!/usr/bin/env python3
"""Inspect branch-site hits for alignment artefacts.

For each orthogroup: reads the codon alignment codeml used (aln.phy) and the BEB sites (P >= 0.95)
from the better Model A run, then shows the protein alignment around the BEB sites for kdesi (S1),
its sister CKG001 (S10), Bd5 (S9) and HD100 (S2), and summarises signs of an artefact:
  - run:   longest run of consecutive BEB sites (real selection rarely hits 4+ adjacent codons)
  - span:  BEB sites lie within a short block (<= 15% of the alignment)
  - kd_id: kdesi identity to CKG001 inside the BEB block vs outside it
           (a block where kdesi alone diverges while the rest is conserved = frameshift/misalignment)
  - gaps:  kdesi gaps or stop codons near the block
Usage: python inspect_beb.py OG0000649 OG0000496 ...   (run inside .../selection/codeml_runs)
       python inspect_beb.py --all                    (every orthogroup with Model A runs)
"""
import os, re, sys, glob

CODE = {}
b = "TCAG"
aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(b):
    for j, y in enumerate(b):
        for k, z in enumerate(b):
            CODE[x + y + z] = aa[16 * i + 4 * j + k]
tr = lambda c: "-" if c == "---" else CODE.get(c.upper(), "X")
SHOW = [("S1", "kdesi"), ("S10", "CKG001"), ("S9", "Bd5"), ("S2", "HD100")]


def lnl(f):
    for line in open(f):
        if line.startswith("lnL"):
            return float(line.split()[-2])
    return None


def best_alt(og):
    runs = [(lnl(f), f) for f in glob.glob(f"{og}/A_alt_w*/out.txt")]
    runs = [r for r in runs if r[0] is not None]
    return max(runs)[1] if runs else None


def beb_sites(f):
    txt = open(f).read()
    i = txt.find("Bayes Empirical Bayes")
    if i < 0:
        return []
    sec = txt[i:txt.find("The grid", i) if txt.find("The grid", i) > 0 else None]
    return [int(m.group(1)) for m in re.finditer(r"^\s*(\d+) [A-Z\-*] +[0-9.]+\*", sec, re.M)]


def read_phy(f):
    lines = [l.split() for l in open(f) if l.strip()]
    return {l[0]: l[1] for l in lines[1:]}


def ident(a, b, idx):
    pairs = [(a[i], b[i]) for i in idx if i < len(a) and a[i] != "-" and b[i] != "-"]
    return sum(x == y for x, y in pairs) / len(pairs) if pairs else float("nan")


def longest_run(s):
    best = cur = 1 if s else 0
    for x, y in zip(s, s[1:]):
        cur = cur + 1 if y == x + 1 else 1
        best = max(best, cur)
    return best


ogs = sys.argv[1:]
if ogs == ["--all"]:
    ogs = sorted({p.split("/")[0] for p in glob.glob("*/A_alt_w0.5/out.txt")})
rows = []
for og in ogs:
    f = best_alt(og)
    if not f:
        print(f"{og}: no finished Model A run"); continue
    sites = sorted(beb_sites(f))
    seq = read_phy(os.path.join(os.path.dirname(f), "aln.phy"))
    prot = {k: "".join(tr(v[i:i + 3]) for i in range(0, len(v), 3)) for k, v in seq.items()}
    L = len(prot["S1"])
    if not sites:
        rows.append((og, L, 0, 0, 0, float("nan"), float("nan"), "", "no BEB sites"))
        continue
    lo, hi = min(sites) - 1, max(sites) - 1
    block = list(range(lo, hi + 1)); out = [i for i in range(L) if i < lo or i > hi]
    kin, kout = ident(prot["S1"], prot["S10"], block), ident(prot["S1"], prot["S10"], out)
    oin = ident(prot["S9"], prot["S10"], block)
    w = prot["S1"][max(0, lo - 10):hi + 11]
    flags = []
    run = longest_run(sites)
    if run >= 4: flags.append(f"{run} adjacent BEB codons")
    if (hi - lo + 1) <= 0.15 * L and len(sites) >= 5: flags.append("BEB sites in one short block")
    if kin < 0.5 and kout > 0.8: flags.append("kdesi diverges only inside block")
    if oin == oin and oin > 0.8 and kin < 0.5: flags.append("block conserved in Bd5 but not kdesi")
    if "*" in prot["S1"][:-1]: flags.append("internal stop in kdesi")
    if w.count("-") > 3: flags.append("kdesi gaps at block edge")
    if lo < 0.05 * L or hi > 0.95 * L: flags.append("block at alignment end")
    rows.append((og, L, len(sites), run, hi - lo + 1, kin, kout, f"{lo+1}-{hi+1}",
                 "; ".join(flags) if flags else "no artefact signs"))
    print(f"\n=== {og}  ({L} codons; BEB95 sites {len(sites)}: {sites[0]}-{sites[-1]})")
    a, z = max(0, lo - 10), min(L, hi + 11)
    for s, name in SHOW:
        if s in prot:
            print(f"  {name:7s} {a+1:5d} {prot[s][a:z]}")
    mark = "".join("^" if (i + 1) in sites else " " for i in range(a, z))
    print(f"  {'BEB':7s} {'':5s} {mark}")

print("\northogroup\tcodons\tBEB95\tmax_adjacent\tblock_len\tkdesi_vs_CKG001_id_in_block\tid_outside\tblock\tassessment")
for r in rows:
    print("\t".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r))
