#!/usr/bin/env python3
"""Do kdesi CheR and CheD proteins diverge more than expected from those of its relatives?

Claim to test (old draft): "CheR2 and CheD showed sequence divergence".
A protein is only 'divergent' if its identity to the matching protein in a relative is unusually low
compared with ALL kdesi proteins against the same relative. For each comparison genome the script:
  1. finds CheR (Pfam PF01739 CheR / PF03705 CheR_N) and CheD (PF03975) proteins in every genome
     from the Pfam domain table, and reports copy numbers across the genus;
  2. runs DIAMOND blastp of every kdesi protein against the comparison proteome (best hit, >= 50% query
     coverage) to get kdesi's genome-wide identity distribution to that relative;
  3. reports each kdesi CheR/CheD protein's identity to its best hit, whether the hit is reciprocal, and
     its percentile within the genome-wide distribution (below the 5th percentile = divergent).

Usage (in the bdello environment):
  python che_divergence.py --out ~/bdello_analysis
"""
import argparse, glob, os, subprocess, tempfile
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True, help="pipeline OUT folder (contains faa/ and pfam/results/)")
ap.add_argument("--foreground", default="kdesi")
ap.add_argument("--compare", default="Bd5,CKG001,HAGR004,HD100",
                help="comma-separated strain names (matched within genome IDs)")
ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args()
out = os.path.expanduser(a.out)
CHE = {"PF01739": "CheR", "PF03705": "CheR", "PF03975": "CheD"}

dom = pd.read_csv(os.path.join(out, "pfam", "results", "protein_domains_filtered.tsv"), sep="\t")
dom["acc"] = dom.pfam_acc.astype(str).str.split(".").str[0]
che = dom[dom.acc.isin(CHE)].copy()
che["family"] = che.acc.map(CHE)
che = che.drop_duplicates(["genome", "protein", "family"])

genomes = sorted({os.path.basename(f)[:-4] for f in glob.glob(os.path.join(out, "faa", "*.faa"))})
def gid(name):
    hits = [g for g in genomes if g.endswith("_" + name) or g == name]
    if not hits:
        hits = [g for g in genomes if name in g]
    if not hits:
        raise SystemExit(f"no genome matching '{name}' in {out}/faa")
    return hits[0]
fg = gid(a.foreground)

# ---- copy numbers across the genus
cn = che.groupby(["genome", "family"]).protein.nunique().unstack(fill_value=0)
cn = cn.reindex(genomes, fill_value=0)
print("CheR / CheD copies per genome (Pfam):")
print(f"  kdesi: CheR {cn.loc[fg].get('CheR', 0)}, CheD {cn.loc[fg].get('CheD', 0)}")
for f in ("CheR", "CheD"):
    if f in cn:
        v = cn[f]
        print(f"  {f}: genus median {v.median():.0f}, range {v.min()}-{v.max()}; "
              f"genomes with fewer than kdesi: {(v < v[fg]).sum()} of {len(v)}")

def read_faa(p):
    seqs, cur = {}, None
    for l in open(p):
        if l.startswith(">"):
            cur = l[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items()}

def best_hits(q_faa, db_faa, tmp):
    db = os.path.join(tmp, os.path.basename(db_faa) + ".dmnd")
    subprocess.run(["diamond", "makedb", "--in", db_faa, "-d", db, "--quiet"], check=True)
    o = os.path.join(tmp, "hits.tsv")
    subprocess.run(["diamond", "blastp", "-q", q_faa, "-d", db, "-o", o, "--quiet", "--more-sensitive",
                    "--max-target-seqs", "1", "--evalue", "1e-5", "--threads", str(a.threads),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore"], check=True)
    h = pd.read_csv(o, sep="\t", names=["q", "s", "pident", "qcov", "bits"])
    h = h.sort_values("bits", ascending=False).drop_duplicates("q")
    return h[h.qcov >= 50]

fg_faa = os.path.join(out, "faa", f"{fg}.faa")
mine = che[che.genome == fg][["protein", "family"]].drop_duplicates("protein")
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for name in [x.strip() for x in a.compare.split(",") if x.strip()]:
        g = gid(name)
        g_faa = os.path.join(out, "faa", f"{g}.faa")
        fwd = best_hits(fg_faa, g_faa, tmp)
        rev = best_hits(g_faa, fg_faa, tmp)
        rbh = dict(zip(rev.q, rev.s))
        dist = fwd.pident
        theirs = set(che[che.genome == g].protein)
        for p, fam in mine.itertuples(index=False):
            r = fwd[fwd.q == p]
            if r.empty:
                rows.append(dict(kdesi_protein=p, family=fam, relative=name, best_hit="none", identity=None,
                                 reciprocal=False, hit_is_same_family=False, percentile=None)); continue
            r = r.iloc[0]
            pct = 100 * (dist < r.pident).mean()
            rows.append(dict(kdesi_protein=p, family=fam, relative=name, best_hit=r.s, identity=round(r.pident, 1),
                             reciprocal=rbh.get(r.s) == p, hit_is_same_family=r.s in theirs,
                             percentile=round(pct, 1)))
        print(f"\nkdesi vs {name}: {len(fwd)} kdesi proteins with a hit; genome-wide identity median "
              f"{dist.median():.1f}% (5th percentile {dist.quantile(0.05):.1f}%)")

R = pd.DataFrame(rows)
R["divergent"] = R.percentile.notna() & (R.percentile < 5)
pd.set_option("display.width", 200)
print("\nkdesi CheR / CheD proteins against each relative:")
print(R.to_string(index=False))
outp = os.path.join(out, "che_divergence.tsv"); R.to_csv(outp, sep="\t", index=False)
close = R[R.relative.isin(["Bd5", "CKG001", "HAGR004"])]
n_div = close.groupby("kdesi_protein").divergent.all().sum()
print(f"\nkdesi CheR/CheD proteins below the 5th percentile against all of its closest relatives: "
      f"{n_div} of {close.kdesi_protein.nunique()}")
print("Keep the 'divergence' claim only for proteins that are divergent against all close relatives;"
      " otherwise drop it.")
print(f"Results: {outp}")
