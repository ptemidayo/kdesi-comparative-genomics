#!/usr/bin/env python3
"""Build an ANI matrix from `skani triangle -E` output and assign genomospecies
(single-linkage clusters at >= cutoff % ANI)."""
import argparse, os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--edges", required=True)
ap.add_argument("--ids", required=True)
ap.add_argument("--cutoff", type=float, default=95.0)
ap.add_argument("--matrix", required=True)
ap.add_argument("--metadata", required=True)
ap.add_argument("--curated", default="", help="curated genomospecies table to use instead of skani clusters")
a = ap.parse_args()

ids = [l.strip() for l in open(a.ids) if l.strip()]
gid = lambda p: os.path.basename(p).rsplit(".", 1)[0]
e = pd.read_csv(a.edges, sep="\t")
e["r"], e["q"] = e["Ref_file"].map(gid), e["Query_file"].map(gid)

M = pd.DataFrame(0.0, index=ids, columns=ids)
for i in ids:
    M.loc[i, i] = 100.0
for r, q, v in zip(e["r"], e["q"], e["ANI"]):
    M.loc[r, q] = M.loc[q, r] = v
M.round(3).to_csv(a.matrix, sep="\t")

# union-find single linkage
parent = {i: i for i in ids}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for r, q, v in zip(e["r"], e["q"], e["ANI"]):
    if v >= a.cutoff:
        parent[find(r)] = find(q)
clusters = {}
for i in ids:
    clusters.setdefault(find(i), []).append(i)
order = sorted(clusters.values(), key=lambda c: (-len(c), sorted(c)[0]))
label = {g: f"GS{n+1}" for n, c in enumerate(order) for g in c}
print(f"{len(order)} genomospecies at >= {a.cutoff}% ANI")
for n, c in enumerate(order):
    print(f"  GS{n+1} ({len(c)}): {', '.join(sorted(c))}")

sk = pd.DataFrame({"genome": ids, "genomospecies": [label[i] for i in ids],
                   "group": [label[i] for i in ids]})
sk.to_csv(os.path.join(os.path.dirname(a.matrix), "skani_genomospecies.tsv"), sep="\t", index=False)

if a.curated:
    from sklearn.metrics import adjusted_rand_score
    cur = pd.read_csv(a.curated, sep="\t", dtype=str)
    miss = sorted(set(ids) - set(cur.genome)); extra = sorted(set(cur.genome) - set(ids))
    if miss or extra:
        raise SystemExit(f"Curated table does not match genome IDs.\n  missing: {miss}\n  not a genome: {extra}")
    cur.to_csv(a.metadata, sep="\t", index=False)
    is_out = cur["outgroup"].astype(str).str.lower() == "yes" if "outgroup" in cur.columns else False
    ing = cur[~is_out].set_index("genome")
    ari = adjusted_rand_score(ing.genomospecies, [label[g] for g in ing.index])
    print(f"\nUsing curated genomospecies from {a.curated}: "
          f"{ing.genomospecies.nunique()} genomospecies in {len(ing)} ingroup genomes")
    print(f"Cross-check, skani {a.cutoff}% single-linkage vs curated (ingroup): ARI = {ari:.3f}")
    for g in sorted(ing.genomospecies.unique()):
        mem = ing.index[ing.genomospecies == g]
        sl = {label[x] for x in mem}
        if len(sl) > 1 or sum(label[x] in sl for x in ing.index) != len(mem):
            print(f"  differs: {g} -> skani clusters {sorted(sl)}")
    print(f"Wrote {a.metadata}")
elif os.path.exists(a.metadata):
    print(f"{a.metadata} exists — not overwritten (delete it to regenerate).")
else:
    sk.to_csv(a.metadata, sep="\t", index=False)
    print(f"Wrote {a.metadata}")
