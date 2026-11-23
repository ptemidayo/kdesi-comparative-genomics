#!/usr/bin/env python3
"""Pick one genome per genomospecies for the codeml branch-site analysis.

Rules (per genomospecies, outgroups and exclude_from_stats genomes skipped):
  1. the foreground genome represents its own genomospecies;
  2. else the type strain (type_strain == yes);
  3. else the most complete genome by CheckM2 (highest completeness, then lowest contamination),
     then the most contiguous assembly (fewest contigs, then largest total length).
Quality columns are read from --metadata (the curated table) when given, else from $OUT/metadata.tsv.
Prints the strain list (foreground first) on stdout; writes a table explaining each choice.
"""
import argparse, os, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--foreground", required=True)
ap.add_argument("--metadata", default="", help="curated table with checkm2_completeness/contamination")
ap.add_argument("--exclude", default="", help="genome IDs to leave out of the selection analysis (space/comma separated)")
a = ap.parse_args()

meta = pd.read_csv(os.path.join(a.out, "metadata.tsv"), sep="\t", dtype=str).fillna("")
if a.metadata and os.path.exists(a.metadata):
    cur = pd.read_csv(a.metadata, sep="\t", dtype=str).fillna("").set_index("genome")
    for c in ("checkm2_completeness", "checkm2_contamination", "type_strain", "exclude_from_stats", "outgroup"):
        if c in cur.columns:
            meta[c] = meta.genome.map(cur[c]).fillna("")
yes = lambda c: meta[c].str.lower().eq("yes") if c in meta.columns else False
drop = {x for x in a.exclude.replace(",", " ").split() if x}
unknown = drop - set(meta.genome)
if unknown:
    sys.exit(f"--exclude names not in metadata.tsv: {', '.join(sorted(unknown))}")
pool = meta[~(yes("outgroup") | yes("exclude_from_stats")) & ~meta.genome.isin(drop)].copy()
if drop:
    print(f"Excluded from the selection analysis by CODEML_EXCLUDE: {', '.join(sorted(drop))}", file=sys.stderr)
if a.foreground not in set(pool.genome):
    sys.exit(f"Foreground {a.foreground} is not an ingroup genome in metadata.tsv")


def assembly_stats(g):
    n, L = 0, 0
    for line in open(os.path.join(a.out, "fna", f"{g}.fna")):
        if line.startswith(">"):
            n += 1
        else:
            L += len(line.strip())
    return n, L


pool[["contigs", "length"]] = [assembly_stats(g) for g in pool.genome]
num = lambda c, default: pd.to_numeric(pool[c], errors="coerce").fillna(default) if c in pool.columns \
    else pd.Series(default, index=pool.index)
pool["completeness"] = num("checkm2_completeness", 0.0)
pool["contamination"] = num("checkm2_contamination", 100.0)
fg_gs = pool.loc[pool.genome == a.foreground, "genomospecies"].iloc[0]
rows = []
for gs, sub in pool.groupby("genomospecies", sort=True):
    if gs == fg_gs:
        pick, why = a.foreground, "foreground"
    elif "type_strain" in sub.columns and (sub.type_strain.str.lower() == "yes").any():
        pick, why = sub[sub.type_strain.str.lower() == "yes"].genome.iloc[0], "type strain"
    elif len(sub) == 1:
        pick, why = sub.genome.iloc[0], "only member"
    else:
        best = sub.sort_values(["completeness", "contamination", "contigs", "length"],
                               ascending=[False, True, True, False]).iloc[0]
        pick, why = best.genome, (f"most complete of {len(sub)} (CheckM2 {best.completeness:.2f}% / "
                                  f"{best.contamination:.2f}%, {best.contigs} contigs)")
    r = sub[sub.genome == pick].iloc[0]
    rows.append(dict(genomospecies=gs, genome=pick, strain=r.get("strain", ""), reason=why,
                     contigs=r.contigs, length_bp=r.length, completeness=r.completeness,
                     contamination=r.contamination, members=len(sub)))
rep = pd.DataFrame(rows)
os.makedirs(os.path.join(a.out, "selection"), exist_ok=True)
rep.to_csv(os.path.join(a.out, "selection", "codeml_representatives.tsv"), sep="\t", index=False)
for r in rep.itertuples():
    print(f"  {r.genomospecies:5s} {r.strain:12s} {r.reason}", file=sys.stderr)
skipped = sorted(set(meta.genome) - set(pool.genome))
print(f"{len(rep)} representatives; skipped (outgroup / excluded): "
      f"{', '.join(meta.set_index('genome').loc[skipped, 'strain']) if skipped else 'none'}", file=sys.stderr)
order = [a.foreground] + [g for g in rep.genome if g != a.foreground]
print(" ".join(order))
