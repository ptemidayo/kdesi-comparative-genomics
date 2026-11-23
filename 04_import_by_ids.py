#!/usr/bin/env python3
"""Import SignalP-6 prediction files whatever their names, by matching protein IDs to genomes.

Each file in --src is assigned to the genome whose $OUT/faa/<genome>.faa shares its protein IDs;
it is imported to $OUT/signalp/<genome>/prediction_results.txt only if every protein of that
genome has exactly one prediction. Files that match no genome (e.g. SignalP run on a different
annotation) are reported.
Usage: python3 scripts/04_import_by_ids.py --out $OUT --src <folder with SignalP files>
"""
import argparse, glob, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import genome_ids

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--src", required=True)
a = ap.parse_args()

ids = genome_ids(a.out)
owner, nprot = {}, {}
for g in ids:
    n = 0
    for l in open(os.path.join(a.out, "faa", f"{g}.faa")):
        if l.startswith(">"):
            owner[l[1:].split()[0]] = g; n += 1
    nprot[g] = n

assigned, problems = {}, []
for f in sorted(glob.glob(os.path.join(a.src, "*"))):
    if not os.path.isfile(f):
        continue
    lines = [l.rstrip("\r\n") for l in open(f, errors="replace")]
    rows = [l for l in lines if l.strip() and not l.startswith("#")]
    if not rows or "\t" not in rows[0]:
        problems.append(f"{os.path.basename(f)}: not a SignalP prediction table"); continue
    pids = [r.split("\t")[0].split()[0] for r in rows]
    hits = {}
    for p in pids:
        g = owner.get(p)
        if g:
            hits[g] = hits.get(g, 0) + 1
    if not hits:
        problems.append(f"{os.path.basename(f)}: {len(pids)} proteins, none match any genome's protein IDs "
                        f"(example ID: {pids[0]})"); continue
    g, k = max(hits.items(), key=lambda x: x[1])
    if k < 0.99 * len(pids):
        problems.append(f"{os.path.basename(f)}: only {k} of {len(pids)} IDs belong to {g}"); continue
    if g in assigned:
        problems.append(f"{os.path.basename(f)}: {g} already matched by {os.path.basename(assigned[g][0])}"); continue
    if len(set(pids)) != nprot[g]:
        problems.append(f"{os.path.basename(f)}: {len(set(pids))} predictions but {g} has {nprot[g]} proteins"); continue
    assigned[g] = (f, lines)

for g, (f, lines) in assigned.items():
    os.makedirs(os.path.join(a.out, "signalp", g), exist_ok=True)
    with open(os.path.join(a.out, "signalp", g, "prediction_results.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  [ok] {g:45s} <- {os.path.basename(f)} ({nprot[g]} proteins)")
for p in problems:
    print(f"  [CHECK] {p}")
have = [g for g in ids if os.path.exists(os.path.join(a.out, "signalp", g, "prediction_results.txt"))]
missing = [g for g in ids if g not in have]
print(f"\nImported {len(assigned)} file(s); genomes with SignalP results: {len(have)} of {len(ids)}")
if missing:
    print("Still missing: " + ", ".join(missing))
