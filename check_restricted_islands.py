import argparse, pandas as pd
ap = argparse.ArgumentParser()
for x in ("--restricted", "--faa", "--gff", "--islands"): ap.add_argument(x, required=True)
ap.add_argument("--flank", type=int, default=5000); ap.add_argument("--keep", default="")
a = ap.parse_args()
def fasta(p):
    d, cur = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): cur = l[1:].split()[0]; d[cur] = []
        elif cur: d[cur].append(l)
    return {k: "".join(v).rstrip("*").upper() for k, v in d.items()}
res, mine = fasta(a.restricted), fasta(a.faa)
if a.keep:
    k = pd.read_csv(a.keep, sep="\t"); keep = set(k.loc[k.status == "restricted_confirmed", "protein"])
    res = {q: s for q, s in res.items() if q in keep}
by_seq = {}
for pid, s in mine.items(): by_seq.setdefault(s, []).append(pid)
pos = {}
for l in open(a.gff):
    if l.startswith("##FASTA"): break
    f = l.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "CDS":
        att = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if "ID" in att: pos[att["ID"]] = (int(f[3]), int(f[4]), f[6], att.get("product", ""))
isl = pd.read_excel(a.islands) if a.islands.endswith((".xls", ".xlsx")) else pd.read_csv(a.islands, sep="\t")
merged = []
for s, e in sorted(set(zip(isl.iloc[:, 0].astype(int), isl.iloc[:, 1].astype(int)))):
    if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
    else: merged.append([s, e])
rows = []
for q, s in res.items():
    hit = by_seq.get(s, [])
    if not hit or hit[0] not in pos:
        rows.append(dict(protein=q, status="no_exact_match")); continue
    st, en, strand, prod = pos[hit[0]]
    best = min(((0 if (st <= b and en >= a0) else min(abs(st - b), abs(a0 - en)), a0, b) for a0, b in merged))
    status = "in_island" if best[0] == 0 else "flanking" if best[0] <= a.flank else "outside"
    rows.append(dict(protein=q, kdesi_locus=hit[0], start=st, end=en, strand=strand, status=status,
                     island=f"{best[1]}-{best[2]}", distance_bp=best[0], product=prod))
df = pd.DataFrame(rows); df.to_csv("restricted_vs_islands.tsv", sep="\t", index=False)
print(f"{len(merged)} island regions; {len(df)} restricted proteins"); print(df.status.value_counts().to_string())
print(df[df.status.isin(["in_island", "flanking"])].to_string(index=False))
