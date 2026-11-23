#!/usr/bin/env python3
"""Can kdesi (and each other genome) make methionine? Checks the methionine biosynthesis pathway in Prokka GFFs.

Steps (any one alternative per step is enough):
  1 homoserine         hom/thrA  homoserine dehydrogenase
  2 acylation          metA      homoserine O-succinyltransferase   | metX  homoserine O-acetyltransferase
  3 sulfur addition    metB+metC cystathionine gamma-synthase + beta-lyase (transsulfuration)
                       | metY/metZ O-acetyl/succinyl-homoserine sulfhydrylase (direct sulfhydrylation)
  4 methylation        metE      cobalamin-independent methionine synthase | metH cobalamin-dependent
  (support)            metF      methylenetetrahydrofolate reductase (supplies the methyl group)
Matches use Prokka gene names and product descriptions, so every hit is listed for inspection:
annotation of the cystathionine gamma-synthase family in particular can be ambiguous.

Usage:  python met_pathway.py --gff-dir ~/bdello_analysis/gff [--focus kdesi]
"""
import argparse, glob, os, re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gff-dir", required=True)
ap.add_argument("--focus", default="kdesi")
ap.add_argument("--out", default=None)
a = ap.parse_args()

GENES = {  # key: (gene-name regex, product regex)
    "hom":  (r"^(hom|thrA)\d*$", r"homoserine dehydrogenase|aspartokinase.*homoserine"),
    "metA": (r"^metA\d*$", r"homoserine O-succinyltransferase"),
    "metX": (r"^metX\d*$", r"homoserine O-acetyltransferase"),
    "metB": (r"^metB\d*$", r"cystathionine gamma-synthase"),
    "metC": (r"^(metC|malY|patB)\d*$", r"cystathionine beta-lyase"),
    "metYZ": (r"^(metY|metZ|cysD)\d*$", r"(O-acetyl|O-succinyl)-?L?-?homoserine (sulfhydrylase|\(thiol\)-lyase)|O-acetylhomoserine"),
    "metE": (r"^metE\d*$", r"5-methyltetrahydropteroyltriglutamate.{0,4}homocysteine methyltransferase|cobalamin-independent methionine synthase"),
    "metH": (r"^metH\d*$", r"^methionine synthase$|cobalamin-dependent methionine synthase"),
    "metF": (r"^metF\d*$", r"methylenetetrahydrofolate reductase"),
}

def parse(gff):
    hits = {k: [] for k in GENES}
    for line in open(gff):
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        name = at.get("gene", at.get("Name", ""))
        prod = at.get("product", "").replace("%2C", ",")
        for k, (gre, pre) in GENES.items():
            if re.search(gre, name, re.I) or re.search(pre, prod, re.I):
                hits[k].append(f"{at.get('locus_tag', at.get('ID', ''))}:{name or '-'}:{prod}")
    return hits

rows, detail = [], {}
for gff in sorted(glob.glob(os.path.join(os.path.expanduser(a.gff_dir), "*.gff"))):
    g = os.path.basename(gff)[:-4]
    h = parse(gff)
    detail[g] = h
    have = {k: len(v) > 0 for k, v in h.items()}
    steps = {"1_hom": have["hom"], "2_acyl": have["metA"] or have["metX"],
             "3_sulfur": (have["metB"] and have["metC"]) or have["metYZ"],
             "4_methyl": have["metE"] or have["metH"]}
    rows.append(dict(genome=g, **{k: len(v) for k, v in h.items()}, **steps,
                     complete=all(steps.values()), missing=",".join(s for s, ok in steps.items() if not ok)))
R = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(R[["genome", "hom", "metA", "metX", "metB", "metC", "metYZ", "metE", "metH", "metF", "complete", "missing"]]
      .to_string(index=False))
print(f"\nComplete methionine pathway: {int(R.complete.sum())} of {len(R)} genomes")
foc = [g for g in detail if a.focus in g]
for g in foc:
    print(f"\n{g}: annotated hits")
    for k, v in detail[g].items():
        print(f"  {k:6s} " + ("; ".join(v) if v else "none"))
out = a.out or os.path.join(os.path.dirname(os.path.expanduser(a.gff_dir.rstrip("/"))), "met_pathway.tsv")
R.to_csv(out, sep="\t", index=False)
print(f"\nResults: {out}")
