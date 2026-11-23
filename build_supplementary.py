#!/usr/bin/env python3
"""Assemble Supplementary Tables S1-S13 into one Excel workbook.

Each table is filled from the pipeline output files on this computer. Files are found automatically under
~/bdello_analysis, ~/bdello_ncbi, ~/bdello_pipeline, ~/bdello_annot and your Downloads folder; any file can
be given explicitly with --file KEY=PATH (keys are listed in the README sheet and in the printout).
A table whose files are not found gets a sheet that says what is missing, so the workbook is always built.

Usage:
  python build_supplementary.py --downloads "/mnt/c/Users/<you>/Downloads" [--file S2_ani=/path/ANIb.tab ...]
"""
import argparse, glob, io, os, re, zipfile
import xml.etree.ElementTree as ET
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

H = os.path.expanduser
ap = argparse.ArgumentParser()
ap.add_argument("--downloads", default="")
ap.add_argument("--file", action="append", default=[], help="KEY=PATH override")
ap.add_argument("--out", default="Supplementary_Tables_S1-S13.xlsx")
a = ap.parse_args()
ROOTS = [H("~/bdello_analysis"), H("~/bdello_ncbi"), H("~/bdello_pipeline"), H("~/bdello_annot")]
if a.downloads:
    ROOTS.append(a.downloads)
OVR = dict(x.split("=", 1) for x in a.file)

def find(key, patterns):
    if key in OVR:
        return OVR[key] if os.path.exists(OVR[key]) else None
    for root in ROOTS:
        for p in patterns:
            hits = sorted(glob.glob(os.path.join(root, "**", p), recursive=True),
                          key=lambda f: ("selection_previous" in f or "codeml_25taxa" in f or "test" in f.lower(),
                                         -os.path.getmtime(f)))
            if hits:
                return hits[0]
    return None

def read_docx_table(path):
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    root = ET.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
    tbl = root.find(f".//{W}tbl")
    rows = [["".join(t.text or "" for t in tc.iter(f"{W}t")).strip() for tc in tr.findall(f"{W}tc")]
            for tr in tbl.findall(f"{W}tr")]
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    return pd.DataFrame(rows[1:], columns=rows[0])

def ezaai_matrix(folder):
    """Combine EzAAI 'calculate' outputs (one file per genome pair, named A__B.tsv) into an AAI matrix,
    limited to the genomes in metadata_genomospecies.tsv when that file is found."""
    vals = {}
    for f in glob.glob(os.path.join(folder, "**", "*__*.tsv"), recursive=True):
        a_, b_ = os.path.basename(f)[:-4].split("__", 1)
        try:
            t = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        col = [c for c in t.columns if c.strip().upper() == "AAI"] or [c for c in t.columns if "AAI" in c.upper()]
        if not col or t.empty:
            continue
        v = float(str(t[col[0]].iloc[0]).rstrip("%"))
        vals[(a_, b_)] = vals[(b_, a_)] = v
    names = sorted({x for p in vals for x in p})
    meta = find("S1_metadata", ["metadata_genomospecies.tsv"])
    if meta:
        keep = set(pd.read_csv(meta, sep="\t").genome)
        names = [n for n in names if n in keep] or names
    m = pd.DataFrame(index=names, columns=names, dtype=float)
    for n in names:
        m.loc[n, n] = 100.0
    for (x, y), v in vals.items():
        if x in m.index and y in m.index:
            m.loc[x, y] = v
    print(f"EzAAI: combined {len(vals)//2} pairs into a {len(names)} x {len(names)} AAI matrix")
    return m.reset_index().rename(columns={"index": "genome"})

def read(path):
    if os.path.isdir(path):
        return ezaai_matrix(path)
    if path.endswith(".docx"):
        return read_docx_table(path)
    if path.endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    if os.path.basename(path).startswith("ggdc_results"):
        return read_ggdc(path)
    sep = "," if path.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, comment=None)

found, missing = [], []
def block(key, label, patterns, transform=None):
    p = find(key, patterns)
    if not p:
        missing.append((key, label, patterns)); return (label, None, key)
    try:
        df = read(p)
        if transform:
            df = transform(df)
    except Exception as e:
        missing.append((key, f"{label} (could not read {p}: {e})", patterns)); return (label, None, key)
    found.append((key, label, p)); return (label, df, p)

# ---------------- helpers for merged columns
def gff_products():
    g = find("kdesi_gff", ["Bdellovibrio_bacteriovorus_kdesi.gff"])
    prod = {}
    if g:
        for line in open(g):
            if line.startswith("##FASTA"):
                break
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] == "CDS":
                at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                prod[at.get("locus_tag", at.get("ID", ""))] = at.get("product", "").replace("%2C", ",")
    return prod

def faa_products(path):
    prod = {}
    if path:
        for l in open(path):
            if l.startswith(">"):
                x = l[1:].strip().split(None, 1); prod[x[0]] = x[1] if len(x) > 1 else ""
    return prod

# assessment of the seven significant orthogroups (inspection of BEB sites and alignments, Nov 2026)
ASSESS = {
    "OG0000360": "robust: 4 scattered BEB sites; kdesi 88% identical to CKG001 in the region",
    "OG0001892": "robust: 3 scattered BEB sites; 90% identity in the region",
    "OG0001943": "robust: 2 widely separated BEB sites; 93% identity in the region",
    "OG0000426": "weak: single BEB site; q = 0.047",
    "OG0000919": "weak: no BEB site > 0.95; maximum branch dS 2.96, close to the saturation cut-off",
    "OG0000505": "artefact: BEB sites in a C-terminal frameshift of the kdesi sequence",
    "OG0001909": "artefact: BEB sites flank a 190-codon region absent from kdesi",
}
def s13(df):
    prod = gff_products()
    df.insert(df.columns.get_loc("foreground_protein") + 1, "kdesi_product",
              df.foreground_protein.map(prod).astype(object).fillna(""))
    df["assessment"] = df.orthogroup.map(ASSESS).fillna("")
    key = df.positive_selection_foreground.astype(str).str.lower() == "true"
    return pd.concat([df[key].sort_values("q_branchsite_BH"), df[~key].sort_values("p_branchsite")])

def s9(df):
    q = find("restricted_faa", ["kdesi_restricted_all.faa"])
    df.insert(1, "product", df.protein.map(faa_products(q)).fillna(""))
    return df

def kdesi_rows(df):
    m = df.astype(str).apply(lambda c: c.str.contains("kdesi|G17_", regex=True)).any(axis=1)
    return df[m]

def tidy_matrix(df):
    df = df.rename(columns={df.columns[0]: "genome"})
    fix = {"SSB219315": "SSB218315", "SKB1291224": "SKB1291214", "Bdellovibrio sp 110": "Bdellovibrio sp HCB-110",
           "Kdesi": "kdesi"}
    def f(x):
        x = str(x)
        for k, v in fix.items():
            x = x.replace(k, v)
        return x
    df.columns = [f(c) for c in df.columns]
    df["genome"] = df["genome"].map(f)
    num = df.columns[1:]
    df[num] = df[num].apply(pd.to_numeric, errors="coerce").round(2)
    return df

def ggdc_kdesi(df):
    """GGDC web results (two header rows): keep query, reference and formula 2 columns."""
    raw = df
    if "Query genome" not in raw.columns:
        raw.columns = raw.iloc[0]; raw = raw.iloc[1:]
    cols = list(raw.columns)
    out = pd.DataFrame({
        "query genome": raw.iloc[:, 0], "reference genome": raw.iloc[:, 1],
        "dDDH formula 2 (%)": pd.to_numeric(raw.iloc[:, 6], errors="coerce"),
        "model 95% CI": raw.iloc[:, 7],
        "probability dDDH >= 70% (%)": pd.to_numeric(raw.iloc[:, 9], errors="coerce"),
        "G+C difference (%)": pd.to_numeric(raw.iloc[:, 14], errors="coerce") if len(cols) > 14 else None})
    return out.sort_values("dDDH formula 2 (%)", ascending=False)

def read_ggdc(path):
    return pd.read_csv(path, header=None, skiprows=1)

# ---------------- table definitions
S5_ROWS = pd.DataFrame([
    ["Curated genomes", 38, 26, 67.6, "39-165", 38, 0.50],
    ["Curated + 8 public isolates", 46, 31, 80.5, "47-181", 39, 0.50],
    ["Curated + isolates + 13 MAGs", 59, 42, 128.5, "74-273", 33, 0.56],
    ["Conservative: added genomes with >=70% ANI to a curated genome", 48, 31, 68.7, "44-144", 45, 0.45],
    ["Curated, without HCB337", 37, 25, 62.2, "36-151", 40, 0.48],
    ["Without HCB337, + isolates + MAGs", 58, 41, 122.2, "71-259", 34, 0.55],
], columns=["dataset", "genomes", "genomospecies_observed", "Chao1", "Chao1_95pct_CI", "percent_sampled",
            "probability_next_genome_new"])

S10_ROWS = pd.DataFrame([
    ["G17_01061", "restriction endonuclease", 222, "+", "1,039,648-1,041,791 (pair 1)", "4 bp with G17_01062",
     "IPR059604 (NspV)", "Type II restriction enzyme NspV (UniProt P35677)",
     "Csp68KII 62.8% / 93%; NspV 61.8% / 94%", "TTCGAA", "yes (putative type II)"],
    ["G17_01062", "DNA methyltransferase", 492, "+", "1,039,648-1,041,791 (pair 1)", "4 bp with G17_01061",
     "IPR029063, IPR041698, IPR050953 (generic SAM-dependent / N4/N6-adenine methyltransferase)", "",
     "M.XfaTI 42.3% / 94%; M.Csp68KII 42.1% / 94%; M.NspV 44.5% / 84%", "TTCGAA", "yes (putative type II)"],
    ["G17_02633", "restriction endonuclease", 354, "-", "2,657,766-2,660,355 (pair 2)", "11 bp with G17_02632",
     "IPR019054 (AccI)", "Type II restriction enzyme AccI (UniProt P24546)",
     "AccI 46.9% / 97%", "GTMKAC", "no"],
    ["G17_02632", "DNA methyltransferase", 511, "-", "2,657,766-2,660,355 (pair 2)", "11 bp with G17_02633",
     "IPR002052, IPR003356 (N6-adenine DNA methylase)", "",
     "M.AccI 38.6% / 92%", "GTMKAC", "yes (putative type II)"],
], columns=["locus_tag", "role", "length_aa", "strand", "location", "reading-frame overlap",
            "InterProScan entries (EMBL-EBI web service)", "Prokka name",
            "closest REBASE Gold Standard enzymes (identity / query coverage)",
            "recognition site of the closest characterized enzymes", "listed for kdesi (CP102930) in REBASE"])

TABLES = [
 ("S1", "S1 Genomes", "Table S1. The 41 genomes analysed: accession, genomospecies, assembly statistics, CheckM2 "
  "completeness and contamination, and isolation source.",
  [lambda: block("S1_docx", "Genome table (from Table S1 Word file)", ["Table_S1_genomes_41*.docx"]),
   lambda: block("S1_metadata", "Pipeline metadata (genomospecies, CheckM2, isolation source)", ["metadata_genomospecies.tsv"])]),
 ("S2", "S2 Pairwise relatedness", "Table S2. Pairwise genome relatedness among the 41 genomes: ANIb (pyani, mean of both directions), AAI "
  "(EzAAI) and dDDH (GBDP formula d4 = GGDC formula 2, TYGS; all 820 pairs), followed by the dDDH of kdesi against "
  "the other 40 genomes with model confidence intervals (GGDC).",
  [lambda: block("S2_ani", "ANIb (%)", ["*ANIb_percentage_identity*", "*ANIb*.tab", "*anib*matrix*.tsv"], tidy_matrix),
   lambda: block("S2_aai", "AAI (%)", ["*AAI*matrix*", "*aai*matrix*.tsv", "*ezaai*.tsv"], tidy_matrix),
   lambda: block("S2_ddh", "dDDH (%, formula d4 / GGDC formula 2; TYGS, all pairs)", ["TYGS_dDDH_d4_all_pairs_matrix.tsv", "*TYGS*dDDH*"], tidy_matrix),
   lambda: block("S2_ddh_kdesi", "dDDH of kdesi against the other 40 genomes (GGDC formula 2)",
                 ["ggdc_results_26-10-08-22-29*.csv"], ggdc_kdesi)]),
 ("S3", "S3 NCBI screen", "Table S3. Classification of the 194 GenBank Bdellovibrio assemblies: relationship to the "
  "curated genomes, GTDB taxonomy, CheckM2 quality and inclusion decision.",
  [lambda: block("S3_screen", "Screen of 194 assemblies (skani, GTDB)", ["ncbi_screen*_genomes_gtdb.tsv", "ncbi_screen*_genomes.tsv"]),
   lambda: block("S3_checkm2", "CheckM2 quality of the candidate genomes", ["checkm2_candidates*.tsv", "quality_report.tsv"])]),
 ("S4", "S4 Assignment evidence", "Table S4. Evidence for the genomospecies assignment of the 21 added public genomes.",
  [lambda: block("S4", "Assignment evidence", ["assignment_evidence.tsv"])]),
 ("S5", "S5 Richness", "Table S5. Genomospecies richness (bias-corrected Chao1 with log-normal 95% CI, share of "
  "estimated richness sampled, and probability that the next genome belongs to an unsampled genomospecies).",
  [lambda: ("Richness estimates (from fig4.py, run on 1 Nov 2026)", S5_ROWS, "values from fig4.py output")]),
 ("S6", "S6 Species-level traits", "Table S6. Share of variance between genomospecies (eta squared) and permutation "
  "tests (9,999 permutations, Benjamini-Hochberg) for 57 functional traits.",
  [lambda: block("S6", "Species-level traits", ["species_level_traits.tsv"])]),
 ("S7", "S7 kdesi vs genus", "Table S7. Functional profile of kdesi compared with the 36 other ingroup genomes: Pfam "
  "categories, chemoreceptors, secretome, CheR/CheD divergence and the methionine pathway.",
  [lambda: block("S7_sizeadj", "Pfam categories, size-adjusted test", ["foreground_vs_others_size_adjusted.tsv"]),
   lambda: block("S7_ch", "Pfam categories, Crawford-Howell test (percent of proteome)", ["foreground_vs_others_crawford_howell.tsv"]),
   lambda: block("S7_mcp", "Chemoreceptor (MCP) domain profile per genome", ["mcp_domain_profile.tsv"]),
   lambda: block("S7_secretome", "Secretome tests", ["foreground_secretome_tests.tsv"]),
   lambda: block("S7_che", "CheR/CheD divergence", ["che_divergence.tsv"]),
   lambda: block("S7_met", "Methionine biosynthesis pathway per genome", ["met_pathway.tsv"])]),
 ("S8", "S8 Genomic islands", "Table S8. Genomic islands of kdesi and HD100 (IslandViewer 4), with gene content, tRNA "
  "and tmRNA associations and candidate att repeats.",
  [lambda: block("S8_ktrna", "kdesi islands and tRNA/tmRNA genes", ["kdesi_trna.tsv"]),
   lambda: block("S8_htrna", "HD100 islands and tRNA/tmRNA genes", ["HD100_trna.tsv"]),
   lambda: block("S8_katt", "kdesi candidate att repeats", ["kdesi_att.tsv"]),
   lambda: block("S8_hatt", "HD100 candidate att repeats", ["HD100_att.tsv"]),
   lambda: block("S8_genes", "Island genes", ["island_genes_all.tsv"])]),
 ("S9", "S9 Restricted proteins", "Table S9. The 62 kdesi proteins without BLASTp homologues in the other genomes, "
  "with tblastn results against all 40 genome sequences.",
  [lambda: block("S9", "Restricted proteins and tblastn results", ["restricted_tblastn_summary.tsv"], s9),
   lambda: block("S9_islands", "Restricted proteins and genomic islands", ["restricted_vs_islands.tsv"])]),
 ("S10", "S10 RM systems", "Table S10. The two type II restriction-modification pairs unique to kdesi: gene "
  "organization, InterProScan entries (EMBL-EBI web service), closest REBASE Gold Standard enzymes (BLASTp, accessed "
  "17 October 2026) and DefenseFinder systems. Recognition sites are those of the closest characterized enzymes and "
  "were not determined for kdesi.",
  [lambda: ("RM proteins: InterProScan and REBASE (values as reported in the manuscript)", S10_ROWS,
            "InterProScan (EMBL-EBI) and REBASE web searches, accessed 17 October 2026"),
   lambda: block("S10", "DefenseFinder systems in kdesi", ["all_systems.tsv"], kdesi_rows)]),
 ("S11", "S11 Pfam categories", "Table S11. The 166 Pfam accessions and their functional categories.",
  [lambda: block("S11", "Pfam accession to category map", ["pfam_categories.tsv"])]),
 ("S12", "S12 Codon masking", "Table S12. Masking of kdesi codons in each orthogroup before the branch-site analysis.",
  [lambda: block("S12", "Foreground masking", ["foreground_masking.tsv"])]),
 ("S13", "S13 Branch-site tests", "Table S13. Branch-site test results for the 1,454 orthogroups tested on the kdesi "
  "lineage: M0 omega, screening and full-model likelihood-ratio statistics, BH-corrected q-values, foreground "
  "omega2, Bayes empirical Bayes sites and the assessment of each significant result.",
  [lambda: block("S13", "Branch-site results (significant orthogroups first)", ["codeml_branchsite_results.tsv"], s13)]),
]

# ---------------- write workbook
FONT = "Arial"
hdr_fill = PatternFill("solid", fgColor="E7E6E1")
wb = Workbook(); readme = wb.active; readme.title = "README"
readme["A1"] = "Supplementary Tables S1-S13"; readme["A1"].font = Font(name=FONT, bold=True, size=12)
readme["A2"] = ("kdesi comparative genomics manuscript. One sheet per table; each sheet starts with its legend, "
                "followed by one or more blocks with their source file.")
readme["A2"].font = Font(name=FONT, italic=True, size=9)
readme.append([]); readme.append(["Table", "Sheet", "Legend", "Status"])
for c in readme[4]:
    c.font = Font(name=FONT, bold=True); c.fill = hdr_fill

for tid, sheet, legend, parts in TABLES:
    ws = wb.create_sheet(sheet)
    ws["A1"] = legend; ws["A1"].font = Font(name=FONT, bold=True, size=11)
    ws["A1"].alignment = Alignment(wrap_text=False)
    r = 3; status = []; widths = {}; first_hdr = None
    for part in parts:
        label, df, src = part()
        ws.cell(r, 1, label).font = Font(name=FONT, bold=True, size=10)
        ws.cell(r + 1, 1, f"Source: {src}" if df is not None else
                f"MISSING - provide with --file {src}=PATH").font = Font(name=FONT, italic=True, size=8,
                                                                         color="000000" if df is not None else "C00000")
        r += 2
        if df is None:
            status.append(f"{label}: missing"); r += 1; continue
        df = df.copy()
        for j, col in enumerate(df.columns, 1):
            c = ws.cell(r, j, str(col)); c.font = Font(name=FONT, bold=True, size=9); c.fill = hdr_fill
            widths[j] = max(widths.get(j, 8), min(len(str(col)) + 2, 45))
        first_hdr = first_hdr or r
        for row in df.itertuples(index=False):
            r += 1
            for j, v in enumerate(row, 1):
                if isinstance(v, float):
                    v = None if pd.isna(v) else (round(v, 6) if abs(v) >= 1e-3 or v == 0 else float(f"{v:.3e}"))
                c = ws.cell(r, j, v); c.font = Font(name=FONT, size=9)
                widths[j] = max(widths.get(j, 8), min(len(str(v)) + 2, 45))
        status.append(f"{label}: {len(df)} rows"); r += 3
    for j, w in widths.items():
        ws.column_dimensions[get_column_letter(j)].width = w
    if first_hdr and len(parts) == 1:
        ws.freeze_panes = ws.cell(first_hdr + 1, 1)
    readme.append([tid, sheet, legend, "; ".join(status)])
    for c in readme[readme.max_row]:
        c.font = Font(name=FONT, size=9); c.alignment = Alignment(wrap_text=True, vertical="top")

readme.append([]); readme.append(["Supplementary figures (supplied as separate files)"])
readme.cell(readme.max_row, 1).font = Font(name=FONT, bold=True)
for f in ["Figure S1. Reference-based core-SNP phylogeny of the 41 genomes (REALPHY).",
          "Figure S2. skani ANI among the 38 curated genomes and the 21 added public genomes (skani_groups_heatmap.pdf).",
          "Figure S3. Principal component analysis of the functional profiles of the 37 ingroup genomes."]:
    readme.append([f]); readme.cell(readme.max_row, 1).font = Font(name=FONT, size=9)
for col, w in zip("ABCD", (8, 26, 90, 60)):
    readme.column_dimensions[col].width = w
wb.save(a.out)

print(f"Saved {os.path.abspath(a.out)}\n\nFound:")
for k, l, p in found:
    print(f"  [{k}] {l}\n        {p}")
print("\nMissing (give with --file KEY=PATH, then rerun):" if missing else "\nNothing missing.")
for k, l, pats in missing:
    print(f"  [{k}] {l}   looked for: {', '.join(pats)}")
