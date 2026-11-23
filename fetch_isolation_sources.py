#!/usr/bin/env python3
"""Retrieve isolation source and related sample metadata from NCBI BioSample for every genome.

For each genome file in --genomes, the first sequence accession in the FASTA headers
(e.g. CP102930.1, NZ_..., JAxxxx01000001.1) is linked nuccore -> BioSample with NCBI E-utilities.
Genomes without an accession in their headers (or without a nuccore -> BioSample link) are
searched in BioSample by strain name; a single hit is used and flagged, several are listed for
you to choose. Nothing in the pipeline metadata is changed: review the table first.
"""
import argparse, glob, json, os, re, sys, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genomes", required=True, help="folder with the genome FASTA files")
ap.add_argument("--metadata", required=True, help="table with genome and strain columns")
ap.add_argument("--out", default="isolation_sources.tsv")
ap.add_argument("--email", default="", help="optional: NCBI asks for a contact email")
ap.add_argument("--api-key", default="", help="optional NCBI API key (faster)")
ap.add_argument("--debug-dir", default="", help="optional: save raw NCBI responses here")
a = ap.parse_args()

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DELAY = 0.12 if a.api_key else 0.4          # NCBI limit: 3 requests/s (10 with a key)
FIELDS = ["isolation_source", "env_broad_scale", "env_local_scale", "env_medium", "host",
          "sample_type", "geo_loc_name", "lat_lon", "collection_date"]
ACC = re.compile(r"\b((?:[A-Z]{2}_)?[A-Z]{1,6}\d{5,}(?:\.\d+)?)\b")


def get(util, **params):
    params.update(tool="bdello_pipeline")
    if a.email:
        params["email"] = a.email
    if a.api_key:
        params["api_key"] = a.api_key
    url = EU + util + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            time.sleep(DELAY)
            with urllib.request.urlopen(url, timeout=60) as r:
                txt = r.read().decode("utf-8", "replace")
            if a.debug_dir:
                os.makedirs(a.debug_dir, exist_ok=True)
                name = re.sub(r"[^A-Za-z0-9]+", "_", f"{util}_{params.get('id', params.get('term', ''))}")[:120]
                open(os.path.join(a.debug_dir, name + ".txt"), "w").write(txt)
            return txt
        except Exception as e:                      # network hiccup or HTTP 429: wait and retry
            err = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"NCBI request failed after retries: {err}")


def first_accession(path):
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            m = ACC.search(line[1:])
            if m:
                return m.group(1)
    return ""


def nuccore_to_biosample(acc):
    ids = json.loads(get("esearch.fcgi", db="nuccore", term=f"{acc}[accn]", retmode="json"))
    ids = ids.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return [], "accession not found in nuccore"
    link = json.loads(get("elink.fcgi", dbfrom="nuccore", db="biosample", id=ids[0], retmode="json"))
    out = []
    for ls in link.get("linksets", []):
        for db in ls.get("linksetdbs", []):
            if db.get("dbto") == "biosample":
                out += db.get("links", [])
    return out, "" if out else "no BioSample linked to this sequence"


def biosample_by_strain(strain):
    term = (f'"{strain}" AND (Bdellovibrio[Organism] OR Pseudobdellovibrio[Organism] '
            f'OR Bdellovibrionaceae[Organism] OR Bdellovibrionota[Organism])')
    r = json.loads(get("esearch.fcgi", db="biosample", term=term, retmode="json", retmax=20))
    return r.get("esearchresult", {}).get("idlist", [])


def biosample_record(uid):
    root = ET.fromstring(get("efetch.fcgi", db="biosample", id=uid, retmode="xml"))
    bs = root.find(".//BioSample")
    if bs is None:
        return {}
    rec = {"biosample": bs.get("accession", ""), "title": (bs.findtext("./Description/Title") or "").strip(),
           "organism": (bs.find("./Description/Organism").get("taxonomy_name", "")
                        if bs.find("./Description/Organism") is not None else "")}
    for at in bs.findall("./Attributes/Attribute"):
        key = (at.get("harmonized_name") or at.get("attribute_name") or "").strip().lower().replace(" ", "_")
        val = (at.text or "").strip()
        if key in FIELDS and val and key not in rec:
            rec[key] = val
        if key == "strain" and val:
            rec["biosample_strain"] = val
    return rec


meta = pd.read_csv(a.metadata, sep="\t", dtype=str).fillna("").set_index("genome")
files = {}
for ext in ("fna", "fa", "fasta", "fas", "fsa"):
    for f in glob.glob(os.path.join(a.genomes, f"*.{ext}")):
        files[os.path.basename(f)[: -len(ext) - 1]] = f
rows = []
for g in meta.index:
    strain = meta.loc[g, "strain"] if "strain" in meta.columns else ""
    row = dict(genome=g, strain=strain, accession="", biosample="", matched_by="", note="")
    f = files.get(g)
    if f is None:
        row["note"] = "genome file not found"
    else:
        row["accession"] = first_accession(f)
    uids, why = [], ""
    try:
        if row["accession"]:
            uids, why = nuccore_to_biosample(row["accession"])
            if uids:
                row["matched_by"] = "sequence accession"
        if not uids and strain:
            hits = biosample_by_strain(strain)
            if len(hits) == 1:
                uids, row["matched_by"] = hits, "strain name (check)"
            elif len(hits) > 1:
                accs = [biosample_record(h).get("biosample", h) for h in hits[:8]]
                why = (why + "; " if why else "") + f"{len(hits)} BioSamples match strain name: {', '.join(accs)}"
            else:
                why = (why + "; " if why else "") + "no BioSample found by strain name"
        if uids:
            rec = biosample_record(uids[0])
            row.update({k: v for k, v in rec.items() if k in FIELDS + ["biosample", "title"]})
            if rec.get("biosample_strain") and strain and \
                    re.sub(r"[^a-z0-9]", "", rec["biosample_strain"].lower()) != re.sub(r"[^a-z0-9]", "", strain.lower()):
                why = (why + "; " if why else "") + f"BioSample strain is '{rec['biosample_strain']}'"
            if len(uids) > 1:
                why = (why + "; " if why else "") + f"{len(uids)} BioSamples linked, first used"
    except Exception as e:
        why = (why + "; " if why else "") + f"error: {e}"
    row["note"] = "; ".join(x for x in (row["note"], why) if x)
    rows.append(row)
    print(f"  {g:45s} {row.get('biosample', '') or '-':14s} {row.get('isolation_source', '') or '-'}", flush=True)

cols = ["genome", "strain", "accession", "biosample", "matched_by"] + FIELDS + ["title", "note"]
df = pd.DataFrame(rows).reindex(columns=cols).astype(object).fillna("")
df.to_csv(a.out, sep="\t", index=False)
print(f"\n{(df.biosample != '').sum()} of {len(df)} genomes linked to a BioSample; "
      f"isolation_source given for {(df.isolation_source != '').sum()}. Table: {a.out}")
need = df[(df.biosample == "") | df.matched_by.str.contains("check")]
if len(need):
    print("Check these by hand:\n" + need[["genome", "accession", "matched_by", "note"]].to_string(index=False))
