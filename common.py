"""Shared helpers for the Bdellovibrio pipeline."""
import os
import numpy as np


def read_fasta(path):
    """Return dict {id(first token): sequence} and dict {id: full header}."""
    seqs, heads, cur = {}, {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                heads[cur] = line[1:]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}, heads


def genome_ids(out_dir):
    with open(os.path.join(out_dir, "genome_ids.txt")) as fh:
        return [l.strip() for l in fh if l.strip()]


def bh_fdr(p):
    """Benjamini–Hochberg adjusted p-values (NaN-safe)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def parse_domtbl(path, max_evalue=1e-5, min_cov=0.30):
    """Parse hmmscan --domtblout. Keeps domains with i-Evalue <= max_evalue and
    HMM coverage ((hmm_to-hmm_from+1)/model_length) >= min_cov.
    Yields (protein_id, pfam_acc_unversioned, pfam_name, i_evalue, coverage)."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            name, acc, tlen, query = f[0], f[1].split(".")[0], int(f[2]), f[3]
            i_eval = float(f[12])
            hmm_from, hmm_to = int(f[15]), int(f[16])
            cov = (hmm_to - hmm_from + 1) / tlen
            if i_eval <= max_evalue and cov >= min_cov:
                yield query, acc, name, i_eval, cov


def load_metadata(out_dir, ids, exclude_outgroups=True):
    """Return (meta indexed by genome, list of genomes used for statistics).
    Genomes with outgroup == 'yes' are left out of statistics when exclude_outgroups."""
    import pandas as pd
    f = os.path.join(out_dir, "metadata.tsv")
    if not os.path.exists(f):
        return None, list(ids)
    meta = pd.read_csv(f, sep="\t", dtype=str).set_index("genome")
    missing = [g for g in ids if g not in meta.index]
    if missing:
        raise SystemExit(f"metadata.tsv has no row for: {', '.join(missing)}")
    meta = meta.reindex(ids)
    def flag(col):
        return meta[col].fillna("").astype(str).str.lower().eq("yes") if col in meta.columns \
            else pd.Series(False, index=meta.index)
    drop = flag("exclude_from_stats")              # e.g. incomplete genomes (always dropped)
    if exclude_outgroups:
        drop = drop | flag("outgroup")
    stat = [g for g in ids if not drop[g]]
    return meta, stat


def single_case_test(case_values, ref_df):
    """Compare one genome with a reference group, column by column.
    Crawford & Howell (1998) modified t-test: t = (x - mean) / (sd * sqrt((n + 1) / n)), df = n - 1.
    Returns a DataFrame with the case value, reference mean/sd/range, z, percentile, t, two-sided p,
    Bonferroni and Benjamini–Hochberg adjusted p."""
    import pandas as pd
    from scipy.stats import t as tdist
    rows = []
    for c in ref_df.columns:
        r = ref_df[c].astype(float).dropna(); x = float(case_values[c]); n = len(r)
        m, sd = r.mean(), r.std(ddof=1)
        if n < 3 or sd == 0 or np.isnan(sd):
            tval, p = np.nan, np.nan
        else:
            tval = (x - m) / (sd * np.sqrt((n + 1) / n))
            p = 2 * tdist.sf(abs(tval), n - 1)
        rows.append(dict(feature=c, foreground=x, ref_n=n, ref_mean=m, ref_sd=sd, ref_min=r.min(),
                         ref_max=r.max(), z=(x - m) / sd if sd else np.nan,
                         percentile=100 * ((r < x).sum() + 0.5 * (r == x).sum()) / n,
                         t_crawford_howell=tval, p=p))
    df = pd.DataFrame(rows)
    df["p_bonferroni"] = np.minimum(df.p * df.p.notna().sum(), 1.0)
    df["q_BH"] = bh_fdr(df.p)
    df["direction"] = np.where(df.foreground > df.ref_mean, "higher", "lower")
    return df


def tip_matcher(strain_of):
    """Return match(tip) -> genome ID, None, or 'AMBIGUOUS:...'. strain_of = {genome: strain}.
    A tip matches a genome if it equals the genome ID or contains the strain name as a whole
    token (case-insensitive; '-', '_', '.' or space allowed between characters); longest wins."""
    import re
    norm = lambda x: re.sub(r"[^a-z0-9]", "", x.lower())
    pats = {g: re.compile(r"(?<![A-Za-z0-9])" + r"[-_. ]?".join(map(re.escape, norm(st)))
                          + r"(?![A-Za-z0-9])", re.I) for g, st in strain_of.items() if norm(st)}
    def match(tip):
        if tip in strain_of:
            return tip
        hits = sorted(((len(norm(strain_of[g])), g) for g, p in pats.items() if p.search(tip)), reverse=True)
        if not hits:
            return None
        if len(hits) > 1 and hits[0][0] == hits[1][0]:
            return "AMBIGUOUS:" + ",".join(h[1] for h in hits if h[0] == hits[0][0])
        return hits[0][1]
    return match


def single_case_regression(x_case, y_case, x_ref, y_ref):
    """Size-adjusted single-case test: is the case's y higher/lower than expected from x?
    Ordinary least squares y ~ x fitted on the reference genomes only; the case's residual is tested
    against the prediction interval (Crawford & Garthwaite 2006):
        t = (y - y_hat) / (s_yx * sqrt(1 + 1/n + (x - x_mean)^2 / ((n - 1) * s_x^2))),  df = n - 2.
    Returns dict(expected, observed, residual, t, p, slope, r2)."""
    from scipy.stats import t as tdist
    x_ref = np.asarray(x_ref, float); y_ref = np.asarray(y_ref, float); n = len(x_ref)
    if n < 4 or np.var(x_ref) == 0:
        return dict(expected=np.nan, observed=y_case, residual=np.nan, t=np.nan, p=np.nan, slope=np.nan, r2=np.nan)
    b, a0 = np.polyfit(x_ref, y_ref, 1)
    fit = a0 + b * x_ref
    ss_res = float(((y_ref - fit) ** 2).sum()); ss_tot = float(((y_ref - y_ref.mean()) ** 2).sum())
    s_yx = np.sqrt(ss_res / (n - 2))
    y_hat = a0 + b * x_case
    if s_yx == 0:
        return dict(expected=y_hat, observed=y_case, residual=y_case - y_hat, t=np.nan, p=np.nan, slope=b,
                    r2=1.0 if ss_tot else np.nan)
    se = s_yx * np.sqrt(1 + 1 / n + (x_case - x_ref.mean()) ** 2 / ((n - 1) * x_ref.var(ddof=1)))
    tval = (y_case - y_hat) / se
    return dict(expected=y_hat, observed=y_case, residual=y_case - y_hat, t=tval,
                p=2 * tdist.sf(abs(tval), n - 2), slope=b, r2=1 - ss_res / ss_tot if ss_tot else np.nan)
