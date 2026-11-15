"""
units.py — the data as one file per UNIT (a blood lane, a skin sample, or a sample of any other
dataset), never merged.

Why: the earlier scripts loaded a whole tissue into memory, so memory grew with every lane
added and the 9-lane build was killed at 6 GB. Here every step reads one unit at a time, so
peak memory is set by the largest single lane, not by how many lanes are used.

Layout under WORK (env SC_IDS_WORK, default results/00_units):
  <tissue>/<unit>.h5ad      raw counts (CSR float32) + obs: donor, subtype, disease_fine
  cells_<tissue>.csv.gz     one row per cell: cell, unit, row, donor, subtype, disease_fine, total_counts
  genes_<tissue>.txt        gene order shared by every unit of the tissue
"<tissue>" is simply the dataset's name. `blood` and `skin` are the two GEO datasets of the catalogue;
any other name is a dataset added with build_dataset() below (docs/new_dataset_guide.md).

Scoring is two passes over the units and reproduces scanpy's `sc.tl.score_genes` on the
merged object: pass 1 gets each gene's mean log-expression over ALL cells, the control genes
are then picked once by scanpy's own binning code from those means, and pass 2 scores each
unit against the same controls.
"""
import gzip, json, os, re
import numpy as np
import pandas as pd
import scipy.sparse as sp
from . import paths

WORK = paths.UNITS
TISSUES = ("blood", "skin")
# disease label shared by the two tissues; the source label is kept as disease_fine
DISEASE = {"Healthy": "healthy", "PSO": "psoriasis", "PSA": "PsA", "NS": "healthy", "PN": "psoriasis", "PP": "psoriasis"}


def unit_path(tissue, unit):
    return os.path.join(WORK, tissue, f"{unit}.h5ad")


def cells(tissue):
    """Cell table of a dataset. Labels are read as text, so a donor called 007 or a disease coded 1 stays as written."""
    f = os.path.join(WORK, f"cells_{tissue}.csv.gz")
    if not os.path.exists(f):
        raise FileNotFoundError(f"no dataset called '{tissue}' under {WORK} (expected {f}). "
                                "Build it first with scripts/02_ids_programs/build_units.py")
    text = {k: str for k in ("cell", "donor", "subtype", "disease_fine", "unit")}
    return pd.read_csv(f, dtype=text, keep_default_na=False)


def genes(tissue):
    return [g.rstrip("\n") for g in open(os.path.join(WORK, f"genes_{tissue}.txt"))]


def read_unit(tissue, unit):
    import anndata as ad
    return ad.read_h5ad(unit_path(tissue, unit))


def lognorm(X, target=1e4):
    """normalize_total(target) over ALL genes, then log1p. CSR float32 in, CSR float32 out."""
    X = sp.csr_matrix(X, dtype=np.float32, copy=True)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    X.data *= np.repeat((target / tot).astype(np.float32), np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


def read_cells(tissue, sel, gene_list=None, norm=None):
    """Counts for selected cells. `sel` is a slice of cells(tissue) (needs columns unit,row);
    rows come back in the order of `sel`. norm='total' divides by each cell's total over ALL
    genes (x1e4) before the gene subset is taken; None returns raw counts."""
    g_all = genes(tissue)
    cols = None if gene_list is None else pd.Index(g_all).get_indexer(gene_list)
    parts, order = [], []
    for unit, s in sel.groupby("unit", sort=False):
        a = read_unit(tissue, unit)
        X = sp.csr_matrix(a.X, dtype=np.float32)[s["row"].values]
        if norm == "total":
            tot = np.asarray(X.sum(axis=1)).ravel(); tot[tot == 0] = 1.0
            X = sp.diags((1e4 / tot).astype(np.float32)) @ X
        parts.append(sp.csr_matrix(X[:, cols]) if cols is not None else X)
        order.append(s.index.values)
        del a
    X = sp.vstack(parts, format="csr")
    pos = pd.Series(np.arange(len(sel)), index=np.concatenate(order)).loc[sel.index.values].values
    return X[pos]


# ---------------------------------------------------------------- scoring
def gene_means(tissue, units=None, cache=None):
    """Pass 1: mean log-normalized expression of every gene over all cells of the units."""
    if cache and os.path.exists(cache):
        z = np.load(cache); return z["mean"], int(z["n"])
    c = cells(tissue)
    units = list(units) if units is not None else list(c["unit"].unique())
    s, n = None, 0
    for u in units:
        X = lognorm(read_unit(tissue, u).X)
        v = np.asarray(X.sum(axis=0, dtype=np.float64)).ravel()
        s = v if s is None else s + v
        n += X.shape[0]
    m = s / n
    if cache:
        np.savez(cache, mean=m, n=n)
    return m, n


def control_sets(var_names, mean, programs, ctrl_size=50, n_bins=25, random_state=0):
    """Control genes per program, chosen by scanpy's own binning code from global means.
    `programs` is {id: [genes]}. Returns {id: (gene_idx, ctrl_idx)}; programs with <2 genes
    present are skipped, as in the earlier scripts."""
    from scanpy.tools._score_genes import _score_genes_bins
    vn = pd.Index(var_names)
    pool = vn.astype("string")
    row = np.asarray(mean, dtype=np.float64)[None, :]
    out = {}
    for pid, gl in programs.items():
        gl = pd.Index([g for g in gl if g in vn])
        if len(gl) < 2:
            continue
        gl = pd.Index(gl).intersection(vn)
        np.random.seed(random_state)
        ctrl = pd.Index([], dtype="string")
        for r in _score_genes_bins(gl, pool, ctrl_as_ref=True, ctrl_size=ctrl_size,
                                   n_bins=n_bins, get_subset=lambda genes: row):
            ctrl = ctrl.union(r)
        out[pid] = (vn.get_indexer(gl), vn.get_indexer(pd.Index(ctrl.astype(str))))
    return out


def score_matrix(X_lognorm, ctrl):
    """Pass 2 for one unit: cells x programs score = mean(program genes) - mean(controls)."""
    G = X_lognorm.shape[1]
    ids = list(ctrl)
    r, c, v = [], [], []
    for j, pid in enumerate(ids):
        gi, ci = ctrl[pid]
        r += list(gi) + list(ci); c += [j] * (len(gi) + len(ci))
        v += [1.0 / len(gi)] * len(gi) + [-1.0 / len(ci)] * len(ci)
    W = sp.csc_matrix((v, (r, c)), shape=(G, len(ids)), dtype=np.float64)
    return np.asarray((X_lognorm.astype(np.float64) @ W).todense()), ids


def stage_units(tissue, lanes=None):
    """Units a catalogue stage works on: the nine October blood lanes unless lanes='all'; every skin sample."""
    from .datasets.gse194315 import LANES_USED
    if tissue == "blood" and lanes != "all":
        return list(lanes or LANES_USED)
    return list(cells(tissue)["unit"].unique())


def iter_scores(tissue, programs, unit_list=None, cache=None):
    """Score `programs` ({id: [genes]}) in every cell, one unit at a time.
    Yields (unit, obs, S, ids): S is cells x programs for that unit, equal to what sc.tl.score_genes
    (ctrl_size 50, random_state 0) gives on the merged object of all `unit_list` units."""
    unit_list = list(unit_list) if unit_list is not None else stage_units(tissue)
    mean, _ = gene_means(tissue, unit_list, cache=cache)
    ctrl = control_sets(genes(tissue), mean, programs)
    for u in unit_list:
        a = read_unit(tissue, u)
        S, ids = score_matrix(lognorm(a.X), ctrl)
        yield u, a.obs, S, ids


# ---------------------------------------------------------------- writing units
def qc(X, var_names, obs, min_genes=0, max_mt=100.0):
    """Add total_counts, n_genes, pct_mt to obs and keep cells with >= min_genes genes and <= max_mt % mitochondrial."""
    mt = np.array([g.upper().startswith("MT-") for g in var_names])
    tot = np.asarray(X.sum(axis=1)).ravel()
    ng = np.asarray((X > 0).sum(axis=1)).ravel()
    pct = 100 * np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(tot, 1)
    keep = (ng >= min_genes) & (pct <= max_mt)
    obs = obs.assign(total_counts=tot, n_genes=ng, pct_mt=pct)
    return X[keep], obs[keep]


def write_unit(tissue, unit, X, var_names, obs):
    import anndata as ad
    p = unit_path(tissue, unit)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    a = ad.AnnData(X=sp.csr_matrix(X, dtype=np.float32), obs=obs,
                   var=pd.DataFrame(index=pd.Index(np.asarray(var_names, dtype=str), name="gene")))
    a.write_h5ad(p + ".tmp", compression="lzf"); os.replace(p + ".tmp", p)


def write_index(tissue, order):
    """cells_<tissue>.csv.gz and genes_<tissue>.txt from the unit files that exist, in the given unit order."""
    import anndata as ad
    rows, g0 = [], None
    for u in order:
        p = unit_path(tissue, u)
        if not os.path.exists(p):
            continue
        a = ad.read_h5ad(p, backed="r")
        o = a.obs[["donor", "subtype", "disease_fine", "unit", "total_counts"]].copy()
        o.insert(0, "row", np.arange(len(o))); o.insert(0, "cell", a.obs_names.values)
        rows.append(o)
        g = list(a.var_names)
        assert g0 is None or g == g0, f"{u}: gene order differs"
        g0 = g
        a.file.close()
    c = pd.concat(rows, ignore_index=True)
    c.to_csv(os.path.join(WORK, f"cells_{tissue}.csv.gz"), index=False)
    open(os.path.join(WORK, f"genes_{tissue}.txt"), "w").write("\n".join(g0) + "\n")
    print(f"[index] {tissue}: {len(c):,} cells in {c.unit.nunique()} units, {len(g0):,} genes")
    return c


def _find(folder, *names):
    for n in names:
        for ext in ("", ".gz"):
            if os.path.exists(os.path.join(folder, n + ext)):
                return os.path.join(folder, n + ext)
    raise FileNotFoundError(f"{folder}: none of {', '.join(names)} (plain or .gz)")


def build_dataset(name, matrix_dir, metadata, sample_col="patientID", celltype_col="simplified_subtype",
                  disease_col="disease", barcode_col=None, min_genes=0, max_mt=100.0):
    """A new dataset from ONE count matrix plus a metadata table, written as one file per sample.

    matrix_dir   folder with matrix.mtx, barcodes.tsv, features.tsv (or genes.tsv); plain or .gz
    metadata     .csv (or .tsv/.txt, tab-separated), one row per cell. Cells are matched on the exact
                 barcode text, e.g. AAACCTGAGAAACCAT-3; the barcodes are expected in `barcode_col`
                 (default: the first column, which is where R's write.csv puts the row names)
    sample_col, celltype_col, disease_col   the three metadata columns the pipeline needs
    min_genes, max_mt   QC is OFF by default (0 and 100): the data are taken as already filtered

    Cells with no metadata row, or with one of the three labels missing, are left out and counted.
    The whole matrix is read once, so this step needs memory for the full matrix; every later
    step reads one sample at a time.
    """
    from anndata.utils import make_index_unique
    if name in TISSUES:
        raise ValueError(f"'{name}' is reserved for the catalogue data; choose another dataset name")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise ValueError("dataset name: letters, digits, '.', '_' and '-' only (it becomes part of file names)")
    # --- small files first, so a wrong column name fails before the slow matrix read
    ft = pd.read_csv(_find(matrix_dir, "features.tsv", "genes.tsv"), header=None, sep="\t", dtype=str)
    bc = pd.read_csv(_find(matrix_dir, "barcodes.tsv"), header=None, sep="\t", dtype=str)[0]
    if bc.duplicated().any():
        raise ValueError(f"barcodes.tsv has {int(bc.duplicated().sum()):,} repeated barcodes, e.g. {bc[bc.duplicated()].iloc[0]}")
    m = pd.read_csv(metadata, sep="\t" if re.search(r"\.(tsv|txt)(\.gz)?$", metadata) else ",", dtype=str)
    bcol = barcode_col or m.columns[0]
    absent = [c for c in (bcol, sample_col, celltype_col, disease_col) if c not in m.columns]
    if absent:
        raise ValueError(f"metadata has no column {absent}. Its columns are: {list(m.columns)}")
    if m[bcol].duplicated().any():
        raise ValueError(f"metadata column '{bcol}' has repeated barcodes, e.g. {m[bcol][m[bcol].duplicated()].iloc[0]}")
    m = m.set_index(bcol)[[sample_col, celltype_col, disease_col]]
    hit = bc.isin(m.index).values
    if not hit.any():
        raise ValueError("no barcode in barcodes.tsv matches the metadata.\n"
                         f"  barcodes.tsv: {list(bc[:3])}\n  metadata '{bcol}': {list(m.index[:3])}\n"
                         "They must be the same text. If the barcodes are in another metadata column, pass --barcode-col.")
    lab = m.reindex(bc.values)
    keep = hit & lab.notna().all(axis=1).values
    print(f"[{name}] {len(bc):,} cells in the matrix, {len(m):,} rows in the metadata: {int(hit.sum()):,} matched, "
          f"{int((hit & ~keep).sum()):,} of those dropped for a missing label", flush=True)
    lab = lab[keep]
    per_sample = lab.groupby(sample_col, sort=False)[disease_col].nunique()
    if (per_sample > 1).any():
        print(f"  note: {int((per_sample > 1).sum())} samples carry more than one '{disease_col}' value, "
              f"e.g. {list(per_sample.index[per_sample > 1][:3])}", flush=True)
    samples = list(pd.unique(lab[sample_col]))
    safe = {s: re.sub(r"[^A-Za-z0-9._+-]+", "_", s) for s in samples}
    if len(set(safe.values())) < len(samples):
        raise ValueError("two sample names become the same file name after replacing special characters")

    # --- the matrix
    try:
        import fast_matrix_market as fmm          # seconds
    except ImportError:
        import scipy.io as fmm                    # same result, slower
    f = _find(matrix_dir, "matrix.mtx")
    with (gzip.open(f, "rb") if f.endswith(".gz") else open(f, "rb")) as fh:
        M = fmm.mmread(fh)
    if M.shape == (len(ft), len(bc)):
        X = sp.csr_matrix(M.T, dtype=np.float32)                      # 10x layout: genes x cells
    elif M.shape == (len(bc), len(ft)):
        X = sp.csr_matrix(M, dtype=np.float32)
    else:
        raise ValueError(f"matrix is {M.shape[0]:,} x {M.shape[1]:,} but there are {len(ft):,} features and {len(bc):,} barcodes")
    del M
    if ft.shape[1] > 2 and (ft[2] == "Gene Expression").any():       # drop antibody / other non-RNA rows
        gex = (ft[2] == "Gene Expression").values
        X, ft = X[:, gex], ft[gex]
    var = make_index_unique(pd.Index(ft[1 if ft.shape[1] > 1 else 0].values))
    X = X[keep]
    whole = np.allclose(X.data[:200000], np.round(X.data[:200000]))
    if not whole:
        print("  WARNING: the matrix does not hold whole-number counts. IDS and the scores expect raw counts.", flush=True)
    obs = pd.DataFrame({"donor": lab[sample_col].values, "subtype": lab[celltype_col].values,
                        "disease_fine": lab[disease_col].values,
                        "unit": lab[sample_col].map(safe).values}, index=pd.Index(lab.index.values, name="cell"))
    n_in = len(obs)
    X, obs = qc(X, var, obs, min_genes, max_mt)
    if len(obs) < n_in:
        print(f"  QC (>= {min_genes} genes, <= {max_mt}% mitochondrial) removed {n_in - len(obs):,} cells", flush=True)
    for s in samples:
        rows = np.where(obs.unit.values == safe[s])[0]
        if len(rows):
            write_unit(name, safe[s], X[rows], var, obs.iloc[rows])
    c = write_index(name, [safe[s] for s in samples])
    tab = c.groupby("disease_fine").agg(cells=("cell", "size"), samples=("donor", "nunique"))
    print(tab.to_string(), flush=True)
    print(f"  {c.subtype.nunique()} cell types: {', '.join(sorted(c.subtype.unique())[:12])}"
          + (" ..." if c.subtype.nunique() > 12 else ""), flush=True)
    json.dump({"name": name, "matrix_dir": os.path.abspath(matrix_dir), "metadata": os.path.abspath(metadata),
               "columns": {"barcode": bcol, "sample": sample_col, "cell_type": celltype_col, "disease": disease_col},
               "qc": {"min_genes": min_genes, "max_mt": max_mt}, "whole_number_counts": bool(whole),
               "cells_in_matrix": int(len(bc)), "cells_kept": int(len(c)), "genes": int(len(var)),
               "by_disease": {k: {"cells": int(v.cells), "samples": int(v.samples)} for k, v in tab.iterrows()}},
              open(os.path.join(WORK, f"dataset_{name}.json"), "w"), indent=1)
    return c
