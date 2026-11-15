"""
build_units.py — one file per blood lane, per skin sample, or per sample of any new dataset
(see sc_ids/units.py).

FIRST STEP OF THE PIPELINE: every later script reads these files, one at a time, so memory is
set by the largest single lane (about 1.5 GB) and not by how many lanes are used. The cells are
exactly those of the merged objects this replaced (74,126 blood cells in the nine October lanes,
67,367 skin cells; checked cell for cell).

Filters, applied per unit:
  blood  source study's IncludedInStudy + demuxlet singlets, Status in Healthy/PSO/PSA,
         >= 200 genes, <= 10% mitochondrial, Platelet and Eryth cells removed
  skin   cells present in the source metadata, >= 200 genes, <= 10% mitochondrial
Checkpointed: a unit that already has a file is skipped, so the script can be re-run until done.

  python scripts/02_ids_programs/build_units.py --tissue blood [--units PBMC-01-1,PBMC-03-1] [--budget 150]
  python scripts/02_ids_programs/build_units.py --tissue skin
  python scripts/02_ids_programs/build_units.py --tissue blood --index      # write cells_/genes_ tables

Any other dataset (one matrix.mtx + barcodes.tsv + features.tsv + a metadata .csv; no QC unless asked):
  python scripts/02_ids_programs/build_units.py --dataset NAME --matrix-dir DIR --metadata DIR/metadata.csv \
      --sample-col patientID --celltype-col simplified_subtype --disease-col disease
"""
import argparse, glob, gzip, os, re, sys, time, resource
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units                                              # noqa: E402

KEEP_STATUS = {"Healthy", "PSO", "PSA"}
MIN_GENES, MAX_MT = 200, 10.0
EXCLUDE = {"Platelet", "Eryth"}
ALL_LANES = [f"PBMC-0{r}-{l}" for r in range(1, 8) for l in range(1, 5)]


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def qc(X, var_names, obs):
    return units.qc(X, var_names, obs, MIN_GENES, MAX_MT)


save = units.write_unit


_META = {}
def blood_meta():
    if "b" not in _META:
        m = pd.read_csv(paths.PBMC_META, sep="\t", low_memory=False,
                        usecols=["CellName", "Sample", "Subject", "Status", "DemuxletDropletType",
                                 "IncludedInStudy", "CellType"])
        m = m[m.IncludedInStudy.astype(str).str.upper().isin(["TRUE", "T", "1"])
              & (m.DemuxletDropletType.astype(str) == "SNG") & m.Status.isin(KEEP_STATUS)]
        _META["b"] = m.drop_duplicates("CellName").set_index("CellName")
    return _META["b"]


def build_lane(lane):
    try:
        import fast_matrix_market as fmm          # seconds per lane
    except ImportError:
        import scipy.io as fmm                    # same result, about a minute per lane
    from anndata.utils import make_index_unique
    d = paths.PBMC_DIR
    with gzip.open(f"{d}/{lane}.matrix.mtx.gz", "rb") as f:
        X = sp.csr_matrix(fmm.mmread(f).T, dtype=np.float32)          # cells x genes
    bc = pd.read_csv(f"{d}/{lane}.barcodes.tsv.gz", header=None)[0].astype(str)
    ft = pd.read_csv(f"{d}/{lane}.features.tsv.gz", header=None, sep="\t")
    # CITE-seq lanes also carry antibody rows; keep RNA only (what scanpy's gex_only=True does).
    # Leaving them in inflates total counts, lowers the mitochondrial share and lets extra cells through QC.
    gex = (ft[2].astype(str) == "Gene Expression").values if ft.shape[1] > 2 else np.ones(len(ft), bool)
    X, ft = X[:, gex], ft[gex]
    var = make_index_unique(pd.Index(ft[1].astype(str).values))
    n_raw = X.shape[0]
    key = lane + "_" + bc.str.replace(r"-\d+$", "", regex=True)
    m = blood_meta()
    hit = key.isin(m.index).values
    X, key, bc = X[hit], key[hit], bc[hit]
    mm = m.loc[key.values]
    obs = pd.DataFrame({"donor": mm.Subject.astype(str).values, "subtype": mm.CellType.astype(str).values,
                        "disease_fine": mm.Status.astype(str).values, "unit": lane},
                       index=pd.Index(lane + "_" + bc.values, name="cell"))
    X, obs = qc(X, var, obs)
    ok = ~obs.subtype.isin(EXCLUDE).values
    X, obs = X[ok], obs[ok]
    save("blood", lane, X, var, obs)
    return n_raw, X.shape[0]


def skin_symbols():
    from anndata.utils import make_index_unique
    ft = pd.read_csv(f"{paths.PBMC_DIR}/PBMC-01-1.features.tsv.gz", header=None, sep="\t")
    return dict(zip(ft[0].astype(str), ft[1].astype(str))), make_index_unique


def build_skin_sample(f, meta):
    sample = re.sub(r"^GSM\d+_", "", os.path.basename(f)).replace(".csv.gz", "").replace("-", "_", 1)
    rows = meta[meta["orig.ident"] == sample]
    bare = pd.Index(rows.index).str.replace(r"-\d+$", "", regex=True)
    suffix = str(rows.index[0]).rsplit("-", 1)[-1]
    wanted = set(bare + "-1")
    with gzip.open(f, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
    index_name = header[0] or "Unnamed: 0"
    take = [c for c in header[1:] if c.strip('"') in wanted]
    assert len(take) == len(wanted), f"{sample}: {len(wanted) - len(take)} metadata barcodes absent"
    chunks, idx, colnames = [], [], None
    for ch in pd.read_csv(f, index_col=0, usecols=[index_name] + take,
                          dtype={c: np.int32 for c in take}, chunksize=4000):
        colnames = ch.columns if colnames is None else colnames
        idx.append(ch.index.astype(str).to_numpy())
        chunks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32)))
    ens = np.concatenate(idx)
    X = sp.vstack(chunks, format="csr").T.tocsr()                         # cells x genes
    keys = pd.Index(colnames).str.replace(r"-\d+$", "", regex=True) + "-" + suffix
    sym, uniq = skin_symbols()
    var = uniq(pd.Index([sym.get(e) or e for e in ens]))
    mm = meta.loc[keys]
    obs = pd.DataFrame({"donor": mm["donor"].astype(str).values, "subtype": mm["fullsubtype"].astype(str).values,
                        "disease_fine": mm["disease"].astype(str).values, "unit": sample},
                       index=pd.Index(keys, name="cell"))
    n_raw = X.shape[0]
    X, obs = qc(X, var, obs)
    ok = ~obs.subtype.isin(EXCLUDE).values
    save("skin", sample, X[ok], var, obs[ok])
    return sample, n_raw, int(ok.sum())


write_index = units.write_index


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", choices=["blood", "skin"], help="one of the two catalogue datasets")
    ap.add_argument("--units", default=None)
    ap.add_argument("--budget", type=float, default=150, help="stop starting new units after this many seconds")
    ap.add_argument("--index", action="store_true")
    g = ap.add_argument_group("any other dataset (docs/new_dataset_guide.md)")
    g.add_argument("--dataset", help="a name you choose for the new dataset")
    g.add_argument("--matrix-dir", help="folder with matrix.mtx, barcodes.tsv, features.tsv (plain or .gz)")
    g.add_argument("--metadata", help="metadata .csv, one row per cell, barcodes in the first column")
    g.add_argument("--sample-col", default="patientID")
    g.add_argument("--celltype-col", default="simplified_subtype")
    g.add_argument("--disease-col", default="disease")
    g.add_argument("--barcode-col", default=None, help="metadata column holding the barcodes (default: the first column)")
    g.add_argument("--min-genes", type=int, default=0, help="QC is off unless set")
    g.add_argument("--max-mt", type=float, default=100.0, help="QC is off unless set")
    a = ap.parse_args()
    t0 = time.time()
    if a.dataset:
        if not (a.matrix_dir and a.metadata):
            ap.error("--dataset needs --matrix-dir and --metadata")
        units.build_dataset(a.dataset, a.matrix_dir, a.metadata, a.sample_col, a.celltype_col, a.disease_col,
                            a.barcode_col, a.min_genes, a.max_mt)
        print(f"done in {time.time()-t0:.0f}s, peak memory {rss():.1f} GB. Files: {os.path.join(units.WORK, a.dataset)}")
        return
    if not a.tissue:
        ap.error("give --tissue blood|skin, or --dataset NAME with --matrix-dir and --metadata")
    if a.tissue == "blood":
        todo = a.units.split(",") if a.units else ALL_LANES
        if not a.index:
            for u in todo:
                if os.path.exists(units.unit_path("blood", u)):
                    continue
                if time.time() - t0 > a.budget:
                    print("[budget] stopping; re-run to continue"); break
                t = time.time(); n0, n1 = build_lane(u)
                print(f"[blood] {u}: {n0:,} droplets -> {n1:,} cells  {time.time()-t:.0f}s  peak {rss():.2f} GB", flush=True)
        else:
            write_index("blood", ALL_LANES)
    else:
        files = sorted(glob.glob(os.path.join(paths.SKIN_DIR, "*.csv.gz")))
        names = [re.sub(r"^GSM\d+_", "", os.path.basename(f)).replace(".csv.gz", "").replace("-", "_", 1) for f in files]
        if not a.index:
            meta = pd.read_csv(paths.SKIN_META, sep="\t", index_col=0, low_memory=False)
            meta.index = meta.index.astype(str)
            for f, nm in zip(files, names):
                if os.path.exists(units.unit_path("skin", nm)):
                    continue
                if time.time() - t0 > a.budget:
                    print("[budget] stopping; re-run to continue"); break
                t = time.time(); s, n0, n1 = build_skin_sample(f, meta)
                print(f"[skin] {s}: {n0:,} -> {n1:,} cells  {time.time()-t:.0f}s  peak {rss():.2f} GB", flush=True)
        else:
            write_index("skin", names)


if __name__ == "__main__":
    main()
