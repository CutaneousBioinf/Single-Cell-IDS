"""
program_activity.py — where each program is active, for any dataset.

Scores every program in every cell, one sample at a time, and writes the summary by cell type,
disease value and donor. The per-cell scores are never kept.

  score of a cell = mean log-expression of the program's genes minus the mean of 50 control genes
                    with matched expression (scanpy's score_genes; controls are chosen once from the
                    dataset's own gene means)
  "expressing"    = score above --cutoff (default 0)

Programs come from a discovery run folder (the programs_<stratum>.csv files written by
ids_programs_disease.py), or from any CSV with a `genes` column and the program's name in
`program_id` or in the first column (--programs), which is how an existing program list, such as
the catalogue's 57 distinct programs, is scored in a new dataset.

  python scripts/02_ids_programs/program_activity.py --run results/02_ids_programs/ids-NAME
  python scripts/02_ids_programs/program_activity.py --programs results/03_catalogue/distinct_programs.csv \
      --datasets NAME --out results/02_ids_programs/catalogue-in-NAME

Writes into the run folder (or --out):
  detections.csv             one row per program found, with its genes. `distinct_program` gives the same
                             label to detections that share at least half of their genes (Jaccard >= 0.5),
                             i.e. the same program found in more than one condition. Run-folder mode only.
  activity_by_donor.csv.gz   dataset, program_id, disease, cell_type, donor: n_cells, mean_score, prop_expressing
  activity_by_celltype.csv   the same pooled over donors, plus n_donors

Checkpointed per sample; with --budget SECONDS it stops early and the same command continues it.
"""
import argparse, glob, hashlib, json, os, sys, time
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import units                                                     # noqa: E402

JACCARD = 0.5


def detections(run):
    """All programs of a run folder, with duplicates across strata given one `distinct_program` label."""
    t = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(run, "programs_*.csv")))
                   if not os.path.basename(f).startswith("programs_dropped")], ignore_index=True)
    if t.empty:
        raise SystemExit(f"{run}: no programs found. Run ids_programs_disease.py with --min-ids first.")
    t.insert(0, "program_id", t.stratum + "|" + t.program)
    G = [set(x.strip() for x in g.split(",")) for g in t.genes]
    parent = list(range(len(t)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i in range(len(t)):
        for j in range(i + 1, len(t)):
            if len(G[i] & G[j]) / len(G[i] | G[j]) >= JACCARD:
                parent[root(i)] = root(j)
    grp = pd.Series([root(i) for i in range(len(t))])
    order = {g: k for k, g in enumerate(pd.unique(grp), 1)}
    t.insert(1, "distinct_program", ["D%02d" % order[g] for g in grp])
    t.insert(2, "found_in", t.groupby("distinct_program").stratum.transform(lambda s: ", ".join(sorted(set(s)))))
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None, help="discovery run folder (programs are read from it, tables written to it)")
    ap.add_argument("--programs", default=None, help="CSV with program_id, genes: score these instead of the run's programs")
    ap.add_argument("--datasets", default=None, help="comma-separated; default: the datasets of the run")
    ap.add_argument("--out", default=None, help="default: the run folder")
    ap.add_argument("--results-dir", default=None,
                    help="read the datasets from <results-dir>/00_units. Not needed with --run: the run folder "
                         "remembers where its datasets are")
    ap.add_argument("--cutoff", type=float, default=0.0, help="a cell expresses a program when its score is above this")
    ap.add_argument("--lanes", default="october", choices=["october", "all"], help="blood only: nine October lanes or all 28")
    ap.add_argument("--budget", type=float, default=None)
    a = ap.parse_args()
    t0 = time.time()
    out = a.out or a.run
    if not out or not (a.run or a.programs):
        ap.error("give --run, or --programs with --datasets and --out")
    os.makedirs(out, exist_ok=True)
    if a.programs:
        det = pd.read_csv(a.programs)
        if "program_id" not in det.columns:                   # e.g. distinct_programs.csv, whose first column is program_no
            det = det.rename(columns={det.columns[0]: "program_id"})
        det["program_id"] = det.program_id.astype(str)
    else:
        det = detections(a.run)
        det.to_csv(os.path.join(out, "detections.csv"), index=False)
        print(f"{len(det)} detections, {det.distinct_program.nunique()} distinct programs -> detections.csv", flush=True)
    progs = {r.program_id: [g.strip() for g in r.genes.split(",")] for r in det.itertuples()}
    cfg = os.path.join(a.run or "", "run_config.json")
    cfg = json.load(open(cfg)) if os.path.exists(cfg) else {}
    if a.results_dir:
        units.use_results_dir(a.results_dir)
    elif cfg.get("units_dir") and os.path.isdir(cfg["units_dir"]):
        units.WORK = cfg["units_dir"]
    names = a.datasets.split(",") if a.datasets else cfg.get("datasets")
    if not names:
        ap.error("--datasets is needed (the run folder has no run_config.json)")

    parts = []
    for ds in names:
        us = units.stage_units(ds, "all" if a.lanes == "all" else None)
        key = hashlib.md5(json.dumps([ds, us, sorted(progs.items()), a.cutoff]).encode()).hexdigest()[:10]
        tmp = os.path.join(units.WORK, "scores", f"activity_{ds}_{key}"); os.makedirs(tmp, exist_ok=True)
        mean, n = units.gene_means(ds, us, cache=os.path.join(tmp, "_gene_means.npz"))
        ctrl = units.control_sets(units.genes(ds), mean, progs)
        print(f"[{ds}] {n:,} cells in {len(us)} samples; {len(ctrl)} of {len(progs)} programs have at least 2 genes in the data "
              f"({time.time()-t0:.0f}s)", flush=True)
        for u in us:
            f = os.path.join(tmp, f"{u}.pkl")
            if not os.path.exists(f):
                if a.budget is not None and time.time() - t0 > a.budget:
                    print("[budget] stopping; run the same command again to continue"); return
                ad_ = units.read_unit(ds, u)
                S, ids = units.score_matrix(units.lognorm(ad_.X), ctrl)
                o = ad_.obs
                code, grp = pd.factorize(pd.MultiIndex.from_arrays(
                    [o["disease_fine"].astype(str), o["subtype"].astype(str), o["donor"].astype(str)]))
                M = sp.csr_matrix((np.ones(len(code)), (code, np.arange(len(code)))), shape=(len(grp), len(code)))
                idx = pd.MultiIndex.from_tuples(list(grp), names=["disease", "cell_type", "donor"])
                long = pd.concat({"sum": pd.DataFrame(M @ S, index=idx, columns=ids).stack(),
                                  "pos": pd.DataFrame(M @ (S > a.cutoff).astype(float), index=idx, columns=ids).stack()}, axis=1)
                long["n_cells"] = np.repeat(np.asarray(M.sum(axis=1)).ravel(), len(ids))
                long.index.names = ["disease", "cell_type", "donor", "program_id"]
                long.to_pickle(f + ".tmp"); os.replace(f + ".tmp", f)
            parts.append(pd.read_pickle(f).reset_index().assign(dataset=ds))

    L = pd.concat(parts, ignore_index=True)
    by = ["dataset", "program_id", "disease", "cell_type"]
    D = L.groupby(by + ["donor"], as_index=False, sort=True)[["n_cells", "sum", "pos"]].sum()      # a donor can be in several lanes
    C = D.groupby(by, as_index=False, sort=True).agg(n_cells=("n_cells", "sum"), n_donors=("donor", "nunique"),
                                                     sum=("sum", "sum"), pos=("pos", "sum"))
    for T in (D, C):
        T["mean_score"] = T.pop("sum") / T.n_cells
        T["prop_expressing"] = T.pop("pos") / T.n_cells
        T["n_cells"] = T.n_cells.astype(int)
    D.to_csv(os.path.join(out, "activity_by_donor.csv.gz"), index=False, compression="gzip")
    C.to_csv(os.path.join(out, "activity_by_celltype.csv"), index=False)
    print(f"done ({time.time()-t0:.0f}s): {out}\n  activity_by_celltype.csv  {len(C):,} rows\n  activity_by_donor.csv.gz  {len(D):,} rows")
    if "distinct_program" in det.columns:
        top = (C[C.n_cells >= 50].sort_values("mean_score", ascending=False).drop_duplicates("program_id")
               .merge(det[["program_id", "distinct_program", "n_genes", "genes"]], on="program_id")
               .sort_values("mean_score", ascending=False))
        top["genes"] = top.genes.str.slice(0, 60)
        print("\nwhere each detection scores highest (groups of at least 50 cells):")
        print(top[["program_id", "distinct_program", "disease", "cell_type", "n_cells", "mean_score", "prop_expressing", "genes"]]
              .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
