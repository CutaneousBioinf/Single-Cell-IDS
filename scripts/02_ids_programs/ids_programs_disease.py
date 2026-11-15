"""
ids_programs_disease.py — IDS gene programs per (tissue x disease) stratum.

Six strata: PBMC x {Healthy, PSO, PSA} and skin x {NS, PN, PP}. Each gets its own gene-gene IDS
matrix, so a program can later be called tissue-shared, tissue-specific, disease-shared or
disease-specific by comparing strata. The work itself is in src/sc_ids/discovery.py; cells are
read from the per-lane files written by build_units.py, a few thousand at a time.

TWO DESIGN CHOICES THAT MATTER, both about comparability:
  1. ONE FIXED GENE PANEL for every stratum (union of per-stratum highly variable genes plus a
     curated marker list). Per-stratum gene selection would make a program look tissue-specific
     merely because its genes were never offered to the other tissue.
  2. A COMMON CELL CAP per stratum, drawn stratified by cell subtype so rare types survive.
     Note the cap does not equalise small strata: inside one lineage some skin strata have only
     400-900 cells, and small samples push IDS upward.

Pass 1 computes and caches the matrices and prints each stratum's off-diagonal IDS percentiles.
The edge threshold is read off those percentiles, then pass 2 extracts programs with --reuse.
Thresholds used for the catalogue: whole tissue 0.50, T+NK 0.35, Myeloid 0.40, Keratinocyte 0.50.

  python scripts/02_ids_programs/ids_programs_disease.py                      # pass 1, whole tissue
  python scripts/02_ids_programs/ids_programs_disease.py --reuse --min-ids 0.5
  python scripts/02_ids_programs/ids_programs_disease.py --lineage T+NK       # then --reuse --min-ids 0.35
Checkpointed (panel and each stratum's matrix are cached); with --budget SECONDS the script stops early
and --reuse continues it.

ANY OTHER DATASET (built with build_units.py --dataset NAME; see docs/new_dataset_guide.md): pass
--datasets NAME. The strata are then one per disease value found in that dataset's metadata, and the
results go to results/02_ids_programs/ids-NAME/. --celltypes "A,B" restricts every stratum to those
cell-type labels (the dataset's own labels), the equivalent of --lineage.

  python scripts/02_ids_programs/ids_programs_disease.py --datasets NAME                 # pass 1
  python scripts/02_ids_programs/ids_programs_disease.py --datasets NAME --reuse --min-ids 0.5
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from sc_ids import paths, units, programs as P, discovery as D               # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="default: the run folder registered in sc_ids/paths.py, "
                                                "or results/02_ids_programs/ids-<datasets> with --datasets")
    ap.add_argument("--datasets", default=None,
                    help="comma-separated dataset names built with build_units.py --dataset; one stratum per "
                         "disease value in each. Default: the six catalogue strata (blood and skin)")
    ap.add_argument("--celltypes", default=None,
                    help="comma-separated cell-type labels, exactly as in the data: restrict every stratum to them")
    ap.add_argument("--no-curated", action="store_true",
                    help="gene panel = highly variable genes only, without the marker list in sc_ids/discovery.py")
    ap.add_argument("--lineage", default=None,
                    help="restrict every stratum to one coarse lineage (T+NK, Myeloid, Keratinocyte). Panel, "
                         "detection floor and IDS are then all computed WITHIN the lineage, so a T-cell "
                         "cytokine is no longer diluted below the floor by keratinocytes.")
    ap.add_argument("--n-hvg", type=int, default=1200)
    ap.add_argument("--cell-cap", type=int, default=4000)
    ap.add_argument("--block-size", type=int, default=1024)
    ap.add_argument("--min-cells", type=int, default=400, help="skip a stratum with fewer cells than this")
    ap.add_argument("--min-cells-frac", type=float, default=0.02, help="detection floor, applied within each stratum")
    ap.add_argument("--seed", type=int, default=0, help="which random cells are drawn when a stratum exceeds the cap")
    ap.add_argument("--norm", default="panel", choices=["panel", "total"])
    ap.add_argument("--lanes", default="october", choices=["october", "all"], help="blood lanes: the nine October lanes, or all 28")
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--min-ids", type=float, default=None)
    ap.add_argument("--min-genes", type=int, default=P.PROGRAM_SIZE_MIN)
    ap.add_argument("--max-genes", type=int, default=P.PROGRAM_SIZE_MAX)
    ap.add_argument("--budget", type=float, default=None,
                    help="stop before starting a new stratum after this many seconds (continue with --reuse)")
    a = ap.parse_args()
    if a.lineage and a.celltypes:
        ap.error("use --lineage (coarse lineages of the catalogue data) or --celltypes (exact labels), not both")
    strata = D.dataset_strata(a.datasets.split(",")) if a.datasets else D.STRATA
    which = [t.strip() for t in a.celltypes.split(",")] if a.celltypes else a.lineage
    if a.celltypes:
        have = set().union(*[set(units.cells(t).subtype.unique()) for t in {t for t, _ in strata.values()}])
        if set(which) - have:
            ap.error(f"--celltypes: no cells labelled {sorted(set(which) - have)}. Labels in the data: {sorted(have)}")
    if a.out:
        out = a.out
    elif a.datasets:
        if a.celltypes:
            ap.error("--celltypes needs --out, so that runs on different cell types do not overwrite each other")
        out = os.path.join(paths.IDS, "ids-" + a.datasets.replace(",", "-"))
    else:
        out = paths.LINEAGE_RUNS[a.lineage] if a.lineage else paths.DISEASE_RUN
    os.makedirs(out, exist_ok=True)
    budget = a.budget if a.budget is not None else float("inf")
    print(f"strata: {', '.join(sorted(strata))}\nwriting to {out}", flush=True)
    json.dump({"datasets": sorted({t for t, _ in strata.values()}), "strata": {k: list(v) for k, v in strata.items()},
               "lineage": a.lineage, "celltypes": a.celltypes, "lanes": a.lanes, "n_hvg": a.n_hvg,
               "cell_cap": a.cell_cap, "min_cells": a.min_cells, "min_cells_frac": a.min_cells_frac, "seed": a.seed,
               "norm": a.norm, "curated": not a.no_curated, "min_ids": a.min_ids,
               "min_genes": a.min_genes, "max_genes": a.max_genes}, open(os.path.join(out, "run_config.json"), "w"), indent=1)
    O = lambda *p: os.path.join(out, *p)
    lanes = "all" if a.lanes == "all" else None
    t0 = time.time()

    pj = O("gene_panel.json")
    if a.reuse and os.path.exists(pj):
        panel = json.load(open(pj))["panel"]
        print(f"  [panel] reusing {len(panel):,} genes from {pj}", flush=True)
    else:
        panel = D.build_panel(which, a.n_hvg, out, lanes=lanes, strata=strata,
                              curated=[] if a.no_curated else D.CURATED)

    log = json.load(open(O("ids_runlog.json"))) if (a.reuse and os.path.exists(O("ids_runlog.json"))) else {}
    rows = []
    for s in sorted(strata):
        cC, cg = O(f"ids_{s}.npy"), O(f"ids_{s}_genes.json")
        if a.reuse and os.path.exists(cC):
            C, genes = np.load(cC).astype(np.float64), json.load(open(cg))
            n_cells = log.get(s, {}).get("cells", -1)
        else:
            if time.time() - t0 > budget:
                print("[budget] stopping; re-run with --reuse to continue"); return
            r = D.compute_ids(which, s, panel, seed=a.seed, norm=a.norm, cap=a.cell_cap, floor=a.min_cells_frac,
                              min_cells=a.min_cells, block_size=a.block_size, lanes=lanes, strata=strata)
            if "skipped" in r:
                print(f"[{s}] only {r['skipped']:,} cells - below --min-cells {a.min_cells:,}; skipped", flush=True)
                continue
            C, genes, n_cells = r["C"], r["genes"], r["cells"]
            print(f"[{s}] detection floor {a.min_cells_frac:.0%}: {len(panel):,} -> {len(genes):,} genes; "
                  f"IDS on {n_cells:,} cells", flush=True)
            np.save(cC, C.astype(np.float32)); json.dump(genes, open(cg, "w"))
            json.dump(r["dropped"], open(O(f"genes_below_floor_{s}.json"), "w"))
            log[s] = {"stratum": s, "cells": n_cells}; json.dump(log, open(O("ids_runlog.json"), "w"))

        off = C[~np.eye(len(C), dtype=bool)]
        q = {f"p{v}": float(np.percentile(off, v)) for v in (50, 90, 99, 99.9)}
        r = {"stratum": s, "cells": n_cells, "genes": len(genes),
             **{k: round(v, 4) for k, v in q.items()}, "max": round(float(off.max()), 4)}
        rows.append(r)
        print(f"  [{s}] IDS off-diagonal " + "  ".join(f"{k}={v:.3f}" for k, v in q.items())
              + f"  max={off.max():.3f}", flush=True)

        if a.min_ids is not None:
            tbl, dropped = D.programs_from(C, genes, a.min_ids, a.min_genes, a.max_genes, with_dropped=True)
            tbl.insert(0, "stratum", s)
            # how many cells actually support each program: a high median IDS says the genes move
            # together, not how many cells that claim rests on
            tbl = D.support(which, s, tbl, lanes=lanes, strata=strata)
            tbl.to_csv(O(f"programs_{s}.csv"), index=False)
            pd.DataFrame([{"stratum": s, "reason": k, "n_genes": n, "first_genes": ", ".join(g)}
                          for k, items in dropped.items() for n, g in items]
                         ).to_csv(O(f"programs_dropped_{s}.csv"), index=False)
            r["n_programs"], r["n_dropped_large"] = len(tbl), len(dropped["too_large"])
            print(f"  [{s}] {len(tbl)} programs kept ({len(dropped['too_small'])} too small, "
                  f"{len(dropped['too_large'])} too large)", flush=True)
        del C

    df = pd.DataFrame(rows)
    df.to_csv(O("ids_summary.csv"), index=False)
    json.dump({r["stratum"]: r for r in rows}, open(O("ids_runlog.json"), "w"))
    print(); print(df.to_string(index=False), flush=True)
    if a.min_ids is None:
        print("\n--min-ids not given. Read it off the percentiles above (above the background bulk, below the\n"
              "program mode), then re-run with --reuse --min-ids <value>.")


if __name__ == "__main__":
    main()
