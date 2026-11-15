"""
program_proportions.py — what share of each cell type's cells express each program, by disease state.

Follow-up to the November meeting: instead of asking only whether a program's MEAN score differs
between disease states, ask how many cells of a cell type express it ("dominance, not presence").

For each of the 57 distinct programs (scored on the union of its member detections' genes):
  1. score every cell, lane by lane (sc_ids.units; same score as everywhere else)
  2. per donor and cell type: cells, mean score, and cells with a score above each cutoff in CUTOFFS
     (a donor x cell type x lane needs >= 20 cells, then lanes are pooled; same rule as the disease tests)
  3. compare disease states across donors, only where both states share lanes:
       blood  Healthy vs PSO in runs 02-05, PSO vs PsA in runs 06-07
       skin   NS vs PP, NS vs PN, and PN vs PP paired within donor
     Mann-Whitney (paired Wilcoxon for PN vs PP), Benjamini-Hochberg within comparison and measure.

Two cutoffs are carried side by side because the choice matters and neither is obviously right.

RANDOM GENE SETS. Every step is repeated on random gene sets (N_RANDOM per program, same number of genes,
drawn from genes of matched expression level). They answer two questions the programs alone cannot:
how many cells "express" a set of unrelated genes at each cutoff, and how many disease differences a
set of unrelated genes shows. A program result only means something where it exceeds that.

  python scripts/05_cell_proportions/program_proportions.py            # checkpointed; re-run until "done"

Writes results/05_cell_proportions/:
  donor_proportions.csv.gz    comparison, cell type, donor, disease, program: cells, mean score, share at each cutoff
  proportion_tests.csv        every test: group means, difference, p, adjusted p, for mean score and each cutoff
  summary.csv                 per comparison and measure: tests, significant, significant with a difference >= MIN_DIFF,
                              for the 57 programs and for the random gene sets
  random_shares.csv           tissue x cell type: share of cells expressing a random gene set, at each cutoff
  celltype_shares.csv         tissue x program x cell type: share of cells expressing, all donors pooled
  bars_<comparison>.pdf       grouped barplots, one page per cell type, one row per cutoff
  heatmap_<tissue>.png        program x cell type share, one panel per cutoff
"""
import argparse, os, sys, time
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "04_method_checks"))
from sc_ids import paths, units                                               # noqa: E402
from sc_ids.lineage import LINEAGE                                           # noqa: E402
from disease_tests import contrast, lanes                                     # noqa: E402

OUT = os.path.join(paths.RESULTS, "05_cell_proportions")
CUTOFFS = (0.0, 0.2)
MIN_CELLS = 20
N_RANDOM = 5
MIN_DIFF = {"mean_score": 0.2, "share_0": 0.10, "share_0.2": 0.10}            # smallest difference called a change
MEASURES = list(MIN_DIFF)
PLAN = [  # tissue, lanes label, lanes kept, state a, state b, paired
    ("blood", "runs 02-05 (16 lanes)", lanes(2, 3, 4, 5), "Healthy", "PSO", False),
    ("blood", "runs 06-07 (8 lanes)", lanes(6, 7), "PSO", "PSA", False),
    ("skin", "all samples", None, "NS", "PP", False),
    ("skin", "all samples", None, "NS", "PN", False),
    ("skin", "donors with both", None, "PN", "PP", True),
]
SURFACE, INK, MUTED = "#fcfcfb", "#0b0b0b", "#52514e"
COLOR = {"Healthy": "#2a78d6", "NS": "#2a78d6", "PSO": "#eb6834", "PP": "#eb6834", "PSA": "#1baf7a", "PN": "#1baf7a"}
NAME = {"Healthy": "Healthy", "PSO": "Psoriasis", "PSA": "PsA", "NS": "Healthy skin (NS)", "PN": "Non-lesional (PN)", "PP": "Lesional (PP)"}


def programs():
    d = pd.read_csv(os.path.join(paths.CATALOGUE, "distinct_programs.csv"))
    return {r.program_no: [g.strip() for g in r.genes.split(",")] for r in d.itertuples()}, d


def random_sets(tissue, progs, n=N_RANDOM, seed=0):
    """n random gene sets per program: same size, drawn from that program's own matched control genes."""
    g = np.array(units.genes(tissue))
    mean, _ = units.gene_means(tissue, units.stage_units(tissue, "all"),
                               cache=os.path.join(units.WORK, "scores", f"proportions_{tissue}", "_gene_means.npz"))
    rng = np.random.default_rng(seed)
    return {f"{pid}|r{r}": list(g[rng.choice(ci, size=len(gi), replace=False)])
            for pid, (gi, ci) in units.control_sets(g, mean, progs).items() for r in range(n)}


def score(tissue, progs, budget, t0, tag=""):
    """Per lane x donor x subtype x disease: cells, summed score and cells above each cutoff, for every program."""
    us = units.stage_units(tissue, "all")
    tmp = os.path.join(units.WORK, "scores", f"proportions_{tissue}"); os.makedirs(tmp + tag, exist_ok=True)
    mean, _ = units.gene_means(tissue, us, cache=os.path.join(tmp, "_gene_means.npz"))
    ctrl = units.control_sets(units.genes(tissue), mean, progs)
    tmp += tag
    parts = []
    for u in us:
        f = os.path.join(tmp, f"{u}.pkl")
        if not os.path.exists(f):
            if time.time() - t0 > budget:
                return None
            a = units.read_unit(tissue, u)
            S, ids = units.score_matrix(units.lognorm(a.X), ctrl)
            o = a.obs
            code, grp = pd.factorize(pd.MultiIndex.from_arrays([o["donor"].astype(str), o["subtype"].astype(str), o["disease_fine"].astype(str)]))
            M = sp.csr_matrix((np.ones(len(code)), (code, np.arange(len(code)))), shape=(len(grp), len(code)))
            idx = pd.MultiIndex.from_tuples(list(grp), names=["donor", "subtype", "disease"])
            cols = {"sum": pd.DataFrame(M @ S, index=idx, columns=ids).stack()}
            for c in CUTOFFS:
                cols[f"k_{c:g}"] = pd.DataFrame(M @ (S > c).astype(float), index=idx, columns=ids).stack()
            L = pd.concat(cols, axis=1); L.index.names = ["donor", "subtype", "disease", "program_id"]
            L["n"] = np.repeat(np.asarray(M.sum(axis=1)).ravel(), len(ids))
            L = L.reset_index(); L.insert(0, "unit", u)
            L.to_pickle(f + ".tmp"); os.replace(f + ".tmp", f)
        parts.append(pd.read_pickle(f))
    return pd.concat(parts, ignore_index=True)


def donor_table(L, keep):
    """One row per donor x cell type x disease x program, lanes pooled (>= MIN_CELLS cells per lane first)."""
    L = L if keep is None else L[L.unit.isin(keep)]
    L = L.assign(group=L.subtype.map(LINEAGE).fillna("Other"))
    L = L[L.group != "Other"]
    val = ["n", "sum"] + [f"k_{c:g}" for c in CUTOFFS]
    a = L.groupby(["unit", "donor", "group", "disease", "program_id"], observed=True)[val].sum()
    a = a[a.n >= MIN_CELLS]
    D = a.groupby(["donor", "group", "disease", "program_id"], observed=True).sum().reset_index()
    D["mean_score"] = D["sum"] / D.n
    for c in CUTOFFS:
        D[f"share_{c:g}"] = D[f"k_{c:g}"] / D.n
    return D.rename(columns={"n": "n_cells"})[["donor", "group", "disease", "program_id", "n_cells"] + MEASURES]


def figure_bars(D, R, tissue, a_, b_, label, ids, pdf):
    import matplotlib.pyplot as plt
    x = np.arange(len(ids)); w = 0.4
    for grp in sorted(R.group.unique()):
        s = D[D.group == grp]
        fig, axes = plt.subplots(len(CUTOFFS), 1, figsize=(15, 6.6), sharex=True, facecolor=SURFACE)
        for ax, c in zip(axes, CUTOFFS):
            m = f"share_{c:g}"
            r = R[(R.group == grp) & (R.measure == m)].set_index("program_id").reindex(ids)
            for k, dz in enumerate((a_, b_)):
                P = s[s.disease == dz].pivot(index="donor", columns="program_id", values=m).reindex(columns=ids) * 100
                if grp in set(R.group) and R.attrs.get("paired"):
                    both = set(s[s.disease == a_].donor) & set(s[s.disease == b_].donor); P = P[P.index.isin(both)]
                pos = x + (k - 0.5) * w
                ax.bar(pos, P.mean(), width=w - 0.05, color=COLOR[dz], label=f"{NAME[dz]} ({len(P)} donors)", linewidth=0, zorder=2)
                ax.vlines(pos, P.quantile(0.25), P.quantile(0.75), color=INK, linewidth=0.7, zorder=3)
            sig = (r.p_adj < 0.05) & (r.delta.abs() >= MIN_DIFF[m])
            top = np.fmax(r.mean_a, r.mean_b).values * 100
            for j in np.where(sig.values)[0]:
                ax.text(x[j], min(top[j] + 9, 103), "*", ha="center", va="center", fontsize=13, color=INK)
            ax.set_ylim(0, 108); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_facecolor(SURFACE)
            ax.grid(axis="y", color="#e4e3df", linewidth=0.6, zorder=0); ax.set_axisbelow(True)
            for sp_ in ("top", "right", "left"):
                ax.spines[sp_].set_visible(False)
            ax.spines["bottom"].set_color("#b9b8b2"); ax.tick_params(colors=MUTED, length=0, labelsize=8)
            ax.set_ylabel("Cells expressing, per 100", color=MUTED, fontsize=9)
            ax.set_title(f"Expressing = score above {c:g}      * = adjusted p < 0.05 and a difference of at least 10 cells per 100: "
                         f"{int(sig.sum())} of {int(r.p.notna().sum())} programs", loc="left", fontsize=9.5, color=INK)
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper right", bbox_to_anchor=(0.995, 0.995), ncol=2, frameon=False, fontsize=9.5, labelcolor=INK)
        axes[-1].set_xticks(x); axes[-1].set_xticklabels(ids, rotation=90, fontsize=7.5, color=MUTED)
        axes[-1].set_xlim(-0.8, len(ids) - 0.2)
        fig.suptitle(f"{tissue.capitalize()}, {grp} cells: {NAME[a_]} vs {NAME[b_]}, {label}", x=0.008, y=0.985, ha="left", fontsize=12.5, color=INK, weight="bold")
        fig.text(0.008, 0.005, "Bar = mean over donors of the share of that donor's cells expressing the program. Line = middle half of donors.",
                 fontsize=8, color=MUTED)
        fig.tight_layout(rect=(0, 0.02, 1, 0.96)); pdf.savefig(fig, facecolor=SURFACE); plt.close(fig)


def figure_heatmap(C, tissue, ids):
    import matplotlib.pyplot as plt
    s = C[C.tissue == tissue]
    cols = list(s.groupby("group").n_cells.max().sort_values(ascending=False).index)
    fig, axes = plt.subplots(1, len(CUTOFFS), figsize=(3.2 + 0.52 * len(cols) * len(CUTOFFS), 12.5), facecolor=SURFACE)
    for ax, c in zip(axes, CUTOFFS):
        M = s.pivot(index="program_id", columns="group", values=f"share_{c:g}").reindex(index=ids, columns=cols) * 100
        im = ax.imshow(M.values, cmap="Blues", vmin=0, vmax=100, aspect="auto")
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8, color=MUTED)
        ax.set_yticks(range(len(ids))); ax.set_yticklabels(ids, fontsize=7, color=MUTED)
        ax.tick_params(length=0)
        for sp_ in ax.spines.values():
            sp_.set_visible(False)
        ax.set_xticks(np.arange(-.5, len(cols)), minor=True); ax.set_yticks(np.arange(-.5, len(ids)), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=1.2); ax.tick_params(which="minor", length=0)
        ax.set_title(f"score above {c:g}", fontsize=10, color=INK, loc="left")
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02); cb.outline.set_visible(False)
    cb.set_label("Cells expressing, per 100 (all donors pooled)", color=MUTED, fontsize=9); cb.ax.tick_params(colors=MUTED, length=0, labelsize=8)
    fig.suptitle(f"{tissue.capitalize()}: share of each cell type's cells expressing each program", x=0.02, y=0.915, va="bottom", ha="left", fontsize=12, color=INK, weight="bold")
    fig.savefig(os.path.join(OUT, f"heatmap_{tissue}.png"), dpi=150, facecolor=SURFACE, bbox_inches="tight"); plt.close(fig)


def smallest_p(na, nb, paired):
    from math import comb
    return 2 / 2 ** na if paired else 2 / comb(na + nb, na)


def compare(L, ids, sets, figures=False):
    """Donor tables, tests and one summary row per comparison x measure for the gene sets `ids`."""
    donors, tests, summ = [], [], []
    for tissue, label, keep, a_, b_, paired in PLAN:
        use = [i for i in ids[tissue]]
        D = donor_table(L[tissue], keep); D = D[D.disease.isin([a_, b_])]
        comp = f"{a_} vs {b_}"
        Rs = []
        for m in MEASURES:
            W = D.pivot(index=["donor", "group", "disease"], columns="program_id", values=m).reset_index()
            R = contrast(W, use, a_, b_, paired); R.insert(2, "measure", m); Rs.append(R)
            sig = R[R.p_adj < 0.05]; big = sig[sig.delta.abs() >= MIN_DIFF[m]]
            na, nb = int(R.n_a.max()), int(R.n_b.max())
            summ.append(dict(gene_sets=sets, tissue=tissue, comparison=comp, lanes=label, measure=m, donors_a=na, donors_b=nb,
                             cell_types=R.group.nunique(), tests=len(R), raw_p_below_0_05=int((R.p < 0.05).sum()),
                             significant=len(sig), significant_and_min_difference=len(big),
                             per_100_tests=round(100 * len(big) / len(R), 1), min_difference=MIN_DIFF[m],
                             largest_abs_difference=round(float(R.delta.abs().max()), 3),
                             smallest_possible_p=smallest_p(na, nb, paired)))
        R = pd.concat(Rs, ignore_index=True); R.attrs["paired"] = paired
        if figures:
            from matplotlib.backends.backend_pdf import PdfPages
            with PdfPages(os.path.join(OUT, f"bars_{tissue}_{a_}_vs_{b_}.pdf")) as pdf:
                figure_bars(D, R, tissue, a_, b_, label, use, pdf)
        for T, store in ((D, donors), (R, tests)):
            T = T.copy(); T.insert(0, "lanes", label); T.insert(0, "comparison", comp); T.insert(0, "tissue", tissue); store.append(T)
    return pd.concat(donors, ignore_index=True), pd.concat(tests, ignore_index=True), pd.DataFrame(summ)


def shares(L):
    """Every cell of each tissue, donors and disease states pooled: share expressing, per cell type and gene set."""
    cell = []
    for tissue in L:
        x = L[tissue].assign(group=L[tissue].subtype.map(LINEAGE).fillna("Other"))
        g = x[x.group != "Other"].groupby(["group", "program_id"], observed=True)[["n", "sum"] + [f"k_{c:g}" for c in CUTOFFS]].sum().reset_index()
        g["mean_score"] = g["sum"] / g.n
        for c in CUTOFFS:
            g[f"share_{c:g}"] = g[f"k_{c:g}"] / g.n
        g.insert(0, "tissue", tissue); cell.append(g.rename(columns={"n": "n_cells"})[["tissue", "group", "program_id", "n_cells"] + MEASURES])
    C = pd.concat(cell, ignore_index=True); C["n_cells"] = C.n_cells.astype(int)
    return C


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--budget", type=float, default=130); a = ap.parse_args()
    t0 = time.time(); os.makedirs(OUT, exist_ok=True)
    progs, meta = programs()
    L, Lr, rnd = {}, {}, {}
    for tissue in ("skin", "blood"):
        L[tissue] = score(tissue, progs, a.budget, t0)
        if L[tissue] is not None:
            rnd[tissue] = random_sets(tissue, progs)
            Lr[tissue] = score(tissue, rnd[tissue], a.budget, t0, tag="_random")
        if L[tissue] is None or Lr[tissue] is None:
            print(f"[budget] {tissue} scoring not finished; run again to continue"); return
        print(f"[score] {tissue}: {len(progs)} programs and {len(rnd[tissue])} random gene sets scored ({time.time()-t0:.0f}s)", flush=True)

    import matplotlib; matplotlib.use("Agg")
    ids = {t: list(progs) for t in L}
    D, T, S = compare(L, ids, "57 programs", figures=True)
    _, Tr, Sr = compare(Lr, {t: list(rnd[t]) for t in Lr}, f"{N_RANDOM * len(progs)} random gene sets")
    D.round(5).to_csv(os.path.join(OUT, "donor_proportions.csv.gz"), index=False, compression="gzip")
    T = T.merge(meta[["program_no", "label_rule"]].rename(columns={"program_no": "program_id", "label_rule": "label"}), on="program_id")
    T.round(6).to_csv(os.path.join(OUT, "proportion_tests.csv"), index=False)
    S = pd.concat([S, Sr], ignore_index=True).sort_values(["tissue", "comparison", "measure", "gene_sets"])
    S.to_csv(os.path.join(OUT, "summary.csv"), index=False)
    C = shares(L); C.round(5).to_csv(os.path.join(OUT, "celltype_shares.csv"), index=False)
    Cr = shares(Lr)
    q = Cr[Cr.n_cells >= 200].groupby(["tissue", "group"])[[f"share_{c:g}" for c in CUTOFFS]].describe(percentiles=[.5, .95])
    q = q.loc[:, (slice(None), ["50%", "95%", "max"])]; q.columns = [f"{m}_{k.replace('%', 'th')}" for m, k in q.columns]
    q.round(4).reset_index().to_csv(os.path.join(OUT, "random_shares.csv"), index=False)
    for tissue in L:
        figure_heatmap(C, tissue, list(progs))
    show = ["gene_sets", "tissue", "comparison", "measure", "donors_a", "donors_b", "tests", "significant", "significant_and_min_difference", "per_100_tests", "largest_abs_difference"]
    print(); print(S[show].to_string(index=False)); print(f"\ndone ({time.time()-t0:.0f}s): {OUT}")


if __name__ == "__main__":
    main()
