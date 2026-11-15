"""A new dataset from matrix.mtx + barcodes.tsv + features.tsv + metadata.csv (sc_ids.units.build_dataset).

Makes a small dataset in that format, builds it, and checks that every kept cell comes back with
its own counts and labels, and that the conditions are read from the data.

  python tests/test_new_dataset.py
"""
import os, sys, tempfile
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io

tmp = tempfile.mkdtemp()
os.environ["SC_IDS_WORK"] = os.path.join(tmp, "units")          # must be set before sc_ids is imported
os.makedirs(os.environ["SC_IDS_WORK"])
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from sc_ids import units, discovery                              # noqa: E402


def make(folder, n_cells=600, n_genes=80, seed=0):
    rng = np.random.default_rng(seed)
    X = sp.random(n_cells, n_genes, density=0.3, random_state=seed, format="csr")
    X.data = np.ceil(X.data * 9)                                               # whole-number counts
    patients = ["007", "P2", "P3", "P4", "P5", "P 6"]                          # a numeric-looking name and one with a space
    disease = dict(zip(patients, ["PsD", "PsD", "RA", "RA", "Control", "Control"]))
    pat = rng.choice(patients, n_cells)
    barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-{patients.index(p) + 1}" for p in pat]
    genes = [f"G{i}" for i in range(n_genes)]; genes[1] = "G0"                 # a repeated gene symbol
    feat = pd.DataFrame({0: [f"ENSG{i:05d}" for i in range(n_genes)] + ["ab1", "ab2"], 1: genes + ["CD3", "CD4"],
                         2: ["Gene Expression"] * n_genes + ["Antibody Capture"] * 2})
    full = sp.hstack([X, sp.csr_matrix(np.ones((n_cells, 2)) * 500)]).T        # genes x cells, with 2 antibody rows
    scipy.io.mmwrite(os.path.join(folder, "matrix.mtx"), sp.coo_matrix(full), field="integer")
    pd.Series(barcodes).to_csv(os.path.join(folder, "barcodes.tsv"), index=False, header=False)
    feat.to_csv(os.path.join(folder, "features.tsv"), sep="\t", index=False, header=False)
    meta = pd.DataFrame({"patientID": pat, "simplified_subtype": rng.choice(["T cell", "B cell", "Monocyte"], n_cells),
                         "disease": [disease[p] for p in pat], "other": 1}, index=barcodes)
    meta = meta.iloc[4:].copy()                                                # cells 0-3 are absent from the metadata
    meta.iloc[:3, 2] = np.nan                                                  # cells 4-6 have no disease label
    meta.loc["TTTTTTTTTTTTTTTT-9"] = ["P9", "T cell", "RA", 1]                 # a metadata row with no cell
    meta.sample(frac=1, random_state=1).to_csv(os.path.join(folder, "metadata.csv"))   # R-style: barcodes are the unnamed first column
    return X, barcodes, pat, disease


def main():
    src = os.path.join(tmp, "raw"); os.makedirs(src)
    X, barcodes, pat, disease = make(src)
    c = units.build_dataset("toy", src, os.path.join(src, "metadata.csv"))
    kept = np.arange(7, len(barcodes))                                         # 4 cells without metadata, 3 without a disease label
    assert len(c) == len(kept), (len(c), len(kept))
    assert set(c.cell) == {barcodes[i] for i in kept}
    assert set(c.unit) == {"007", "P2", "P3", "P4", "P5", "P_6"} and set(c.donor) == set(pat)
    assert units.genes("toy")[:2] == ["G0", "G0-1"] and len(units.genes("toy")) == X.shape[1]      # antibody rows gone
    got = units.read_cells("toy", c)                                           # counts in table order
    pos = {b: i for i, b in enumerate(barcodes)}
    want = X[[pos[b] for b in c.cell]]
    assert abs(got - want).sum() == 0
    assert (c.disease_fine.values == [disease[p] for p in c.donor]).all()
    assert np.allclose(c.total_counts.values, np.asarray(want.sum(axis=1)).ravel())
    s = discovery.dataset_strata(["toy"])
    assert s == {"toy_Control": ("toy", "Control"), "toy_PsD": ("toy", "PsD"), "toy_RA": ("toy", "RA")}, s
    _, cc = discovery.stratum_cells(["T cell", "B cell"], "toy_RA", strata=s)
    assert set(cc.disease_fine) == {"RA"} and set(cc.subtype) == {"T cell", "B cell"}
    for bad, kw in (("no column", dict(sample_col="patient")), ("matches the metadata", dict(barcode_col="patientID2"))):
        m = pd.read_csv(os.path.join(src, "metadata.csv")); m["patientID2"] = m["patientID"] + "x" + m.index.astype(str)
        m.to_csv(os.path.join(src, "metadata2.csv"), index=False)
        try:
            units.build_dataset("toy2", src, os.path.join(src, "metadata2.csv"), **kw); raise SystemExit("should have failed")
        except ValueError as e:
            assert bad in str(e), e
    import gzip, shutil                                                        # the same files gzipped give the same dataset
    for f in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
        with open(os.path.join(src, f), "rb") as i, gzip.open(os.path.join(src, f + ".gz"), "wb") as o:
            shutil.copyfileobj(i, o)
        os.remove(os.path.join(src, f))
    cz = units.build_dataset("toy_gz", src, os.path.join(src, "metadata.csv"))
    assert cz.equals(c) and abs(units.read_cells("toy_gz", cz) - got).sum() == 0
    print("new-dataset build: OK")


if __name__ == "__main__":
    main()
