# Adding a new dataset

For a dataset that comes as one count matrix plus a metadata table.
Nothing about the dataset is written in code: we give the folder, the metadata file and three column names. Run from the repo root.

## What you need

| File                                         | What it must hold                                                                                                                                                                                                                 |
| -------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `matrix.mtx`, `barcodes.tsv`, `features.tsv` | One matrix for all samples, 10x layout (genes × cells), raw counts. Plain or `.gz`. `genes.tsv` is accepted in place of `features.tsv`.                                                                                           |
| `metadata.csv`                               | One row per cell. The barcodes are in the first column (where R's `write.csv` puts row names) and are the same text as in `barcodes.tsv`, including the `-[sample index]` suffix. Three more columns: sample, cell type, disease. Expected in the same folder as the matrix; `.tsv` and `.txt` (tab-separated) are also accepted. |

If the file names share a prefix (`syn_matrix.mtx`, `syn_barcodes.tsv`, `syn_features.tsv`, `syn_metadata.csv`), give it
with `--prefix syn_`, or put it at the end of the folder path: `--matrix-dir /path/to/folder/syn_`.

## The three commands

For a dataset with the sample in `patientID`, the cell type in `simplified_subtype` and the disease in `disease`
(PsD, RA, Control). `psd_ra` is a name you choose; it is used for folder and file names.

```bash
# 1. Build: one file per patient. No QC. metadata.csv is read from the same folder.
python scripts/02_ids_programs/build_units.py --dataset psd_ra --matrix-dir /path/to/folder \
    --sample-col patientID --celltype-col simplified_subtype --disease-col disease

# 2. Find gene programs separately in each disease group (here PsD, RA, Control).
python scripts/02_ids_programs/ids_programs_disease.py --datasets psd_ra --min-ids 0.5

# 3. Score every program in every cell and summarise by cell type, disease and patient.
python scripts/02_ids_programs/program_activity.py --run results/02_ids_programs/ids-psd_ra
```

The three column names shown are also the defaults, so they can be left out. If the metadata is somewhere
else, add `--metadata /other/place/metadata.csv`.

## Where the results go

By default everything is written inside the repo: the dataset to `results/00_units/` and the programs to
`results/02_ids_programs/ids-psd_ra/`. Two options change that.

| Option | On | What it does |
|---|---|---|
| `--results-dir DIR` | steps 1 and 2 | Writes under `DIR` instead of the repo: the dataset to `DIR/00_units/`, the programs to `DIR/02_ids_programs/ids-psd_ra/`. Give the same folder to both steps. |
| `--out DIR` | step 2 | The folder of this one run. Use a different folder for each run you want to keep, for example one per threshold. Without it, a second run on the same dataset writes into the same folder as the first. |

Step 3 needs neither: `--run` points at a run folder, and the run folder remembers where its dataset is.

```bash
python scripts/02_ids_programs/build_units.py --dataset psd_ra --matrix-dir /path/to/folder --results-dir /proj/ids
python scripts/02_ids_programs/ids_programs_disease.py --datasets psd_ra --results-dir /proj/ids \
    --out /proj/ids/02_ids_programs/psd_ra_minids050 --min-ids 0.5
python scripts/02_ids_programs/ids_programs_disease.py --datasets psd_ra --results-dir /proj/ids \
    --out /proj/ids/02_ids_programs/psd_ra_minids045 --min-ids 0.45
python scripts/02_ids_programs/program_activity.py --run /proj/ids/02_ids_programs/psd_ra_minids050
```

The dataset is built once and shared by every run. Each run folder holds a `run_config.json` with its settings.

## What each step writes, and what to check

**Step 1** writes `results/00_units/psd_ra/<patient>.h5ad` and two index tables next to that folder.
It prints how many cells of the matrix were matched to a metadata row, then cells and patients per disease value.
Check that the matched count is the number of cells you expect. If no barcode matches, it stops and shows
three barcodes from each file so the difference is visible. Cells with no metadata row or with a missing
label are left out and counted.

**Step 2** writes to `results/02_ids_programs/ids-psd_ra/`:

| File                                    | Content                                                                                |
| --------------------------------------- | -------------------------------------------------------------------------------------- |
| `ids_summary.csv`                       | Per disease group: cells used, genes tested, the spread of IDS values, programs found. |
| `programs_psd_ra_<disease>.csv`         | The programs found in that group: genes, median IDS, how many cells detect the genes.  |
| `programs_dropped_psd_ra_<disease>.csv` | Gene groups left out for size (under 3 or over 50 genes).                              |
| `gene_panel.csv`                        | The genes that were tested (same list for every group).                                |

The edge threshold `--min-ids` is the one setting that is chosen by hand. Genes are linked when their IDS is at
least this value. 0.5 is what the whole-tissue catalogue runs used and is a reasonable first look. To try another
value without recomputing anything: add `--reuse`, e.g. `--reuse --min-ids 0.45` (seconds). Signs it is too low:
many groups in `programs_dropped_*` are "too_large". Too high: few programs, most with 3 genes.

**Step 3** writes to the same folder:

| File                       | Content                                                                                                                                                                                                                      |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `detections.csv`           | Every program found, one row each. `distinct_program` gives the same label to detections that share at least half of their genes, i.e. the same program found in more than one disease group; `found_in` lists those groups. |
| `activity_by_celltype.csv` | Program × disease × cell type: cells, patients, mean score, proportion of cells expressing the program.                                                                                                                      |
| `activity_by_donor.csv.gz` | The same per patient.                                                                                                                                                                                                        |

A cell's score is the mean log-expression of the program's genes minus the mean of 50 control genes with matched
expression. "Expressing" means a score above 0; change it with `--cutoff`.

What every column in these files means, and how to read them: [`output_reference.md`](output_reference.md).

## Optional

Programs inside one group of cell types (the equivalent of the T+NK, myeloid and keratinocyte runs), using the
dataset's own labels:

```bash
python scripts/02_ids_programs/ids_programs_disease.py --datasets psd_ra --celltypes "CD4 T,CD8 T,NK" \
    --out results/02_ids_programs/ids-psd_ra-TNK --min-ids 0.4
```

Score an existing program list in the new dataset, for example the 57 distinct programs of the psoriasis catalogue:

```bash
python scripts/02_ids_programs/program_activity.py --programs results/03_catalogue/distinct_programs.csv \
    --datasets psd_ra --out results/02_ids_programs/catalogue-in-psd_ra
```

QC, if it is ever wanted: add `--min-genes 200 --max-mt 10` to step 1.

## Settings that are applied

| Setting                                 | Value                                                                                          | Change with                  |
| --------------------------------------- | ---------------------------------------------------------------------------------------------- | ---------------------------- |
| Cells used to find programs             | at most 4,000 per disease group, drawn evenly across cell types                                | `--cell-cap`                 |
| Smallest disease group that is analysed | 400 cells (smaller groups are skipped and reported)                                            | `--min-cells`                |
| Genes tested                            | 1,200 most variable genes per group, pooled, plus a fixed list of immune and skin marker genes | `--n-hvg`, `--no-curated`    |
| Genes dropped per group                 | detected in under 2% of that group's sampled cells                                             | `--min-cells-frac`           |
| Program size                            | 3 to 50 genes                                                                                  | `--min-genes`, `--max-genes` |

Scoring (step 3) uses every cell, not the 4,000.

## Limits

- **The file is per patient.** If one patient has several samples, their cells go into one file.
- **Raw counts are expected.** Step 1 warns if the matrix holds non-whole numbers.

## How this path was checked

- 29,034 skin cells (15 samples) were exported to exactly this format, with the cells in shuffled order, and rebuilt
  with step 1: same cells, counts and labels as the original per-sample files.
- Step 2 run on the skin data through this path gives the catalogue's skin programs: same genes, values equal to
  14 decimal places. The default run for blood and skin is unchanged.
- Step 3 run on all 67,367 skin cells gives the catalogue's cell counts exactly and its mean scores to within 0.0000001.
- `python tests/test_new_dataset.py` builds a small made-up dataset (PsD, RA, Control) and checks the cells, counts,
  labels and error messages.
