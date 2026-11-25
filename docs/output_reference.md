# Output files: what each one is and what its columns mean

Reference for the files written by the three new-dataset commands (`docs/new_dataset_guide.md`).
They are in `results/02_ids_programs/ids-<dataset>/`, or in the folder given with `--out` or `--results-dir`. The same column definitions apply to the
catalogue runs for blood and skin.

| File | Written by | One row is | Use it to |
|---|---|---|---|
| `programs_<dataset>_<disease>.csv` | step 2 | one program found in that disease group | see which genes go together in that group and how tightly |
| `detections.csv` | step 3 | one program found in one disease group, all groups stacked | see which programs are the same one found in several groups |
| `activity_by_celltype.csv` | step 3 | one program in one disease group and cell type | see where each program is active |
| `activity_by_donor.csv.gz` | step 3 | the same, per patient | compare disease groups across patients |

Smaller files are listed at the end.

## How a program is built

Three steps turn pairs of genes into programs and their statistics.

1. **Pairs to links.** IDS is computed for every pair of genes, on at most 4,000 cells of the disease group.
   It is a number from 0 to 1: how strongly the two genes go on and off together across cells. A pair is
   **linked** when its IDS is at least the threshold `--min-ids`.
2. **Links to programs.** The linked genes form a network. It is split into groups with many links inside
   and few links between (greedy modularity). Each group of 3 to 50 genes is a program. A gene does not have
   to be linked to every other gene of its program.
3. **Programs to statistics.** For each program, the IDS of every pair of its genes is looked up again,
   including pairs below the threshold. The statistics summarise that list.

Example: a program with 4 genes has 6 pairs. Threshold 0.5.

| Pair | IDS | Linked? |
|---|---|---|
| A–B | 0.70 | yes |
| A–C | 0.60 | yes |
| B–C | 0.55 | yes |
| C–D | 0.52 | yes |
| A–D | 0.30 | no |
| B–D | 0.20 | no |

`median_ids` 0.535, `mean_ids` 0.48, `q25_ids` 0.36, `min_pair_ids` 0.20, `density` 4 of 6 = 0.67.
Gene D is in the program only through its link to C.

## `programs_<dataset>_<disease>.csv`

The programs found in one disease group, ranked.

| Column | Meaning |
|---|---|
| `stratum` | The group the programs were found in: dataset name, underscore, disease value. |
| `rank` | Order by `median_ids`, highest first (ties: higher `density` first). |
| `program` | `P` + rank + the first gene alphabetically. The gene in the name is not the most important one. |
| `n_genes` | Number of genes in the program (3 to 50). |
| `median_ids` | Middle value of the IDS over all gene pairs of the program. Higher = the genes go together more tightly. |
| `q25_ids` | A quarter of the gene pairs are below this value. |
| `mean_ids` | Average over the gene pairs. |
| `min_pair_ids` | The weakest gene pair. |
| `density` | Share of gene pairs that are linked (IDS at least `--min-ids`). 1 = every gene is linked to every other gene. |
| `recurrence` | Empty. Reserved for the number of groups a program recurs in; not filled by this step. |
| `genes` | The genes, alphabetical. |
| `cells_all_genes` | Cells of the whole disease group in which every gene of the program is detected (count above 0). |
| `cells_half_genes` | Cells in which at least half of the genes are detected (at least 2). |
| `frac_cells_half_genes` | `cells_half_genes` divided by `stratum_cells`. |
| `stratum_cells` | All cells of the disease group. |

How to read a row:

- **Tight or loose.** High `density` and a `min_pair_ids` near the threshold: every gene goes with every other.
  Low `density` and a low `min_pair_ids`: many genes are connected only through other genes, so the gene list
  deserves less confidence. Large programs are usually loose.
- **How many cells it rests on.** Use `cells_half_genes`. `cells_all_genes` is too strict for large programs:
  a 28-gene program can have a handful of cells with all 28 genes detected and thousands with half of them.
- **The IDS columns and the cell columns are separate.** The IDS columns come from the 4,000 sampled cells.
  The cell columns count over every cell of the group.

## `detections.csv`

The `programs_*.csv` files of every disease group stacked into one table, with three columns added in front.
One row is one **detection**: one program found in one disease group.

| Column | Meaning |
|---|---|
| `program_id` | `stratum`, a vertical bar, then `program` (for example `mydata_RA\|P01_C1QA`). Unique per row. It is the key used in the activity tables. |
| `distinct_program` | Detections with the same label are treated as the same program. Two detections are joined when they share at least half of their genes: shared genes divided by all genes in either list (Jaccard) at least 0.5. Joining is chained: if A joins B and B joins C, all three get one label. |
| `found_in` | The disease groups in which that distinct program was found. |
| other columns | As in `programs_<dataset>_<disease>.csv`. |

Cautions:

- **"Found in one group" does not mean absent in the others.** It means the genes did not form a group above
  the threshold there. Check `activity_by_celltype.csv`: if the program scores the same in the same cell type
  in every disease group, it is present in all of them.
- **A group with few cells of a cell type cannot find that cell type's programs.** A plasma-cell program will
  be "missing" from a group that has almost no plasma cells.
- **The half-overlap rule keeps related programs apart.** A 13-gene program that sits entirely inside a
  28-gene program overlaps by 13 of 28, below half, and stays separate. The number of distinct programs is
  therefore an upper count. Compare the gene lists of same-named detections before calling a program
  specific to one disease.
- **This is the unfiltered list.** The catalogue's later steps (removing technical programs such as
  mitochondrial and ribosomal ones, the permutation test, the activity cutoff) are not applied.

## `activity_by_celltype.csv`

Every detection scored in every cell of the dataset, then summarised by disease group and cell type.
All cells are used, not the 4,000, and not only the group where the program was found.

Score of a cell = mean log-expression of the program's genes minus the mean of 50 control genes with matched
expression. 0 = no different from background.

| Column | Meaning |
|---|---|
| `dataset` | Dataset name. |
| `program_id` | The detection, as in `detections.csv`. |
| `disease` | Disease group of the cells summarised in this row (not where the program was found). |
| `cell_type` | Cell type of those cells, from the metadata column given at build time. |
| `n_cells` | Cells of that disease group and cell type. |
| `n_donors` | Patients those cells come from. |
| `mean_score` | Average score of those cells. |
| `prop_expressing` | Share of those cells with a score above `--cutoff`. 0.80 = 80 of 100 cells. |

How to use it:

- **What cell type is this program?** Filter to one `program_id`, sort by `mean_score`.
- **Is it different between disease groups?** Compare rows with the same `program_id` and `cell_type`.

Cautions:

- **The default cutoff is 0, which is loose.** In the psoriasis data, unrelated random genes counted as
  "expressed" in about 30 of 100 cells at cutoff 0 and in 1 to 8 at cutoff 0.2
  (`results/05_cell_proportions/`). Prefer `mean_score`, or run step 3 with `--cutoff 0.2`.
- **Check `n_cells` and `n_donors` first.** A row with 3 cells from 1 patient describes one person.
  `n_donors` = 1 means disease and patient cannot be told apart in that row.
- **Cells are pooled across patients.** A patient with many cells counts more. For comparisons across
  patients use `activity_by_donor.csv.gz`.
- **No test is applied.** The table describes; it does not say whether a difference is real.

## `activity_by_donor.csv.gz`

Same columns as `activity_by_celltype.csv` without `n_donors`, plus `donor`: one row per program, disease
group, cell type and patient. `n_cells`, `mean_score` and `prop_expressing` are for that patient's cells only.

## Smaller files

| File | Content |
|---|---|
| `ids_summary.csv` | One row per disease group. `cells`: cells used for IDS (at most 4,000). `genes`: genes tested after the detection floor. `p50`, `p90`, `p99`, `p99.9`, `max`: spread of IDS over all gene pairs, used to judge the threshold. `n_programs`: programs kept. `n_dropped_large`: gene groups over 50 genes. |
| `programs_dropped_<stratum>.csv` | Gene groups left out for size. `reason`: `too_small` (under 3 genes) or `too_large` (over 50). `n_genes`, `first_genes`. |
| `gene_panel.csv` | The genes offered to every group. `source`: `HVG` (highly variable in at least one group) or `curated` (from the fixed marker list). |
| `genes_below_floor_<stratum>.json` | Panel genes dropped in that group because under 2% of its sampled cells detect them. |
| `run_config.json` | The settings of the run: datasets, groups, cell cap, thresholds. |
| `ids_<stratum>.npy`, `ids_<stratum>_genes.json` | The IDS matrix and its gene order. Kept so the threshold can be changed with `--reuse` without recomputing. Not committed. |
