# Refactor briefing: LC-NE_BARseq_MAT-RDS_conversion

## Goal

Refactor this Code Ocean capsule so that it (a) runs only on the two relevant brains, (b) consumes the new bar-seq-only input data assets, and (c) produces output folders named per AIND derived-asset conventions. After this PR is merged and a new capsule release is cut, two derived assets (one per subject) will be saved from the released capsule's reproducible run.

## Repo

- GitHub: https://github.com/AllenNeuralDynamics/LC-NE_BARseq_MAT-RDS_conversion
- Default branch: `master`
- License: MIT (already present)

You should have write access. If `git push` fails, ask the human to grant access before continuing.

## Workflow (do these in order)

### 1. Create a Code Ocean capsule linked to the original repo

In Code Ocean, click **New Capsule → from Git repository** and point it at `AllenNeuralDynamics/LC-NE_BARseq_MAT-RDS_conversion`. Do NOT use the "Edit this capsule" button on the original capsule — that creates a detached duplicate, not a git-linked working capsule.

### 2. Create a working branch from master

From inside the new capsule's git panel (or via `gh` CLI), create a branch off `master`. Suggested name: `refactor/per-subject-outputs`.

Switch to that branch. All commits below should land on this branch.

### 3. Attach the new input data assets to the capsule

Two new bar-seq input assets need to be mounted:

- `barseq_780345_2025-02-24_12-00-00` (id `a9c7cc08-a4d8-458e-a7ad-47fad643e7a5`)
- `barseq_780346_2025-06-13_12-00-00` (id: ask the human, or look up via Code Ocean's data asset search by name)

Both are in `aind-open-data` already. Each contains a `BARseq/` subdirectory with `combined_neurons_clust_CCFv2.mat`.

Detach the existing mounted assets:
- `BARseq669594` (brain1 — being removed entirely)
- `780345_2025-02-20_00-00-00` (old specimen-level asset bundling MAPseq + BARseq)
- `780346_2025-06-11_00-00-00` (same — old combined asset)

Attach the two new bar-seq-only assets at their default mount names. This will update `.codeocean/datasets.json`. Commit that change with a message like `chore: swap to bar-seq-only input assets`.

### 4. Code changes

#### 4a. Delete the brain1 conversion entirely

- Delete `code/01_BarSeq_RDSconvert_brain1_v2.R`
- In `code/run`, delete the line `render_html 01_BarSeq_RDSconvert_brain1_v2.R` and its preceding comment line

Brain1 (specimen 669594) does not appear in the manuscript and its converted output is not consumed by any downstream capsule. The `BARseq669594` data asset (already detached in step 3) is no longer referenced.

#### 4b. Update brain3 script for new input asset and dynamic output folder

In `code/01_BarSeq_RDSconvert_brain3_v2.R`:

- **Input path** (currently line 14 in `convert()`): change from
  `/data/780345_2025-02-20_00-00-00/BARseq/combined_neurons_clust_CCFv2.mat`
  to
  `/data/barseq_780345_2025-02-24_12-00-00/BARseq/combined_neurons_clust_CCFv2.mat`

- **Output folder**: at the top of the script, after the libraries are loaded, define an output folder variable using the input asset name + process name + timestamp:

  ```r
  INPUT_ASSET <- "barseq_780345_2025-02-24_12-00-00"
  PROCESS_NAME <- "processed_MAT2RDS"
  TIMESTAMP <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  OUT_DIR <- file.path("/results", paste0(INPUT_ASSET, "_", PROCESS_NAME, "_", TIMESTAMP))
  dir.create(OUT_DIR, recursive = TRUE)
  ```

  Then replace every hardcoded `/results/BARseq_780345/` path in the script with `OUT_DIR`. There are 5 such hardcoded paths (3 saveRDS calls and the `dir.create` plus the one `readRDS` for sanity-checking). Use `file.path(OUT_DIR, "combined_neurons_clust_CCFv2.rds")` etc.

- **Latent crash fix** (around line 137 in the current script): brain3 currently has

  ```r
  # Drop the rows with higher indices
  sce <- sce[-rows_to_drop, ]
  ```

  Replace with the guarded version that brain4 already uses:

  ```r
  # Drop the rows with higher indices (only if duplicates exist)
  if (length(rows_to_drop) > 0) {
    sce <- sce[-rows_to_drop, ]
  }
  ```

  Without this guard, the script crashes if there are no duplicate hybridization-cycle genes in the input.

- **Top comment** (line 1): currently reads `######## Brain #4 - 780345 ########` — fix the typo to `Brain #3` so it matches the brain it actually processes.

- **Other comment fix**: line 79 reads `#exectute convert function` — fix to `#execute convert function`.

#### 4c. Update brain4 script for new input asset and dynamic output folder

In `code/01_BarSeq_RDSconvert_brain4_v2.R`:

- **Input path**: change from
  `/data/780346_2025-06-11_00-00-00/BARseq/combined_neurons_clust_CCFv2.mat`
  to
  `/data/barseq_780346_2025-06-13_12-00-00/BARseq/combined_neurons_clust_CCFv2.mat`

- **Output folder**: same pattern as brain3 but with the brain4 asset name:

  ```r
  INPUT_ASSET <- "barseq_780346_2025-06-13_12-00-00"
  PROCESS_NAME <- "processed_MAT2RDS"
  TIMESTAMP <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  OUT_DIR <- file.path("/results", paste0(INPUT_ASSET, "_", PROCESS_NAME, "_", TIMESTAMP))
  dir.create(OUT_DIR, recursive = TRUE)
  ```

  Replace every hardcoded `/results/BARseq_780346/` path with `OUT_DIR` via `file.path(OUT_DIR, ...)`.

- The brain4 script already has the `if (length(rows_to_drop) > 0)` guard, so no fix needed there.

### 5. Repository hygiene

- Delete the `.Trash-0/` directory at the repo root. Add `.Trash-0/` to `.gitignore` so it won't reappear if anyone deletes a file through the Code Ocean UI in the future. (Code Ocean uses `.Trash-0/` as a soft-delete location — files deleted through the UI move there rather than being removed from the underlying git repo.)

### 6. Rewrite the README

Delete `readme.md` and create a new `README.md` (uppercase). Use this content:

```markdown
# LC-NE BARseq MATLAB-to-RDS Conversion

Data preparation capsule that converts BARseq gene-expression data from MATLAB format to R SingleCellExperiment (SCE) objects for downstream analysis, as part of:

> Su, Kosillo, Jung, Chen *et al.* (2026). Topographic structure and function of locus coeruleus norepinephrine neurons. [bioRxiv 2026.04.10.717727](https://www.biorxiv.org/content/10.64898/2026.04.10.717727v1)

This capsule does not produce manuscript figures. Its outputs are saved as two per-subject derived data assets which are consumed by the downstream analysis capsule [LC-NE_BARseq_MAPseq_analyses](https://github.com/AllenNeuralDynamics/LC-NE_BARseq_MAPseq_analyses) ([Code Ocean](https://codeocean.allenneuraldynamics.org/capsule/2195789/tree)), which uses them to generate Figure S5 of the manuscript.

**GitHub:** https://github.com/AllenNeuralDynamics/LC-NE_BARseq_MAT-RDS_conversion  
**Code Ocean:** https://codeocean.allenneuraldynamics.org/capsule/3953531/tree  
**Full collection:** https://codeocean.allenneuraldynamics.org/collections/9cf044ce-93c7-4c7e-bfa1-5d8c37aa42ec

## Code

| File | Description |
|------|-------------|
| `00_env_library_loading.R` | Loads the `r4-base` conda environment and core libraries (`hdf5r`, `Matrix`, `SingleCellExperiment`). Provided as a reference for interactive use; not called directly by the `run` script. |
| `01_BarSeq_RDSconvert_brain3_v2.R` | Converts specimen 780345 (brain 3). |
| `01_BarSeq_RDSconvert_brain4_v2.R` | Converts specimen 780346 (brain 4). |
| `run` | Bash entry point for Reproducible Run. Renders each conversion script to an HTML report via `knitr::spin`. |

Each conversion script performs the same steps:

1. Opens the BARseq MATLAB file (`.mat`, HDF5 v7.3 format) using `hdf5r::H5File` and reads the `filt_neurons` group.
2. Reconstructs the sparse gene-by-cell count matrix from stored CSC components using `Matrix::sparseMatrix`.
3. Extracts per-cell metadata: slice, position, FOV coordinates, angle, depth, barcode status, batch number, CCF coordinates, and CCF annotation.
4. Constructs a unique cell identifier (`uid`) from batch, slice, and cell ID.
5. Assembles into a `SingleCellExperiment` object and saves as `combined_neurons_clust_CCFv2.rds`.
6. Validates uid uniqueness, renames columns by uid, removes placeholder genes (`unused-*`) and duplicate hybridization-cycle genes.
7. Saves the cleaned SCE as `combined_neurons_clust_CCFv2_uid.rds` — this is the file consumed by the downstream analysis capsule.
8. Filters to putative LC-NE neurons (Dbh expression > 2) and saves as `DBHfilteredneurons_clust_CCFv2_uid.rds`.

## Data assets

| Asset | is_public | Description | Source |
|-------|-----------|-------------|--------|
| `barseq_780345_2025-02-24_12-00-00` | true | BARseq data for specimen 780345 (brain 3). Contains `BARseq/combined_neurons_clust_CCFv2.mat`. Bucket: `aind-open-data`. | AIND pipeline |
| `barseq_780346_2025-06-13_12-00-00` | true | BARseq data for specimen 780346 (brain 4). Contains `BARseq/combined_neurons_clust_CCFv2.mat`. Bucket: `aind-open-data`. | AIND pipeline |

## Outputs

Each conversion script writes to a per-subject output folder under `/results/`, named:

```
/results/<input_asset_name>_processed_MAT2RDS_<timestamp>/
```

For example, a run on May 6, 2026 might produce:

- `/results/barseq_780345_2025-02-24_12-00-00_processed_MAT2RDS_2026-05-06_17-30-00/`
- `/results/barseq_780346_2025-06-13_12-00-00_processed_MAT2RDS_2026-05-06_17-30-00/`

Each output folder contains three `.rds` files:

| File | Description | Consumed downstream? |
|------|-------------|---------------------|
| `combined_neurons_clust_CCFv2.rds` | Initial SingleCellExperiment object before duplicate-gene cleanup | No |
| `combined_neurons_clust_CCFv2_uid.rds` | Cleaned SCE with unique cell IDs and `unused-*` / duplicate hybridization-cycle genes removed | **Yes** — read by `LC-NE_BARseq_MAPseq_analyses` |
| `DBHfilteredneurons_clust_CCFv2_uid.rds` | Same as above but filtered to putative LC-NE neurons (Dbh expression > 2) | No |

After a reproducible run from the released capsule, these two output folders are saved as separate AIND-metadata-tagged data assets in `aind-open-data`, with processing JSON pointing back to this capsule. Those published assets — not this capsule's `/results/` directly — are what the downstream analysis capsule mounts.

## How to run

Click **Reproducible Run** in Code Ocean. The `run` script processes both brains sequentially. Runtime is approximately 10 minutes on a large instance.

## Environment

R 4.2.3 in a conda environment (`r4-base`) with `hdf5r`, `Matrix`, and `SingleCellExperiment` as core dependencies. The full environment is defined in `environment/r4-base.yml`.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
```

### 7. Verify with a reproducible run

After all the above changes are committed and pushed to your branch:

1. From the new capsule, trigger a Reproducible Run on your branch.
2. When it completes, check `/results/`:
   - There should be exactly two top-level subfolders matching `barseq_<subject>_<date>_processed_MAT2RDS_<timestamp>/`.
   - Each subfolder should contain three `.rds` files plus an HTML report should also exist at the top level of `/results/`.
   - There should be NO `BARseq_669594/` folder.
3. Open one of the HTML reports and skim for errors. The script's `cat()` outputs should show duplicate-gene handling messages and the final saveRDS confirmations.

If any of those checks fail, fix the issue, commit, and re-run before opening the PR.

### 8. Open a PR against master

Push your branch to GitHub and open a pull request titled something like:

> Refactor: per-subject outputs, new bar-seq-only inputs, remove brain1

The PR description should call out the substantive changes:

- Removed brain1 (specimen 669594) — not used in manuscript or any downstream capsule
- Swapped to new bar-seq-only input data assets
- Output folders now follow `<input>_processed_MAT2RDS_<timestamp>` naming so each output corresponds 1:1 to one input asset
- Fixed brain3's latent crash bug (missing `length(rows_to_drop) > 0` guard)
- Removed `.Trash-0/` directory and added to `.gitignore`
- Rewrote README

### 9. Merge and release

Once the PR is approved and merged into master:

1. In the original Code Ocean capsule (`https://codeocean.allenneuraldynamics.org/capsule/3953531/tree`), pull the latest master.
2. Create a new release of that capsule. This is the version the manuscript collection links to, and it's the version whose reproducible run will produce the published derived data assets.

### 10. After release: produce the published derived assets

This step happens AFTER release, not as part of the PR:

1. From the released capsule, click Reproducible Run once.
2. From the resulting `/results/`, save each `barseq_<subject>_..._processed_MAT2RDS_<timestamp>/` subfolder as a separate data asset (Code Ocean's "save as data asset" lets you scope to a subdirectory).
3. Hand off those two assets to the SciComp colleague responsible for migrating to `aind-open-data` with proper AIND metadata and a processing JSON pointing at this capsule's released commit hash.
4. Once the published derived assets exist in `aind-open-data`, a follow-up PR on the downstream `LC-NE_BARseq_MAPseq_analyses` capsule will swap its `datasets.json` to mount those new per-subject assets in place of the current combined `BARseq_MATtoRDSfiles_brain3_brain4` asset. That follow-up PR is out of scope for this task.

## Notes for the agent

- All file paths in this briefing are relative to the repo root.
- If anything is ambiguous or doesn't match what you see in the repo (e.g., line numbers in the brain3 script have shifted), ask the human before guessing.
- Don't rename `01_BarSeq_RDSconvert_brain3_v2.R` and `_brain4_v2.R` to drop the `_v2` suffix — keeping the existing filenames preserves git history and avoids confusion against the released capsule.
- The `run` script's `render_html` helper invokes `rmarkdown::render` which writes the HTML to `/results/` at the top level. That stays unchanged.
- The new output folder names will include a timestamp, so two runs of the same capsule produce different folder names. This is intentional — it makes provenance traceable for whoever saves the final published assets.
