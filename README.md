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
