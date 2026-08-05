# Data preprocessing

This directory contains the scripts used to preprocess the scRNA-seq expression matrices, calculate expression thresholds and select cells for large-scale model construction. Run the scripts in the order described below.

The preprocessing script requires R with the `SingleCellExperiment`, `Linnorm` and `SAVER` packages. MATLAB is required for the threshold-calculation and cell-sampling scripts.

## 1. `a_data_preprocessing.r`

Performs quality control, normalization and imputation of the original expression matrices. The source datasets can be downloaded from GEO or Figshare using the accession numbers and links listed in the [Datasets](../../README.md#datasets) section of the main README.

The preprocessing workflow produces:

- `R_*.xlsx`: quality-controlled expression matrices before normalization and imputation;
- `LS_*.xlsx`: matrices normalized using Linnorm and subsequently imputed using SAVER.

The processed matrices are stored in `./data/expression_matrices/`, for example `R_LUAD.xlsx` and `LS_LUAD.xlsx`.

## 2. `b_run_expression_threshold_calculation.m`

Calculates the 25th, 50th, 75th and 90th percentile expression thresholds across the filtered genes and cells. These thresholds are used in the subsequent model-construction analyses.

## 3. `c_select_randomSample.m`

Randomly selects cells for large-scale model construction and saves the sampling order or indices to support reproducibility. Sampling indices are stored in `./data/cell_labels/`, for example `randomIndices_LUAD.mat`. The selected expression matrices are stored in `./data/expression_matrices/`, for example `R_Selected_LUAD.xlsx` and `LS_Selected_LUAD.xlsx`.

The PAAD and BRCA datasets retain all available cells; therefore, separate `Selected` expression matrices are not generated for these datasets.

## Output data

The quality-controlled complete and selected expression matrices used in the benchmark are available from [Zenodo](https://doi.org/10.5281/zenodo.21774476). After extraction, place them under `./data/expression_matrices/` in the repository.
