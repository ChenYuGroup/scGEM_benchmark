# Code

This directory contains the MATLAB and R scripts used for data preprocessing, scGEM construction and benchmark analysis. The numbered subdirectories represent the main stages of the workflow and should generally be followed in numerical order.

## Workflow

### `01_data_preprocessing/`

Preprocesses the scRNA-seq expression matrices, calculates expression thresholds and selects cells for large-scale model construction. See [`01_data_preprocessing/README.md`](01_data_preprocessing/README.md) for details.

### `02_model_construction/`

Constructs single-cell genome-scale metabolic models using different combinations of expression preprocessing, threshold settings and metabolic extraction methods. See [`02_model_construction/README.md`](02_model_construction/README.md) for details.

### `03_accuracy_evaluation/`

Evaluates scGEM accuracy using five complementary criteria:

- housekeeping-gene coverage;
- essential metabolic-task completion;
- marker-gene consistency;
- cell-type clustering performance; and
- essential-gene prediction.


### `04_sensitivity_analysis/`

Contains `generate_perturbed_expression_matrices.m`, which generates perturbed expression matrices from the selected LS-LUAD matrix at predefined Spearman correlation levels. The generated matrices are stored under `./data/expression_matrices/LUAD/01_LUAD_sensitivity/`.

The sensitivity-analysis workflow consists of three stages:

1. generate the perturbed expression matrices using `generate_perturbed_expression_matrices.m`;
2. construct scGEMs from each perturbed matrix using the scripts in [`02_model_construction/`](02_model_construction/); and
3. evaluate the reconstructed models using the scripts in [`03_accuracy_evaluation/`](03_accuracy_evaluation/).

### `05_feasibility_analysis/`

Summarizes computational feasibility, including model-construction runtime, timeout events and construction success.

### `06_dataset_dependence_analysis/`

Assesses dataset dependence and cross-dataset agreement in the performance rankings of model-construction strategies.

### `07_model_size_analysis/`

Calculates model size and evaluates its association with scGEM accuracy.

### `08_statistical_analysis/`

Performs statistical comparisons and variance decomposition of the benchmark results.

## General instructions

1. Install and configure the software and dependencies listed in the [Requirements](../README.md#requirements) section of the main README.
2. Download the expression matrices from [Zenodo](https://doi.org/10.5281/zenodo.21774476), extract the archives and place the files under `./data/expression_matrices/`.
3. Before running a script, replace any machine-specific paths with the corresponding paths on the local computer.
4. Adjust dataset names, expression thresholds, input files and output filenames consistently when applying the scripts to other datasets or preprocessing settings.
5. Run the scripts in numerical order unless an analysis uses previously generated intermediate results.

