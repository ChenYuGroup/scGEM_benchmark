# scGEM_benchmark

## Introduction

The `scGEM_benchmark` repository contains the data, code and processed results supporting the study “A benchmarking framework for single-cell genome-scale metabolic model construction.”

## Code

The `./code/` directory contains all scripts used for the analyses reported in this study. These scripts include:

* preprocessing of scRNA-seq data;
* construction of scGEMs using different strategies;
* evaluation of scGEM accuracy;
* sensitivity analysis based on expression perturbation;
* evaluation of computational feasibility;
* analysis of dataset dependence and cross-dataset consistency;
* analysis of model size and its effects on accuracy; and
* statistical analysis of benchmark results.

## Data

The `./data/` directory contains scRNA-seq expression matrices, cell labels, model input files and reference data used in the benchmark. 

## Result

The `./result/` directory contains the processed results generated in this study, including:

* `accuracy_performance.xlsx`: accuracy scores of the 26 construction strategies across datasets;
* `clustering_performance_scRNA.xlsx`: clustering performance calculated from the scRNA-seq datasets;
* `feasibility_performance.xlsx`: model construction runtime, timeout events and success status;
* `Kendall_W_ranking_agreement.xlsx`: cross-dataset agreement in strategy rankings assessed using Kendall’s W;
* `sensitivity_performance.xlsx`: structural and accuracy responses to expression perturbation; and
* `variance_decomposition.xlsx`: contributions of construction factors to variation in model performance.

## Figure

The `./figure/` directory contains the R scripts used to generate the figures in this study. Each script is named according to the corresponding figure panel.

## Requirements

The main software and dependencies include:

* MATLAB R2023a (version 9.14);
* R (version 4.3.2);
* COBRA Toolbox (version 3.6);
* RAVEN Toolbox (version 2.8.4);
* Human-GEM (Human1 version 1.16);
* Gurobi Optimizer 13.0.1; 
* IBM ILOG CPLEX 12.10.





