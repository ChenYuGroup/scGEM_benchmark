# Model construction

This directory contains MATLAB scripts for constructing single-cell genome-scale metabolic models (scGEMs) using GIMME, iMAT, ftINIT and rFASTCORMICS.

The scripts use `LS_Selected_LUAD.xlsx` as an example. For other datasets or preprocessing settings, update the input file, applicable expression thresholds, output directory and output filenames consistently. Local paths to the required toolboxes and data files must also be configured before running the scripts.

## `generate_GIMME_model.m`

Constructs one GIMME-based scGEM for each column of the input expression matrix.

The [`get_GPR.m`](https://github.com/hongzhonglu/yeast_GEM_multi_omics_analysis/blob/main/single_cell/code/get_GPR.m) function from the [yeast_GEM_multi_omics_analysis](https://github.com/hongzhonglu/yeast_GEM_multi_omics_analysis) repository must be added to the MATLAB path. The remaining model-construction functions are provided by the COBRA Toolbox.

### Input

- `Human-GEM.mat`: Human1 v1.16 reference model;
- `LS_Selected_LUAD.xlsx`: quality-controlled LS-LUAD expression matrix, with genes in rows and cells in columns.

### Key settings

- Global expression preprocessing threshold: `0.039` (25th percentile);
- GIMME reaction-expression threshold: `0`;
- metabolic extraction method (MEM): GIMME;


### Output

- `models_GIMME_LS_25th.mat`: reconstructed models;
- `times_GIMME_LS_25th.mat`: construction time per model.

## `generate_iMAT_model.m`

Constructs one iMAT-based scGEM for each column of the input expression matrix.

The [`get_GPR.m`](https://github.com/hongzhonglu/yeast_GEM_multi_omics_analysis/blob/main/single_cell/code/get_GPR.m) function from the [yeast_GEM_multi_omics_analysis](https://github.com/hongzhonglu/yeast_GEM_multi_omics_analysis) repository must be added to the MATLAB path. The remaining model-construction functions are provided by the COBRA Toolbox.

### Input

- `Human-GEM.mat`: Human1 v1.16 reference model;
- `LS_Selected_LUAD.xlsx`: quality-controlled LS-LUAD expression matrix, with genes in rows and cells in columns.

### Key settings

- Global expression preprocessing threshold: `0.039` (25th percentile);
- iMAT lower and upper thresholds: `1e-10`;
- MEM: iMAT;
- optimization solver: Gurobi.

### Output

- `models_iMAT_LS_25th.mat`: reconstructed models;
- `times_iMAT_LS_25th.mat`: construction time per model.

## `generate_ftINIT_model.m`

Constructs one ftINIT-based scGEM for each column of the input expression matrix.

### External requirements and file locations

The external functions required by this script are not redistributed in this repository and should be obtained from the following sources:

- `ftINIT`: provided by the [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) at `INIT/ftINIT.m`;
- `getHumanGEMINITSteps`: provided by [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) at `code/tINIT/getHumanGEMINITSteps.m`;
- `prepHumanModelForftINIT`: provided by Human-GEM at `code/tINIT/prepHumanModelForftINIT.m`.

The supplied `prepDataHumanGEM.mat` file was generated using `MatlabCode/prepHumanModel.m` from the [SingleCellModeling](https://github.com/SysBioChalmers/SingleCellModeling) repository. Model preparation used Human1 v1.16 together with:

- `data/metabolicTasks/metabolicTasks_Essential.txt` from Human-GEM; 

The RAVEN Toolbox, Human-GEM and COBRA Toolbox repositories must be added to the MATLAB path, and Gurobi must be configured as the optimization solver. The preparation step does not need to be repeated when using the supplied `prepDataHumanGEM.mat`.

### Input

- `prepDataHumanGEM.mat`: ftINIT-compatible input generated from Human1 v1.16 using `prepHumanModel.m`;
- `LS_Selected_LUAD.xlsx`: quality-controlled LS-LUAD expression matrix, with genes in rows and cells in columns.

### Key settings

- Global expression threshold: `0.039` (25th percentile);
- INIT steps: `getHumanGEMINITSteps('1+0')`;
- MEM: ftINIT;
- optimization solver: Gurobi.

### Output

- `models_ftINIT_LS_25th.mat`: reconstructed models;
- `times_ftINIT_LS_25th.mat`: construction time per model.

## `generate_rFASTCORMICS_model.m`

Constructs one rFASTCORMICS-based scGEM for each column of the input expression matrix.

The script was adapted from `rFASTCORMICS for RNA-seq data/Scripts/fastcormics_RNAseq.m` in the [rFASTCORMICS](https://github.com/sysbiolux/rFASTCORMICS) repository. rFASTCORMICS must be downloaded and added to the MATLAB path. The COBRA Toolbox and IBM ILOG CPLEX Studio 12.10 are also required.

### Input

- `LS_Selected_LUAD.xlsx`: quality-controlled LS-LUAD expression matrix, with genes in rows and cells in columns;
- `Cmodel_rFASTCM.mat`: flux-consistent Human1 reference model generated using `fastcc_4_rfastcormics` with a tolerance of `1e-4`;
- `dico_rFASTCORMICS.mat`: dictionary used to map expression-matrix gene identifiers to genes in the reference model.

`Human-GEM.mat` (Human1 v1.16) is required when regenerating `Cmodel_rFASTCM.mat`.

### Key settings

- Expression values are rescaled column-wise to a total of $10^6$;
- expression discretization: `discretize_FPKM`;
- MEM: rFASTCORMICS;
- optimization solver: IBM ILOG CPLEX Studio 12.10.

### Output

- `models_rFAST_LS.mat`: reconstructed models;
- `times_rFAST_LS.mat`: construction time per model.
