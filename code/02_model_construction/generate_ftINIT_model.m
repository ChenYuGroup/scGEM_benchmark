%% External dependencies
% ftINIT is implemented in the RAVEN Toolbox:
% https://github.com/SysBioChalmers/RAVEN
%
% getHumanGEMINITSteps and prepHumanModelForftINIT are provided by
% the Human-GEM project: % https://github.com/SysBioChalmers/Human-GEM
%
% Reference model: Human1 version 1.16
% MILP solver: Gurobi 
%
% The external functions are not redistributed in this repository.
% Install the specified RAVEN and Human-GEM versions before running
% this script.

%% Runs prepInitModel for Human-GEM, including filtering of reactions.
cd ..\GitHub\cobratoolbox;
initCobraToolbox(false);
changeCobraSolver ('gurobi', 'all');

% prepDataHumanGEM was generated using prepHumanModel.m from:
% https://github.com/SysBioChalmers/SingleCellModeling
% File: MatlabCode/prepHumanModel.m
% The script was run using Human1 version 1.16 and the corresponding essential metabolic task and reaction
% annotation files.
load('data/data_model/prepDataHumanGEM.mat');

%-------------------- scRNA-Seq data --------------------
cd ..\scGEM_benchmark\data\expression_matrices
[num,txt] = xlsread('LS_Selected_LUAD.xlsx');

arrayData.genes = txt((2:end), 1);
tissues = txt(1, (2:end))';
arrayData.tissues = tissues(1:end,:);
arrayData.levels = num(:,1:end);

% Calculated as the percentile across all filtered genes and cells
arrayData.threshold = 0.039; % 25th 


%% run ftINIT
times = {};
models = cell(length(arrayData.tissues),1);

for i = 1:length(arrayData.tissues)

    N=i
    tic
    mres = ftINIT(prepDataHumanGEM,arrayData.tissues{i},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);

    times{end+1} = toc;

    mres.id = arrayData.tissues{i};
    models{i,1} = mres;

end
times = times';

              
% save results
cd ..\scGEM_benchmark\output_all\LUAD
save('models/models_ftINIT_LS_25th.mat',"models")
save('times/times_ftINIT_LS_25th.mat',"times")
