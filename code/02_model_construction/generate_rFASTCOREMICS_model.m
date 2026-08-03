%% External dependencies
% Optimization solver: IBM ILOG CPLEX Studio 12.10
%
% rFASTCORMICS source:
% https://github.com/sysbiolux/rFASTCORMICS
%
% This script was adapted from:
% rFASTCORMICS for RNA-seq data/Scripts/fastcormics_RNAseq.m

cd ..\GitHub\cobratoolbox;
initCobraToolbox(false);

addpath(genpath(pwd)) % add all subfolders to path in the current dicrectory. 
feature astheightlimit 2000

%-------------------- scRNA-Seq data --------------------
cd ..\scGEM_benchmark\data\expression_matrices
[num,txt] = xlsread('LS_Selected_LUAD.xlsx');

%Data discretization
colnames = txt(1, (2:end)); 
rownames= txt((2:end), 1); 
fpkm=num; 
fpkm = fpkm.*10^6./sum(fpkm,1); 
discretized = discretize_FPKM(fpkm, colnames); % no figures

%Prepare rFASTCORMICS
load('data_model/Human-GEM.mat');
model = ihuman;

% Generate a flux-consistent reference model using FASTCC (Vlassis et al., 2014). The function is provided by rFASTCORMICS:
% ..\GitHub\rFASTCORMICS\rFASTCORMICS for RNA-seq data\Scripts\fastcc_4_rfastcormics.m
% A = fastcc_4_rfastcormics(model, 1e-4,0); % create consistent model by running FASTCC (Vlassis et al., 2014)
% Cmodel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
load('data_model/Cmodel_rFASTCM.mat');

load dico_rFASTCORMICS.mat % dictionary to map the rownname identifier to the genes in the model

% Create models
% single models:
models = cell(size(fpkm,2),1);
times = {};


for i = 1:numel(colnames)

    tic
    N=i
    [ContextModel, A_keep] = fastcormics_RNAseq(Cmodel, discretized(:,i), rownames, dico);

    times{end+1} = toc;
    
    ContextModel.id = colnames{i};
    models{i,1} =ContextModel;
end

times = times';


% save results
cd ..\scGEM_benchmark\output_all\LUAD
save('models/models_rFAST_LS.mat',"models")
save('times/times_rFAST_LS.mat',"times")
