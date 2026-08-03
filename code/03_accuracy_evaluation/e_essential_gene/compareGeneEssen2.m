% Refer to Human-GEM： Human1_Publication_Data_Scripts\tINIT_GEMs
% SOURCE FILE NAME: compare_gene_essentiality_results.m, run_gene_essentiality_analysis.m
%The code compares the essential genes predicted by the model with experimental data from CRISPR screens.

cd ..\GitHub\cobratoolbox;
initCobraToolbox(false);
changeCobraSolver ('gurobi', 'all');
setRavenSolver('gurobi')

%% Load experimental information
cd ..\scGEM_benchmark\data\reference_informantion\essential_genes

load("data_model/expEssential_HCC827.mat"); % 1813
load("data_model/expNonEssential_HCC827.mat"); % 16592 (18333-1813)

%% Load and organize model-predicted gene essentiality results
taskStruct = parseTaskList('..\scGEM_benchmark\data\reference_data\metabolic_tasks\metabolicTasks_Essential.txt');
load('data_model/Human-GEM.mat');
refModel = ihuman;

% Convert refModel genes from Ensembl IDs to gene symbols
[grRules,genes,rxnGeneMat] = translateGrRules(refModel.grRules,'Name','ENSG');
refModel.grRules = grRules;
refModel.genes = genes;
refModel.rxnGeneMat = rxnGeneMat;

% initialize eGenes structure
eGenes = struct();
eGenes.refModel = refModel;
eGenes.taskList = {taskStruct(:).description}';
eGenes.tissues = {};
eGenes.geneList = {};
eGenes.essentialGenes = {};

% Load scGEMs
cd ..\scGEM_benchmark\output_all\LUAD\models;
load('models_ftINIT_LS_25th.mat');

% Select HCC827 models
% Load the cell-to-cell-line mapping information
cellMapping = readtable('cell_cellline_LUAD.xlsx');
% Select cell IDs belonging to the HCC827 cell line
hcc827CellIDs = cellMapping{strcmp(cellMapping{:, 2}, 'HCC827'), 1};
% Select HCC827 models from models
newModels = {};
for i = 1:length(models)
    cellID = models{i}.id;
    if any(strcmp(cellID, hcc827CellIDs))
        newModels{end + 1} = models{i};
    end
end
% newModels=newModels';
% % Select cell IDs
modelIds = cell(length(newModels), 1);
for i = 1:length(newModels)
    modelIds{i} = newModels{i, 1}.id;
end
models =newModels;

% % For rFASTCORMICS models, remove the rules field from each model
% for i = 1:length(models)
%     models{i} = rmfield(models{i},'rules');
% end


%%
tissues = cell(length(models),1);
for i = 1:length(models)
    tissues{i} = models{i}.id;
end

INIT_output.tissues = tissues;
INIT_output.model = models;

% Close exchange reactions
for i=1:length(models)
    id = models{i}.id;
    model = models{i};
    [~, selUpt] = findExcRxns(model);
    selUpt_rxns = model.rxns(selUpt);
    model = changeRxnBounds(model, selUpt_rxns, 0, 'l');
    INIT_output.id{i,1} = id;
    INIT_output.model{i,1} = model;
end

tic
eG = getTaskEssentialGenes(INIT_output, refModel, taskStruct);
toc
eGenes.tissues = [eGenes.tissues; eG.tissues];
eGenes.geneList = [eGenes.geneList; eG.geneList];
eGenes.essentialGenes = [eGenes.essentialGenes; eG.essentialGenes];


% Add the .essentialMat field to the eGenes structure
eGenes.essentialMat = false(numel(eGenes.refModel.genes), numel(eGenes.taskList), numel(eGenes.tissues)); %2864x57x10 logical

for i = 1:numel(eGenes.tissues) % 1:10
    [hasMatch,ind] = ismember(eGenes.refModel.genes,eGenes.geneList{i}); % hasMatch: 2864x1 logical;  ind:2864x1 double
    eGenes.essentialMat(hasMatch,:,i) = eGenes.essentialGenes{i}(ind(hasMatch),:); %2864x57x10 logical
end


%% Specify model-determined essential genes
modelPred = squeeze(any(eGenes.essentialMat,2));  % essential for any task

% Compare model predictions with experimental essentiality data
% calculate true and false positives and negatives
[TP,TN,FP,FN,Penr] = deal(nan(numel(tissues),1));  % initialize variables

for i = 1:numel(tissues) % 1:10

    modelGenes = eGenes.refModel.genes;
    if strcmpi(tissues{i},'all')
        modelEssential = modelGenes(sum(modelPred,2) == size(modelPred,2));
    else
        modelEssential = modelGenes(modelPred(:,i));
    end
    modelNonEssential = setdiff(modelGenes,modelEssential);


    TP(i) = sum(ismember(modelEssential, expEssential));        % true positives
    TN(i) = sum(ismember(modelNonEssential, expNonEssential));  % true negatives
    FP(i) = sum(ismember(modelEssential, expNonEssential));     % false positives
    FN(i) = sum(ismember(modelNonEssential, expEssential));     % false negatives

    % Penr(i) = EnrichmentTest(intersect(modelGenes,expGenes), intersect(modelEssential,expGenes), intersect(expEssential,modelGenes));
end

% calculate some metrics
sensitivity = TP./(TP + FN);
specificity = TN./(TN + FP);
accuracy = (TP + TN)./(TP + TN + FP + FN);
F1 = 2*TP./(2*TP + FP + FN);
MCC = ((TP.*TN) - (FP.*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));  % Matthews correlation coefficient

% Save model-predicted essential genes
for i = 1:numel(tissues) % 1:10

    modelGenes = eGenes.refModel.genes;
    modelEssential2{i,:} = modelGenes(modelPred(:,i));

end

save('../EssentialGenes/HCC827/GeneEssen_ftINIT_LS_25th.mat',"modelEssential2")