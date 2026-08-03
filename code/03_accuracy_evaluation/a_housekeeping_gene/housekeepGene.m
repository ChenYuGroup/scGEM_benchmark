%% Extract genes from the grRules field as model genes
% load scGEMs
cd ..\scGEM_benchmark\output_all\LUAD\models;
load('models_ftINIT_LS_25th.mat');

% Initialize a cell array to store the gene list for each model
model_genes = cell(size(models));
% Iterate over each model
for i = 1:numel(models)
    % Handle empty models
    if isempty(models{i})
        model_genes{i} = {};
        continue;
    end

    % Obtain the grRules field of the current model
    grRules = models{i, 1}.grRules;
    all_genes = {};
    for j = 1:numel(grRules)
        rule = grRules{j}; % Current rule (string)
        % Use a regular expression to extract ENSG gene IDs from the rule
        genes_in_rule = regexp(rule, 'ENSG[0-9]+', 'match');
        % Add the extracted genes to the gene collection
        all_genes = [all_genes; genes_in_rule(:)];
    end
    % Remove duplicates and store the gene list for the current model
    model_genes{i} = unique(all_genes);
end
% Save the results
cd ..\scGEM_benchmark\output_all\LUAD
save('genes/genes_ftINIT_LS_25th.mat',"model_genes")

%% Calculate housekeeping gene coverage in scGEMs
% Check for and remove empty cells
emptyIdx = cellfun(@isempty, model_genes);
fprintf('Number of empty cells = %d\n', sum(emptyIdx));
model_genes(emptyIdx) = [];
fprintf('Number of cells after removal = %d\n', numel(model_genes));

%% 1. sc2019
cd ..\scGEM_benchmark\data\reference_data\housekeeping_genes;
load('gene_HkMe_sc2019.mat');
% Initialize the result arrays
numCommonGenes = zeros(length(model_genes), 1);
ratio1 = zeros(length(model_genes), 1);
% Iterate over each model
for i = 1:length(model_genes)

    numCommonGenes(i) = sum(ismember(model_genes{i}, gene_HkMe));
    ratio1 = numCommonGenes/numel(gene_HkMe)*100;
end
meanratio1 = mean(ratio1);

%% 2. sc2022
cd ..\scGEM_benchmark\data\reference_data\housekeeping_genes;
load('gene_HkMe_sc2022.mat');
% Initialize the result arrays
numCommonGenes = zeros(length(model_genes), 1);
ratio2 = zeros(length(model_genes), 1);
% Iterate over each model
for i = 1:length(model_genes)
    % Use ismember to determine the number of shared genes
    numCommonGenes(i) = sum(ismember(model_genes{i}, gene_HkMe));
    ratio2 = numCommonGenes/numel(gene_HkMe)*100;
end
meanratio2 = mean(ratio2);

%% 3. sc2025
cd ..\scGEM_benchmark\data\reference_data\housekeeping_genes;
load('gene_HkMe_sc2025.mat');
% Initialize the result arrays
numCommonGenes = zeros(length(model_genes), 1);
ratio3 = zeros(length(model_genes), 1);
% Iterate over each model
for i = 1:length(model_genes)
    % Use ismember to determine the number of shared genes
    numCommonGenes(i) = sum(ismember(model_genes{i}, gene_HkMe));
    ratio3 = numCommonGenes/numel(gene_HkMe)*100;
end
meanratio3 = mean(ratio3);
result=[meanratio1,meanratio2,meanratio3];