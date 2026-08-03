%% Model size: numbers of reactions, genes and metabolites
% load scGEMs
cd ..\scGEM_benchmark\output_all\LUAD\models;
load('models_ftINIT_LS_25th.mat');

% Numbers of reactions, genes and metabolites
rxn_counts = [];
gene_counts = [];
met_counts = [];
skip_count = 0;

for i = 1:numel(models)
    model = models{i};
    
    % Check whether the model is valid and contains the rxns, genes and mets fields
    if isstruct(model) && ...
            isfield(model, 'rxns') && ...
            isfield(model, 'genes') && ...
            isfield(model, 'mets') && ...
            ~isempty(model.rxns)

        rxn_counts(end+1,1) = numel(model.rxns);
        gene_counts(end+1,1) = numel(model.genes);
        met_counts(end+1,1) = numel(model.mets);
    else
        skip_count = skip_count + 1;
    end
end

% Calculate the mean values
mean_rxn_num = mean(rxn_counts);
mean_gene_num = mean(gene_counts);
mean_met_num = mean(met_counts);

result = [mean_rxn_num, mean_gene_num, mean_met_num];

% Display the results
fprintf('Average reactions per valid model: %.2f\n', mean_rxn_num);
fprintf('Average genes per valid model: %.2f\n', mean_gene_num);
fprintf('Average metabolites per valid model: %.2f\n', mean_met_num);
fprintf('Skipped models: %d\n', skip_count);