% Calculate marker gene coverage
cd scGEM_benchmark\data\reference_data\marker_genes
load('markerGenes_LUAD.mat');

% Convert marker genes to the required identifier format
data = readtable('GeneMapping.tsv', 'Delimiter', '\t', 'FileType', 'text'); %

marker = string(marker_genes);
[tf, idx] = ismember(upper(marker), upper(string(data.Symbol)));
% 1. Display unmatched gene names
unmatched_genes = marker(~tf);
disp('Unmatched marker genes:')
disp(unmatched_genes)
fprintf('Number of unmatched genes = %d\n', length(unmatched_genes));

% 2. Remove unmatched genes
marker_genes = marker_genes(tf);
marker_genes_ENSG = data.ENSG(idx(tf));

% 3. Display the number of successfully matched genes
fprintf('Number of successfully matched marker genes = %d\n', length(marker_genes_ENSG));



%% % Load model gene lists
cd ..\scGEM_benchmark\output_all\LUAD\genes
load('genes_ftINIT_LS_25th.mat');

% % Extract models by category if required, for example, for the ESC dataset
% cellMapping = readtable('cell_cellline_ESC.xlsx');
% load('randomIndices_ESC.mat') % Positions of the 100 selected models among all cells
% % Obtain the cell-line labels corresponding to the 100 selected cells
% selected_cellline = cellMapping.cellline(randomIndices);
% % Identify the positions of models belonging to H1_ESC and DEC
% idx_H1_ESC = strcmp(selected_cellline, 'H1_ESC');
% idx_DEC    = strcmp(selected_cellline, 'DEC');
% % Extract gene lists by cell line
% model_genes_H1_ESC = model_genes(idx_H1_ESC);
% model_genes_DEC    = model_genes(idx_DEC);


%% Calculate the metric
n_model  = length(model_genes);
n_marker = length(marker_genes_ENSG);
result = zeros(n_model, n_marker);

for i = 1:n_model    
    genes_i = model_genes{i};   % Gene list of the ith model    
    result(i,:) = ismember(marker_genes_ENSG, genes_i)';    
end

% Calculate marker gene coverage for each model
ratio = sum(result, 2) / size(result, 2);
mean_ratio = mean(ratio);
fprintf('Mean marker gene coverage across models = %.4f\n', mean_ratio);
fprintf('Mean marker gene coverage across models = %.2f%%\n', mean_ratio * 100);