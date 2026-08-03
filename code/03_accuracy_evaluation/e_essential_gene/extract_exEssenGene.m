% Refer to Human-GEM: Human1_Publication_Data_Scripts\tINIT_GEMs
% SOURCE FILE NAME: compare_gene_essentiality_results.m

load('data_model/Human-GEM.mat');
refModel = ihuman;

[grRules,genes,rxnGeneMat] = translateGrRules(refModel.grRules,'Name','ENSG');
refModel.grRules = grRules;
refModel.genes = genes;
refModel.rxnGeneMat = rxnGeneMat;

%% Experimental essential gene data
x = readtable('data/Achilles_gene_effect_22Q2.csv');
genes = x.Properties.VariableDescriptions(2:end)'; %depmap
genes = regexprep(genes, '\(.+\)', '');
genes = regexprep(genes, ' ', '');

celltypes = x.DepMap_ID(1:end);


col_types = table2cell(varfun(@class, x))';
cell_col = x.Properties.VariableNames(ismember(col_types,'cell'))';
cell_col = setdiff(cell_col,'Var1');  % skip the gene column
for i = 1:numel(cell_col)
    x.(cell_col{i})(ismember(x.(cell_col{i}), 'NA')) = {'NaN'};
    x.(cell_col{i}) = cellfun(@str2double, x.(cell_col{i}));
end

% convert to matrix and transpose
x = table2array(x(:,2:end))';

% collect results
expdata = {};
expdata.genes = genes;
expdata.tissues = celltypes;
expdata.gene_effect = x;


%% % specify essentiality threshold (for DepMap data only)
DepMapThresh = -0.6;

expdata.essential = double(expdata.gene_effect < DepMapThresh);

% Load
tissue_ind = find(strcmp(expdata.tissues, 'ACH-000551'));

expGenes = expdata.genes(~isnan(expdata.essential(:,tissue_ind)));
expEssential = expdata.genes(expdata.essential(:,tissue_ind) == 1);
expNonEssential = setdiff(expGenes,expEssential);

save('../data_model/expNonEssential_A549.mat', 'expNonEssential');
save('../data_model/expEssential_A549.mat', 'expEssential');