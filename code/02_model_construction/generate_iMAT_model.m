cd ..\GitHub\cobratoolbox;
initCobraToolbox(false);
changeCobraSolver ('gurobi', 'all');

load('data_model/Human-GEM.mat');
model = ihuman;

%-------------------- scRNA-Seq data --------------------
cd ..\scGEM_benchmark\data\expression_matrices
[num,txt] = xlsread('LS_Selected_LUAD.xlsx');


globalth = 0.039; % 25th
num(num <= globalth) = 0;

name = txt(1, (2:end));
exp.gene = txt((2:end), 1);

nModels = size(num,2);
model_indx = 1:nModels;
models = cell(nModels,1);
times = {};

for i = 1:nModels

    N=i
    exp.value = num(:,i);

    % get_GPR was obtained from:
    % https://github.com/hongzhonglu/yeast_GEM_multi_omics_analysis/blob/main/single_cell/code/get_GPR.m
    parsedGPR = get_GPR(model);
    [gene_id, gene_expr] = findUsedGenesLevels(model,exp);
    expressionRxns = selectGeneFromGPR(model, gene_id, gene_expr, parsedGPR, false);
    expressionRxns(isnan(expressionRxns)) = -1;

    options.expressionRxns = expressionRxns;

    options.threshold_lb = 1e-10;
    options.threshold_ub = 1e-10;

    options.solver = 'iMAT';

    tic
    model_p = createTissueSpecificModel(model, options);
    toc

    % [model_p, unusedGenes] = removeUnusedGenes(model_p);

    times{end+1} = toc;

    model_p.id = name{i};
    models{i,1} =model_p;

end

times=times';

% save results
cd ..\scGEM_benchmark\output_all\LUAD
save('models/models_iMAT_LS_25th.mat',"models")
save('times/times_iMAT_LS_25th.mat',"times")
