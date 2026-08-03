%% Perturbation-level experiment
% Objective: Generate perturbed datasets with Spearman correlations of
% approximately 0, 0.1, 0.3, 0.5, 0.7 and 0.9 with the original expression matrix
% Perturbation method: Randomly permute the gene expression profile of each cell
% and apply linear interpolation
% Use a grid search to determine lambda and avoid fzero errors

clear; clc;

%% 1. Load data
cd ..\scGEM_benchmark\data\expression_matrices\01_LUAD_sensitivity

[num, txt] = xlsread('LS_Selected_LUAD.xlsx');

X = num;
if any(isnan(X(:)))
    error('The expression matrix contains NaN values. Please check the input file.');
end

[nGenes, nCells] = size(X);
fprintf('Expression matrix size: %d genes × %d cells\n', nGenes, nCells);

%% 2. Set target correlations
target_rhos = [0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7, 0.8, 0.9];

% Lambda search range
lambda_grid = linspace(0, 1, 1001);

% Set a fixed random seed to ensure reproducibility
rng(2026);

%% 3. Create the output directory
outDir = fullfile(pwd, 'Perturbation_level_data');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end


%% 5. Generate perturbed datasets at different target correlations
for t = 1:length(target_rhos)

    target_rho = target_rhos(t);
    result = zeros(size(X));
    best_lambda = zeros(1, nCells);
    rho_check = zeros(1, nCells);

    fprintf('\nGenerating dataset with target Spearman rho = %.1f\n', target_rho);

    %% Perturb each cell independently
    for i = 1:nCells

        x = X(:, i);

        % If the expression values of the cell are constant, the correlation
        % cannot be calculated; retain the original column
        if std(x) == 0
            result(:, i) = x;
            best_lambda(i) = 0;
            rho_check(i) = NaN;
            continue;
        end

        % Randomly permute the gene expression values of the cell
        r = x(randperm(length(x)));

        rho_grid = nan(size(lambda_grid));

        % Perform a grid search over lambda
        for j = 1:length(lambda_grid)

            lambda = lambda_grid(j);
            y = x + lambda * (r - x);

            rho_grid(j) = corr(x, y, ...
                'Type', 'Spearman', ...
                'Rows', 'complete');
        end

        % Select the lambda that produces the correlation closest to the target
        [~, bestIdx] = min(abs(rho_grid - target_rho));

        lambda_best = lambda_grid(bestIdx);
        y_best = x + lambda_best * (r - x);

        result(:, i) = y_best;
        best_lambda(i) = lambda_best;
        rho_check(i) = rho_grid(bestIdx);

    end

    %% 6. Check the achieved correlations
    mean_rho = mean(rho_check, 'omitnan');
    min_rho  = min(rho_check, [], 'omitnan');
    max_rho  = max(rho_check, [], 'omitnan');

    fprintf('Mean Spearman rho = %.4f\n', mean_rho);
    fprintf('Range = %.4f ~ %.4f\n', min_rho, max_rho);

    %% 7. Save the Excel file while retaining the original gene and cell names
    outCell = txt;
    outCell(2:end, 2:end) = num2cell(result);

    rhoName = strrep(num2str(target_rho), '.', 'p');
    outFile = fullfile(outDir, ['Perturb_rho_', rhoName, '.xlsx']);

    writecell(outCell, outFile);


end