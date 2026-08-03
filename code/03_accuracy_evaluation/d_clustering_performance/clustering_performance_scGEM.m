% load scGEMs
cd ..\scGEM_benchmark\output_all\LUAD\models;
load('models_ftINIT_LS_25th.mat');

% Set the random seed to ensure reproducibility
rng(0); % Fix the random seed

% Step 1: Extract reaction sets from all models and construct the feature matrix
allReactions = strings(0); % Initialize an empty string array
for i = 1:length(models)
    allReactions = union(allReactions, string(models{i}.rxns)); % Convert to a string array and merge
end

numModels = length(models);
numReactions = length(allReactions);
featureMatrix = false(numModels, numReactions);
for i = 1:numModels
    modelRxns = string(models{i}.rxns); % Ensure that reactions are represented as a string array
    for j = 1:numReactions
        if ismember(allReactions(j), modelRxns) % Check whether the current reaction is present in the model
            featureMatrix(i, j) = 1; % Mark as 1
        end
    end
end

% Step 2: Load cell-line information from the Excel table
cellLineData = readtable('cell_cellline_LUAD.xlsx');
cellLineMapping = containers.Map(cellLineData.cell, cellLineData.cellline);

% Obtain the cell-line array and identify unmatched models
cellLines = cell(numModels, 1);
unmatchedModels = {};

for i = 1:numModels
    cellID = models{i}.id;
    if isKey(cellLineMapping, cellID) % Check whether a matching entry exists
        cellLines{i} = cellLineMapping(cellID); % Obtain the corresponding cell line
    else
        cellLines{i} = 'Unknown'; % Label unmatched models as 'Unknown'
        unmatchedModels{end+1} = cellID; % Record unmatched models
    end
end

if ~isempty(unmatchedModels)
    disp('Warning: The following models were not matched to any cell line:');
    disp(unmatchedModels);
end

% Prepare valid data
validIndices = ~strcmp(cellLines, 'Unknown');
validCellLines = cellLines(validIndices);
uniqueCellLines = unique(validCellLines);
numClusters = length(uniqueCellLines);
if numClusters < 1
    warning('No valid cell lines were detected. The number of clusters was set to 1.');
    numClusters = 1;
end
disp(['Number of clusters: ', num2str(numClusters)]);

% Step 3: Prepare the true labels
trueLabels = zeros(size(validCellLines));
for i = 1:length(validCellLines)
    trueLabels(i) = find(strcmp(uniqueCellLines, validCellLines{i}));
end

%%
% Step 4: Perform hierarchical clustering
% Convert featureMatrix to double precision for use with the linkage function
featureMatrixDouble = double(featureMatrix(validIndices, :));

Z = linkage(featureMatrixDouble, 'ward', 'jaccard'); % Use Jaccard distance and Ward's method

clusterLabels = cluster(Z, 'maxclust', numClusters); % Cut the hierarchy according to the expected number of clusters

%% Calculate the ARI (Adjusted Rand Index)
ari = rand_index(trueLabels, clusterLabels);
fprintf('Adjusted Rand Index (ARI): %.4f\n', ari);

% Calculate the NMI (Normalized Mutual Information)
NMI = nmi(trueLabels, clusterLabels);
fprintf('Normalized Mutual Information (NMI): %.4f\n', NMI);

% Calculate the silhouette coefficient
silhouetteCoefficients = silhouette(featureMatrixDouble, clusterLabels, 'jaccard');
averageSilhouette = mean(silhouetteCoefficients);
fprintf('Silhouette Coefficient: %.4f\n', averageSilhouette);

% Calculate the Calinski-Harabasz Index
chObj = evalclusters(featureMatrixDouble, clusterLabels, 'CalinskiHarabasz');
chIndex = chObj.CriterionValues;
fprintf('Calinski - Harabasz Index: %.4f\n', chIndex);

% Calculate the Davies-Bouldin Index
dbObj = evalclusters(featureMatrixDouble, clusterLabels, 'DaviesBouldin');
dbIndex = dbObj.CriterionValues;
fprintf('Davies - Bouldin Index: %.4f\n', dbIndex);

% Calculate compactness
n = size(featureMatrixDouble, 2); % Number of features
compactness = CpS(clusterLabels, featureMatrixDouble, n);
fprintf('Compactness: %.4f\n', compactness);

% Calculate the mean intra-cluster distance and within-cluster variance
totalIntraDistance = 0;
totalWithinVariance = 0;

for k = 1:numClusters
    clusterMembers = featureMatrixDouble(clusterLabels == k, :);
    numMembers = size(clusterMembers, 1);

    % Calculate the mean intra-cluster distance
    if numMembers > 1
        distances = pdist(clusterMembers);
        intraDistance = mean(distances);
        totalIntraDistance = totalIntraDistance + intraDistance;
    end

    % Calculate the within-cluster variance
    clusterCenter = mean(clusterMembers);
    squaredDistances = sum((clusterMembers - repmat(clusterCenter, numMembers, 1)).^2, 2);
    withinVariance = sum(squaredDistances) / numMembers;
    totalWithinVariance = totalWithinVariance + withinVariance;
end

% Calculate the mean intra-cluster distance and mean within-cluster variance across all clusters
averageIntraDistance = totalIntraDistance / numClusters;
averageWithinVariance = totalWithinVariance / numClusters;

fprintf('Mean intra-cluster distance: %.4f\n', averageIntraDistance);
fprintf('Within-cluster variance: %.4f\n', averageWithinVariance);

%% Save the metrics
% Initialize the index variable
index = cell(2, 8);

% Store metric names
index{1, 1} = 'ARI';
index{1, 2} = 'NMI';
index{1, 3} = 'Silhouette';
index{1, 4} = 'chIndex';
index{1, 5} = 'dBIndex';
index{1, 6} = 'Intra-cluster Distance';
index{1, 7} = 'Within-cluster Variance';
index{1, 8} = 'Compactness';

% Store metric values
index{2, 1} = ari; % Adjusted Rand Index
index{2, 2} = NMI; % Normalized mutual information
index{2, 3} = averageSilhouette; % Mean silhouette coefficient
index{2, 4} = chIndex; % Calinski-Harabasz Index
index{2, 5} = dbIndex; % Davies-Bouldin Index
index{2, 6} = averageIntraDistance; % Mean intra-cluster distance
index{2, 7} = averageWithinVariance; % Within-cluster variance
index{2, 8} = compactness; % Compactness

disp('Metric results:');
disp(index);