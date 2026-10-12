%% step1: Calculate the ARI and other metrics using metabolic genes
% Set the random seed to ensure reproducibility
rng(0); % Fix the random seed

% Load the single-cell transcriptomic data
[num, txt] = xlsread('LS_Selected_LUAD.xlsx');

cellIDs = txt(1, (2:end)); % Cell identifiers
geneIDs = txt((2:end), 1); % Gene identifiers
dataMatrix = num; % Data matrix (genes x cells)


% Select genes represented in the metabolic model
load('data_model/Human-GEM.mat');
model = ihuman;
selectedIndices = ismember(geneIDs, model.genes);
selectedGeneIDs = geneIDs(selectedIndices);
selectedDataMatrix = dataMatrix(selectedIndices, :);

% Transpose to cells x genes
selectedDataMatrix = selectedDataMatrix';

% Load the cell-line mapping information
cellMapping = readtable('cell_cellline.xlsx');
uniqueCellLines = unique(cellMapping{:, 2}); % Obtain the unique cell lines

% Create the cell-line label vector
cellLineLabels = cell(1, length(cellIDs));
for i = 1:length(cellIDs)
    cellID = cellIDs{i};
    matchingRow = find(strcmp(cellMapping{:, 1}, cellID), 1);
    if ~isempty(matchingRow)
        cellLineLabels{i} = cellMapping{matchingRow, 2}{1}; % Obtain the corresponding cell line
    else
        cellLineLabels{i} = 'Unknown'; % Label unmatched cells as Unknown
    end
end

% Remove unmatched cells labeled as 'Unknown'
valid_indices = ~strcmp(cellLineLabels, 'Unknown');
selectedDataMatrix = selectedDataMatrix(valid_indices, :);
cellLineLabels = cellLineLabels(valid_indices);

% Convert the cell-line labels to numeric labels
trueLabels = zeros(1, length(cellLineLabels));
for i = 1:length(cellLineLabels)
    trueLabels(i) = find(strcmp(uniqueCellLines, cellLineLabels{i}));
end

featureMatrixDouble=selectedDataMatrix;

%% step2: Hierarchical clustering analysis
Z = linkage(featureMatrixDouble, 'ward'); 

num_clusters = length(unique(trueLabels)); % Set the number of clusters to the number of unique reference labels
clusterLabels = cluster(Z, 'maxclust', num_clusters); % Cut the hierarchy according to the specified number of clusters


% Calculate the ARI (Adjusted Rand Index) 
ari = rand_index(trueLabels, clusterLabels, 'adjusted');
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

for k = 1:num_clusters
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
averageIntraDistance = totalIntraDistance / num_clusters;
averageWithinVariance = totalWithinVariance / num_clusters;

fprintf('Mean intra-cluster distance: %.4f\n', averageIntraDistance);
fprintf('Within-cluster variance: %.4f\n', averageWithinVariance);

%% Save the metrics
% Initialize the index variable
index = cell(2, 8);

% Store the metric names
index{1, 1} = 'ARI';
index{1, 2} = 'NMI';
index{1, 3} = 'Silhouette';
index{1, 4} = 'chIndex';
index{1, 5} = 'dBIndex';
index{1, 6} = 'Intra-cluster Distance';
index{1, 7} = 'Within-cluster Variance';
index{1, 8} = 'Compactness';

% Store the metric values
index{2, 1} = ari; % Adjusted Rand Index
index{2, 2} = NMI; % Normalized mutual information
index{2, 3} = averageSilhouette; % Mean silhouette coefficient
index{2, 4} = chIndex; % Calinski-Harabasz Index
index{2, 5} = dbIndex; % Davies-Bouldin Index
index{2, 6} = averageIntraDistance; % Mean intra-cluster distance
index{2, 7} = averageWithinVariance; % Within-cluster variance
index{2, 8} = compactness; % Compactness

% Display the results
disp('Metric results:');
disp(index);


%% step3: Calculate the ARI and other metrics using all genes
% Fix the random seed
rng(0);

% Load the single-cell transcriptomic data
[num, txt] = xlsread('LS_Selected_LUAD.xlsx');

cellIDs = txt(1, (2:end)); 
geneIDs = txt((2:end), 1); 
% Data matrix (genes x cells)
dataMatrix = num; 

% Transpose the matrix to cells x genes
dataMatrix = dataMatrix'; 
% Load the cell-line mapping information
cellMapping = readtable('cell_cellline.xlsx');
% Obtain the unique cell lines
uniqueCellLines = unique(cellMapping{:, 2}); 

% Create the cell-line label vector
cellLineLabels = cell(1, length(cellIDs));
for i = 1:length(cellIDs)
    cellID = cellIDs{i};
    matchingRow = find(strcmp(cellMapping{:, 1}, cellID), 1);
    if ~isempty(matchingRow)
        cellLineLabels{i} = cellMapping{matchingRow, 2}{1}; % Obtain the corresponding cell line
    else
        cellLineLabels{i} = 'Unknown'; % Label unmatched cells as Unknown
    end
end

% Remove unmatched cells labeled as 'Unknown'
valid_indices = ~strcmp(cellLineLabels, 'Unknown');
dataMatrix = dataMatrix(valid_indices, :);
cellLineLabels = cellLineLabels(valid_indices);

% Convert the cell-line labels to numeric labels
trueLabels = zeros(1, length(cellLineLabels));
for i = 1:length(cellLineLabels)
    trueLabels(i) = find(strcmp(uniqueCellLines, cellLineLabels{i}));
end

featureMatrixDouble=dataMatrix;

%% Continue with the hierarchical clustering analysis in step 2