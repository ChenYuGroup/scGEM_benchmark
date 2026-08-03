%% 1. Randomly select samples
% Read the XLSX file
filename = 'LS_LUAD.xlsx';
data = readtable(filename);

geneIDs = data{:,1};
expressionData = data{:, 2:end};

sampleIDs = data.Properties.VariableNames(2:end); % Assume that the first column contains gene IDs, so samples start from the second column

% Randomly select the column indices of 100 samples
randomIndices = randperm(length(sampleIDs), 100); % Save the selected indices in randomIndices

% Extract the expression data and sample IDs of the selected samples
selectedExpressionData = expressionData(:, randomIndices);
selectedSampleIDs = sampleIDs(randomIndices);

% Create a new table containing gene IDs
T = array2table(selectedExpressionData, 'VariableNames', selectedSampleIDs);
T = [geneIDs, T]; % Add the gene IDs back to the table as the first column

% Define the output XLSX filename
outputFileName = 'LS_Selected_LUAD.xlsx';

% Write the table to a new XLSX file, including gene IDs as the first column
writetable(T, outputFileName, 'WriteVariableNames', true);

%% 2. Read the second file and select the same samples as in the first file
load('randomIndices_LUAD.mat');

filename2 = 'R_LUAD.xlsx';
data2 = readtable(filename2);

geneIDs2 = data2{:,1};
expressionData2 = data2{:, 2:end};

% Extract the expression data and sample IDs using the previously generated random indices
selectedExpressionData2 = expressionData2(:, randomIndices);
selectedSampleIDs2= data2.Properties.VariableNames(randomIndices + 1); % Add 1 because MATLAB indexing starts at 1 and the first column is excluded

% Create a new table containing gene IDs
T2 = array2table(selectedExpressionData2, 'VariableNames', selectedSampleIDs2);
T2 = [array2table(geneIDs2, 'VariableNames', {'GeneID'}), T2]; % Add the gene IDs back to the table as the first column

% Define the output XLSX filename
outputFileName2 = 'R_Selected_LUAD.xlsx';
writetable(T2, outputFileName2, 'WriteVariableNames', true)