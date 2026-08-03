%% Computational feasibility analysis
clear; clc;

% Using the LUAD dataset and the ftINIT-LS-25th strategy as an example.
datasetName = 'LUAD';
memName = 'ftINIT';
strategyName = 'ftINIT-LS-25th';

load('models_ftINIT_LS_25th.mat', 'models', 'times');

% Convert times to a numeric column vector
if iscell(times)
    times = cell2mat(times);
end
times = times(:);

% Identify valid models
validModel = cellfun(@(m) ...
    isstruct(m) && isfield(m,'rxns') && ~isempty(m.rxns), ...
    models);
validModel = validModel(:);

if numel(times) ~= numel(validModel)
    error('The numbers of models and runtime records do not match.');
end

% Set the runtime limit
switch lower(memName)
    case 'ftinit'
        timeoutLimit = 5000;
    case 'imat'
        timeoutLimit = 7200;
    otherwise
        timeoutLimit = Inf;
end

% Calculate feasibility statistics
timeoutOccurred = times >= timeoutLimit;

nModels = numel(models);
nValidModels = sum(validModel);
nEmptyModels = nModels - nValidModels;
nTimeouts = sum(timeoutOccurred);
successRate = nValidModels / nModels;
meanTime = mean(times(validModel), 'omitnan');

summaryTable = table( ...
    string(datasetName), string(strategyName), string(memName), ...
    nModels, nValidModels, nEmptyModels, nTimeouts, ...
    successRate, meanTime, ...
    'VariableNames', { ...
    'Dataset','Strategy','MEM','NumberOfModels', ...
    'NumberOfValidModels','NumberOfEmptyModels', ...
    'NumberOfTimeouts','SuccessRate','MeanTimeSeconds'});

disp(summaryTable);
