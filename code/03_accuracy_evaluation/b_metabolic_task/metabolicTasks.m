% Refer to GitHub\ftINIT-RAVEN
% load scGEMs
cd ..\scGEM_benchmark\output_all\LUAD\models;
load('models_ftINIT_LS_25th.mat');

% Load the essential metabolic task list
essTasks = parseTaskList('..\scGEM_benchmark\data\reference_data\metabolic_tasks\metabolicTasks_Essential.txt');
allTasks = essTasks;

taskReports = cell(length(models),1);
parfor i = 1:length(models)
    disp(['Running: ' num2str(i)])
    [taskReports{i}, ~]=checkTasks(closeModel(models{i}),[],true,false,false,allTasks);
end

%% Calculate the essential metabolic task completion rate
N = length(taskReports);              
results = zeros(numel(allTasks), N);        

% Extract the task completion status for each model
for i = 1:N
    results(:, i) = taskReports{i}.ok;
end

% Calculate the number and percentage of completed tasks for each model
counts = sum(results, 1)';               
percentages = counts / numel(allTasks) * 100;       
tasks_meanper = mean(percentages);