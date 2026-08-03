%% Three-factor ANOVA
% Analyze each dataset separately
MEM = repelem({'ftINIT', 'GIMME', 'iMAT'}, 8)';
Data = repmat([repmat({'LS'},1,4), repmat({'R'},1,4)], 1, 3)';
Threshold = repmat({'25th','50th','75th','90th'}, 1, 6)';

% Enter the scores of 24 strategies for each of the nine datasets separately, excluding rFASTCORMICS strategies
Score = zeros(24,1);

% Construct the table
data = table(MEM, Data, Threshold, Score, ...
    'VariableNames', {'MEM', 'Data', 'Threshold', 'Score'});

% Convert grouping variables to categorical variables for ANOVA
data.MEM = categorical(data.MEM);
data.Data = categorical(data.Data);
data.Threshold = categorical(data.Threshold);

% Fit a main-effects linear model without interaction terms
mdl = fitlm(data, 'Score ~ MEM + Data + Threshold');

[p, tbl, stats] = anovan(Score, ...
    {MEM, Data, Threshold}, ...
    'model', 'linear', ...
    'varnames', {'MEM','Data','Threshold'});

disp('ANOVA table:');
disp(tbl);

%% Post-hoc multiple comparisons
factorNames = {'MEM','Data','Threshold'};
for i = 1:3
    if p(i) < 0.05
        fprintf('\nPost-hoc for %s:\n', factorNames{i});

        [c,~,~,g] = multcompare(stats, ...
            'Dimension', i, ...
            'Display', 'off');

        sig = c(c(:,6) < 0.05, :);

        if isempty(sig)
            disp('  No significant pairwise differences.');
        else
            for k = 1:size(sig,1)
                fprintf('  %s vs %s : p = %.4f\n', ...
                    g{sig(k,1)}, g{sig(k,2)}, sig(k,6));
            end
        end
    end
end


%% Analysis across datasets
Score = zeros(24,9);

MEM = repelem({'ftINIT','GIMME','iMAT'},8)';
Data = repmat([repmat({'LS'},1,4), repmat({'R'},1,4)],1,3)';
Threshold = repmat({'25th','50th','75th','90th'},1,6)';

% Construct the table
data_show = table(MEM, Data, Threshold, Score, ...
    'VariableNames', {'MEM', 'Data', 'Threshold', 'Score'});

% ===== Expand repeated observations =====
nRepeat = size(Score,2);   
Score_long     = reshape(Score.', [], 1);   % Reshape in row-major order
MEM_long       = repelem(MEM, nRepeat);
Data_long      = repelem(Data, nRepeat);
Threshold_long = repelem(Threshold, nRepeat);

% ===== Construct the table =====
data = table(MEM_long, Data_long, Threshold_long, Score_long, ...
    'VariableNames', {'MEM','Data','Threshold','Score'});

data.MEM       = categorical(data.MEM);
data.Data      = categorical(data.Data);
data.Threshold = categorical(data.Threshold);

% ===== Main-effects linear model with three factors =====
mdl = fitlm(data, 'Score ~ MEM + Data + Threshold');

anovaTable = anova(mdl);
disp(anovaTable);

% ===== ANOVA for post-hoc comparisons =====
[p,~,stats] = anovan(data.Score, ...
    {data.MEM, data.Data, data.Threshold}, ...
    'model','linear', ...
    'varnames',{'MEM','Data','Threshold'}, ...
    'display','off');

% ===== Post-hoc comparisons =====
names = {'MEM','Data','Threshold'};

for i = 1:3
    if p(i) < 0.05
        fprintf('\nPost-hoc for %s:\n', names{i});
        [c,~,~,g] = multcompare(stats,'Dimension',i,'Display','off');
        c = c(c(:,6)<0.05,:);
        for k = 1:size(c,1)
            fprintf('%s vs %s : p = %.4f\n', ...
                g{c(k,1)}, g{c(k,2)}, c(k,6));
        end
    end
end
% Save anovaTable