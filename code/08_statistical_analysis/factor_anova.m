%% Three-factor ANOVA
%% 1. 定义24种策略及其顺序
MEM = repelem({'ftINIT', 'GIMME', 'iMAT'}, 8)';
Data = repmat([repmat({'LS'},1,4), repmat({'R'},1,4)], 1, 3)';
Threshold = repmat({'25th','50th','75th','90th'}, 1, 6)';

%% 2. 输入9个数据集的分数
% 每列对应一个数据集；每行顺序与MEM、Data、Threshold一致
% 请将下面的零矩阵替换为实际分数
Score = zeros(24, 9);

datasetNames = compose("Dataset%d", 1:9);

assert(isequal(size(Score), [24, 9]), ...
    'Score必须是24×9矩阵。');

%% 3. 初始化结果
factorOrder = {'MEM', 'Data', 'Threshold', 'Residual'};

ExplainedPercent_matrix = NaN(4, 9);
PartialEtaSquared_matrix = NaN(4, 9);
SumSq_matrix = NaN(4, 9);
R2_percent = NaN(1, 9);

varianceResults_all = cell(1, 9);
anovaTables = cell(1, 9);
models = cell(1, 9);

%% 4. 依次分析9个数据集
for d = 1:9

    fprintf('\n========== %s ==========\n', char(datasetNames(d)));

    % 当前数据集的24个分数
    currentScore = Score(:, d);

    % 排除缺失或非有限数值
    valid = isfinite(currentScore);

    if sum(valid) < 2
        warning('%s 有效分数不足，跳过。', char(datasetNames(d)));
        continue;
    end

    if sum(valid) < 24
        warning('%s 仅有%d个有效分数，因素设计可能不平衡。', ...
            char(datasetNames(d)), sum(valid));
    end

    if all(currentScore(valid) == currentScore(find(valid, 1)))
        warning('%s 分数没有变异，无法计算占比，结果记为NaN。', ...
            char(datasetNames(d)));
        continue;
    end

    try
        %% 构建当前数据集的表
        data = table( ...
            categorical(MEM(valid)), ...
            categorical(Data(valid)), ...
            categorical(Threshold(valid)), ...
            currentScore(valid), ...
            'VariableNames', {'MEM', 'Data', 'Threshold', 'Score'});

        %% 主效应模型：不包含交互项
        mdl = fitlm(data, 'Score ~ MEM + Data + Threshold');

        % 检查因素效应是否可以独立估计
        if mdl.NumEstimatedCoefficients < 7
            error('有效数据不足以独立估计全部因素效应。');
        end

        %% ANOVA：显式使用Type III平方和
        [p, tbl, stats] = anovan(data.Score, ...
            {data.MEM, data.Data, data.Threshold}, ...
            'model', 'linear', ...
            'sstype', 3, ...
            'varnames', {'MEM', 'Data', 'Threshold'}, ...
            'display', 'off');

        disp('ANOVA table:');
        disp(tbl);

        %% 提取平方和
        factorSS = cell2mat(tbl(2:4, 2));
        errorSS  = tbl{5, 2};
        totalSS  = tbl{6, 2};

        assert(totalSS > 0, ...
            'Score没有变异，无法计算解释方差比例。');

        %% 因素和残差占总方差的比例
        Factor = string(factorOrder(:));
        SumSq  = [factorSS; errorSS];
        ExplainedPercent = 100 * SumSq / totalSS;

        % 残差不计算偏eta平方
        PartialEtaSquared = [
            factorSS ./ (factorSS + errorSS);
            NaN
        ];

        varianceResults = table( ...
            Factor, SumSq, ExplainedPercent, PartialEtaSquared);

        disp('各因素及残差占总方差的比例（%）：');
        disp(varianceResults);

        %% 保存至工作区结果变量
        ExplainedPercent_matrix(:, d)  = ExplainedPercent;
        PartialEtaSquared_matrix(:, d) = PartialEtaSquared;
        SumSq_matrix(:, d)             = SumSq;
        R2_percent(d)                 = 100 * mdl.Rsquared.Ordinary;

        varianceResults_all{d} = varianceResults;
        anovaTables{d} = tbl;
        models{d} = mdl;

    catch ME
        warning('%s 分析失败：%s；继续下一数据集。', ...
            char(datasetNames(d)), ME.message);
    end
end

%% 5. 汇总：4行×9列
ExplainedPercent_table = array2table( ...
    ExplainedPercent_matrix, ...
    'VariableNames', cellstr(datasetNames), ...
    'RowNames', factorOrder);

PartialEtaSquared_table = array2table( ...
    PartialEtaSquared_matrix, ...
    'VariableNames', cellstr(datasetNames), ...
    'RowNames', factorOrder);

SumSq_table = array2table( ...
    SumSq_matrix, ...
    'VariableNames', cellstr(datasetNames), ...
    'RowNames', factorOrder);

R2_table = array2table( ...
    R2_percent, ...
    'VariableNames', cellstr(datasetNames), ...
    'RowNames', {'Model_R2_percent'});

fprintf('\n========== 9个数据集的方差占比（%%） ==========\n');
disp(ExplainedPercent_table);

fprintf('\n========== 偏eta平方 ==========\n');
disp(PartialEtaSquared_table);

fprintf('\n========== 模型整体R²（%%） ==========\n');
disp(R2_table);


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