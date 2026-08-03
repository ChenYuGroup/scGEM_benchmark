%% Calculate percentile-based expression cutoffs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
[num,txt] = xlsread('LS_LUAD.xlsx');
allData =num(:);

% Percentiles evaluated in the benchmark.
percentiles = [25, 50, 75, 90]; 

percentileValues = prctile(allData, percentiles);















