
%% load histogram data arrays
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

path_to_data    = 'F:\2024_07_24_Analysis\XCorr_Histogram_Data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;

%% extract bin counts of desired frequency, sign, layer group

superficial = 1:5;
uppermid    = 6:10;
lowermid    = 11:15;
deep        = 16:20;

freq_band       = 'theta';
sign            = 'positive';
pfc_layer_group = superficial;
ac_layer_group  = superficial;
edges           = -200:1:200;

% only prior
array_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with positive peaks and within specified layer group
bins_to_plot = array_to_plot(:,4);

OnlyPrior_hist_counts = histcounts(bins_to_plot, edges);

% only pretone
array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with positive peaks and within specified layer group
bins_to_plot = array_to_plot(:,4);

OnlyPretone_hist_counts = histcounts(bins_to_plot, edges);

%% perform chi square test

% Example data creation
rows = 2; % number of categories for the first variable
columns = 400; % number of categories for the second variable
n = 1000; % number of observations

% Random categorical data
var1 = randi([1, rows], n, 1); % First variable (e.g., 1 or 2)
var2 = randi([1, columns], n, 1); % Second variable (e.g., 1 to 400)

% Create the contingency table
[contingencyTable, ~, ~, labels] = crosstab(var1, var2);

% Compute the row and column sums
rowSums = sum(contingencyTable, 2);
columnSums = sum(contingencyTable, 1);
grandTotal = sum(rowSums);

% Calculate expected frequencies
expected = (rowSums * columnSums) / grandTotal;

% Compute the Chi-Square statistic
chiSquareStat = sum((contingencyTable(:) - expected(:)).^2 ./ expected(:));

% Degrees of freedom
df = (rows - 1) * (columns - 1);

% Compute the p-value
pValue = 1 - chi2cdf(chiSquareStat, df);