
%% load histogram data arrays
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

path_to_data    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;


%%

data1 = bins_to_plot_OnlyPrior;
data2 = bins_to_plot_OnlyPretone;

% Define the bin edges (or let MATLAB determine them)
bin_edges = linspace(min([data1; data2]), max([data1; data2]), 50);

% Compute the histograms
[counts1, edges1] = histcounts(data1, bin_edges, 'Normalization', 'pdf');
[counts2, edges2] = histcounts(data2, bin_edges, 'Normalization', 'pdf');

% Calculate the bin centers for plotting
bin_centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
bin_centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;

figure;
plot(bin_centers1, counts1, '-r', 'LineWidth', 2); % Red line for the first array
hold on;
plot(bin_centers2, counts2, '-b', 'LineWidth', 2); % Blue line for the second array
hold off;

% Add labels and legend
xlabel('Value');
ylabel('Probability Density');
legend('Data 1', 'Data 2');
title('Histogram Shape Comparison');

%% run a statistical test between the two hisgram shapes 

% Define the bin edges
num_bins = 50;
bin_edges = linspace(min([data1; data2]), max([data1; data2]), num_bins);

% Compute the histograms (without normalization)
counts1 = histcounts(data1, bin_edges);
counts2 = histcounts(data2, bin_edges);

% Combine the counts into a single matrix for the test
observed = [counts1; counts2];
observed = observed + 0.5; % Add small counts to avoid zeros

% Perform Chi-Squared test
[h, p_value, stats] = chi2gof(observed(:), 'Expected', mean(observed(:)));

% Display results
if h == 0
    fprintf('The histograms are not significantly different (p-value = %.4f).\n', p_value);
else
    fprintf('The histograms are significantly different (p-value = %.4f).\n', p_value);
end

% Calculate the bin centers for plotting
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% Plot the histogram shapes as lines
figure;
plot(bin_centers, counts1, '-r', 'LineWidth', 2); % Red line for the first array
hold on;
plot(bin_centers, counts2, '-b', 'LineWidth', 2); % Blue line for the second array
hold off;

% Add labels and legend
xlabel('Value');
ylabel('Count');
legend('Data 1', 'Data 2');
title('Histogram Shape Comparison');

%%

% Use MATLAB's ksdensity function to estimate the PDFs
[pdf1, x_values1] = ksdensity(data1);
[pdf2, x_values2] = ksdensity(data2);

figure;
plot(x_values1, pdf1, '-r', 'LineWidth', 2); % Red line for the first array
hold on;
plot(x_values2, pdf2, '-b', 'LineWidth', 2); % Blue line for the second array
hold off;

% Add labels and legend
xlabel('Value');
ylabel('Probability Density');
legend('Data 1', 'Data 2');
title('Kernel Density Estimate Comparison');

%%

% Compute the KDEs and CDFs
[pdf1, x_values1, cdf1] = ksdensity(data1, 'Function', 'cdf');
[pdf2, x_values2, cdf2] = ksdensity(data2, 'Function', 'cdf');

% Perform the KS test
[h, p_value, ks_stat] = kstest2(data1, data2);

% Display results
if h == 0
    fprintf('The two distributions are not significantly different (p-value = %.4f).\n', p_value);
else
    fprintf('The two distributions are significantly different (p-value = %.4f).\n', p_value);
end

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