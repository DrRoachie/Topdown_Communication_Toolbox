
%% Define data

Animal  = 'both';       % 'both' or 'MrCassius' or 'MrM'
freq_band_list  = {'beta'};

%% load histogram data arrays and max channel pairs
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

maxdir = 'F:\2024_07_24_Analysis\XCorr_Histogram_Data';

path_to_data    = 'F:\2024_07_24_Analysis\XCorr_Histogram_Data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;


%% plot channel pairs by max channels selected in Summary_PValue_Blur_Cluster_Detection.m

edges           = -200:5:200;

figure; 

% positive peaks
for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    subplot(2,length(freq_band_list),i); hold on
    
    fName = sprintf('pvalue_max_3_channels_%s_%s_PFC_AC.mat', Animal, freq_band);
    load(fullfile(maxdir, fName));
    max_results_PFC_AC = max_results;
    clear max_results;

    fName = sprintf('pvalue_max_3_channels_%s_%s_AC_PFC.mat', Animal, freq_band);
    load(fullfile(maxdir, fName));
    max_results_AC_PFC = max_results;
    clear max_results;

    max_chans       = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3);  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
                       max_results_AC_PFC(:,3) max_results_AC_PFC(:,2)]; % NOTE: since xcorr is always calculated with PFC first and AC second, the columns are flipped when taking from the AC_PFC array
    
    max_chans       = unique(max_chans,'rows'); % remove duplicate channel pairs

    % only prior
    passing_bins = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = [];

    for ii = 1:length(passing_bins)
        if (passing_bins(ii,6) > 0) && ... % if peak is positive
           (ismember(passing_bins(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
            array_to_plot = [array_to_plot; passing_bins(ii,:)];
        end
    end

    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Lag (ms)');
    ylabel('Bin Count');
    title(['Positive XCorr Peaks in ' freq_band]);
    legend show
    
    % only pretone
    array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with positive peaks and within specified layer group
    bins_to_plot = array_to_plot(:,4);     
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
end

% negative peaks
for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    subplot(2,length(freq_band_list),i+length(freq_band_list)); hold on

    fName = sprintf('pvalue_max_3_channels_%s_%s_PFC_AC.mat', Animal, freq_band);
    load(fullfile(maxdir, fName));
    max_results_PFC_AC = max_results;
    clear max_results;

    fName = sprintf('pvalue_max_3_channels_%s_%s_AC_PFC.mat', Animal, freq_band);
    load(fullfile(maxdir, fName));
    max_results_AC_PFC = max_results;
    clear max_results;

    max_chans       = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3);  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
                       max_results_AC_PFC(:,3) max_results_AC_PFC(:,2)]; % NOTE: since xcorr is always calculated with PFC first and AC second, the columns are flipped when taking from the AC_PFC array
    
    max_chans       = unique(max_chans,'rows'); % remove duplicate channel pairs
    
    % only prior
    passing_bins = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = [];

    for ii = 1:length(passing_bins)
        if (passing_bins(ii,6) < 0) && ... % if peak is negative
           (ismember(passing_bins(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
            array_to_plot = [array_to_plot; passing_bins(ii,:)];
        end
    end

    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Lag (ms)');
    ylabel('Bin Count');
    title(['Negative XCorr Peaks in ' freq_band]);
    legend show
    
    % only pretone
    array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) < 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with negative peaks and within specified layer group
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
end