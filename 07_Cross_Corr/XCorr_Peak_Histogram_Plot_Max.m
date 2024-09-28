
%% Define data

Animal  = 'MrCassius';       % 'both' or 'MrCassius' or 'MrM'
freq_band_list  = {'beta'};

%% load histogram data arrays and max channel pairs
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

maxdir          = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_24_Analysis\XCorr_Histogram_Data\COHERENCE_testTone_xcorr_histogram_data';
path_to_data    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_24_Analysis\XCorr_Histogram_Data\COHERENCE_testTone_xcorr_histogram_data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;


%%

if strcmp(Statistic, 'Coherence')
    edges = -200:5:200;
    
    figure;
    
    for i = 1:length(freq_band_list)
        freq_band = freq_band_list{i};
        
        subplot(1, length(freq_band_list), i); hold on
        
        fName = sprintf('pvalue_max_channels_%s_%s_Coh.mat', Animal, freq_band);
        load(fullfile(maxdir, fName));
        max_results_Coh = max_results;
        max_results_Coh = max_results_Coh(max_results_Coh(:, 1) == 1, :);
        
        max_chans = [max_results_Coh(:,2) max_results_Coh(:,3)];
        
        % Only prior
        passing_bins_OnlyPrior = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
        array_to_plot_prior = [];
        
        for ii = 1:length(passing_bins_OnlyPrior)
            if passing_bins_OnlyPrior(ii, 6) && ...
               (ismember(passing_bins_OnlyPrior(ii, 2:3), max_chans, 'rows'))
                array_to_plot_prior = [array_to_plot_prior; passing_bins_OnlyPrior(ii,:)];
            end
        end
        
        bins_to_plot_OnlyPrior = array_to_plot_prior(:, 4);
        
        % Only pretone
        passing_bins_OnlyPretone = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
        array_to_plot_pretone = [];
        
        for ii = 1:length(passing_bins_OnlyPretone)
            if passing_bins_OnlyPretone(ii, 6) && ...
               (ismember(passing_bins_OnlyPretone(ii, 2:3), max_chans, 'rows'))
                array_to_plot_pretone = [array_to_plot_pretone; passing_bins_OnlyPretone(ii,:)];
            end
        end
        
        bins_to_plot_OnlyPretone = array_to_plot_pretone(:, 4);
        
        % Calculate histogram counts without plotting
        counts_prior = histcounts(bins_to_plot_OnlyPrior, edges);
        counts_pretone = histcounts(bins_to_plot_OnlyPretone, edges);
        
        % Find the maximum count value between the two histograms
        max_count = max([max(counts_prior), max(counts_pretone)]);
        
        % Normalize both histograms to the maximum count
        norm_counts_prior = (counts_prior / max_count) * 100;
        norm_counts_pretone = (counts_pretone / max_count) * 100;
        
        % Plot the normalized histograms using bar
        % PriorOnly in red
        bar(edges(1:end-1), norm_counts_prior, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
        
        % PretoneOnly in blue
        bar(edges(1:end-1), norm_counts_pretone, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
        
        % Formatting
        xlim([-200 200]);
        xlabel('Lag (ms)');
        ylabel('Normalized Percentage');
        title([Animal ' Absolute XCorr Peaks in ' freq_band]);
        legend show
    end
end


%% Collects XCorr Peaks from the combined Granger clusters. Basically smooshes the strongest PFC2AC and AC2PFC clusters together and drops repeats. 

if strcmp(Statistic, 'Granger')

        edges           = -200:5:200;
        
        figure; 
        
        
        for i = 1:length(freq_band_list)
        
            freq_band = freq_band_list{i};
        
            subplot(1, length(freq_band_list),i); hold on

        
            fName = sprintf('pvalue_max_channels_%s_%s_%s_PFC_AC.mat', Statistic, Animal, freq_band);
            load(fullfile(maxdir, fName));
            max_results_PFC_AC = max_results;                                          % top three clusters 
            max_results_PFC_AC = max_results_PFC_AC(max_results_PFC_AC(:, 1)==1, :);   % the max or strongest connection
            clear max_results;               
        
            fName = sprintf('pvalue_max_channels_%s_%s_%s_AC_PFC.mat', Statistic, Animal, freq_band);
            load(fullfile(maxdir, fName));
            max_results_AC_PFC = max_results;                                          % top three clusters 
            max_results_AC_PFC = max_results_AC_PFC(max_results_AC_PFC(:, 1)==1, :);   % the max or strongest connection
            clear max_results;
        
            max_chans       = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3);  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
                               max_results_AC_PFC(:,3) max_results_AC_PFC(:,2)]; % IMPORTANT NOTE: since xcorr is always calculated with PFC first and AC second, the columns are flipped when taking from the AC_PFC array
        
            max_chans       = unique(max_chans,'rows'); % remove duplicate channel pairs
        
            % only prior
            passing_bins_OnlyPrior = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
            array_to_plot = [];
        
            for ii = 1:length(passing_bins_OnlyPrior)
                if passing_bins_OnlyPrior(ii,6) && ... % get all peaks 
                   (ismember(passing_bins_OnlyPrior(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
                    array_to_plot = [array_to_plot; passing_bins_OnlyPrior(ii,:)];
                end
            end
        
            bins_to_plot_OnlyPrior = array_to_plot(:,4);
        
            histogram( bins_to_plot_OnlyPrior, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
            xlim([-200 200]);
            xlabel('Lag (ms)');
            ylabel('Bin Count');
            title(['Positive XCorr Peaks in ' freq_band]);
            legend show
        
            % only pretone
        
            passing_bins_OnlyPretone = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
            array_to_plot = [];   % find bins with positive peaks and within specified layer group
        
            for ii = 1:length(passing_bins_OnlyPretone)
                if passing_bins_OnlyPretone(ii,6) && ... % get all peaks 
                   (ismember(passing_bins_OnlyPretone(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
                    array_to_plot = [array_to_plot; passing_bins_OnlyPretone(ii,:)];
                end
            end
        
            bins_to_plot_OnlyPretone = array_to_plot(:,4);
        
            histogram(bins_to_plot_OnlyPretone, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
            xlim([-200 200]);
            xlabel('Lag (ms)');
            ylabel('Bin Count');
            title([Animal ' Absolute XCorr Peaks in ' freq_band]);
            legend show
        
        end

end
%% Optional chunk for separating the positive and negative peaks 

% %% positive peaks
% 
% figure(1)
% for i = 1:length(freq_band_list)
% 
%     freq_band = freq_band_list{i};
% 
%     subplot(2,length(freq_band_list),i); hold on
% 
%     fName = sprintf('pvalue_max_channels_%s_%s_PFC_AC.mat', Animal, freq_band);
%     load(fullfile(maxdir, fName));
%     max_results_PFC_AC = max_results;
%     clear max_results;
% 
%     fName = sprintf('pvalue_max_channels_%s_%s_AC_PFC.mat', Animal, freq_band);
%     load(fullfile(maxdir, fName));
%     max_results_AC_PFC = max_results;
%     clear max_results;
% 
%     max_chans       = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3);  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
%                        max_results_AC_PFC(:,3) max_results_AC_PFC(:,2)]; % IMPORTANT NOTE: since xcorr is always calculated with PFC first and AC second, the columns are flipped when taking from the AC_PFC array
% 
%     max_chans       = unique(max_chans,'rows'); % remove duplicate channel pairs
% 
%     % only prior
%     passing_bins_OnlyPrior = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
%     array_to_plot = [];
% 
%     for ii = 1:length(passing_bins_OnlyPrior)
%         if (passing_bins_OnlyPrior(ii,6) > 0) && ... % if peak is positive
%            (ismember(passing_bins_OnlyPrior(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
%             array_to_plot = [array_to_plot; passing_bins_OnlyPrior(ii,:)];
%         end
%     end
% 
%     bins_to_plot_OnlyPrior = array_to_plot(:,4);
% 
%     histogram( bins_to_plot_OnlyPrior, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
%     xlim([-200 200]);
%     xlabel('Lag (ms)');
%     ylabel('Bin Count');
%     title(['Positive XCorr Peaks in ' freq_band]);
%     legend show
% 
%     % only pretone
% 
%     passing_bins_OnlyPretone = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
%     array_to_plot = [];   % find bins with positive peaks and within specified layer group
% 
%     for ii = 1:length(passing_bins_OnlyPretone)
%        if (passing_bins_OnlyPretone(ii,6) > 0) && ... % if peak is positive
%            (ismember(passing_bins_OnlyPretone(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
%             array_to_plot = [array_to_plot; passing_bins_OnlyPretone(ii,:)];
%         end
%     end
% 
%     bins_to_plot_OnlyPretone = array_to_plot(:,4);
% 
%     histogram(bins_to_plot_OnlyPretone, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
%     xlim([-200 200]);
%     xlabel('Lag (ms)');
%     ylabel('Bin Count');
%     title(['Positive XCorr Peaks in ' freq_band]);
%     legend show
% 
% end
% 
% 
% 
% %% negative peaks
% 
% figure(2)
% 
% for i = 1:length(freq_band_list)
% 
%     freq_band = freq_band_list{i};
% 
%     subplot(2,length(freq_band_list),i); hold on
% 
%     fName = sprintf('pvalue_max_channels_%s_%s_PFC_AC.mat', Animal, freq_band);
%     load(fullfile(maxdir, fName));
%     max_results_PFC_AC = max_results;
%     clear max_results;
% 
%     fName = sprintf('pvalue_max_channels_%s_%s_AC_PFC.mat', Animal, freq_band);
%     load(fullfile(maxdir, fName));
%     max_results_AC_PFC = max_results;
%     clear max_results;
% 
%     max_chans       = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3);  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
%                        max_results_AC_PFC(:,3) max_results_AC_PFC(:,2)]; % IMPORTANT NOTE: since xcorr is always calculated with PFC first and AC second, the columns are flipped when taking from the AC_PFC array
% 
%     max_chans       = unique(max_chans,'rows'); % remove duplicate channel pairs
% 
%     % only prior
%     passing_bins_OnlyPrior = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
%     array_to_plot = [];
% 
%     for ii = 1:length(passing_bins_OnlyPrior)
%         if (passing_bins_OnlyPrior(ii,6) < 0) && ... % if peak is negative
%            (ismember(passing_bins_OnlyPrior(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
%             array_to_plot = [array_to_plot; passing_bins_OnlyPrior(ii,:)];
%         end
%     end
% 
%     bins_to_plot_OnlyPrior = array_to_plot(:,4);
% 
%     histogram( bins_to_plot_OnlyPrior, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
%     xlim([-200 200]);
%     xlabel('Lag (ms)');
%     ylabel('Bin Count');
%     title(['Positive XCorr Peaks in ' freq_band]);
%     legend show
% 
%     % only pretone
% 
%     passing_bins_OnlyPretone = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
%     array_to_plot = [];   % find bins with positive peaks and within specified layer group
% 
%     for ii = 1:length(passing_bins_OnlyPretone)
%         if (passing_bins_OnlyPretone(ii,6) < 0) && ... % if peak is negative
%            (ismember(passing_bins_OnlyPretone(ii,2:3), max_chans, 'rows')) % if channel is in max_chans array
%             array_to_plot = [array_to_plot; passing_bins_OnlyPretone(ii,:)];
%         end
%     end
% 
%     bins_to_plot_OnlyPretone = array_to_plot(:,4);
% 
%     histogram(bins_to_plot_OnlyPretone, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
%     xlim([-200 200]);
%     xlabel('Lag (ms)');
%     ylabel('Bin Count');
%     title(['Negative XCorr Peaks in ' freq_band]);
%     legend show
% 
% end