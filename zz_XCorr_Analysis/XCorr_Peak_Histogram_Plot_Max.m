
%% Define data

Animal  = 'MrCassius';       % 'both' or 'MrCassius' or 'MrM'
freq_band_list  = {'theta'};

%% load histogram data arrays and max channel pairs
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

maxdir          = 'D:\2024_09_27_Analysis\XCorr_Histogram_Data\Coherence_preCue_xcorr_histogram_data';
path_to_data    = 'D:\2024_09_27_Analysis\XCorr_Histogram_Data\Coherence_preCue_xcorr_histogram_data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;


%%

if strcmp(Statistic, 'Coherence')
    edges = -200:5:200;
    
    % figure (1);
  
    % 
    % for i = 1:length(freq_band_list)
    %     freq_band = freq_band_list{i};
    % 
    %     subplot(1, length(freq_band_list), i); hold on
    % 
    %     fName = sprintf('pvalue_max_channels_%s_%s_Coh.mat', Animal, freq_band);
    %     load(fullfile(maxdir, fName));
    %     max_results_Coh = max_results;
    %     max_results_Coh = max_results_Coh(max_results_Coh(:, 1) == 1, :);
    % 
    %     max_chans = [max_results_Coh(:,2) max_results_Coh(:,3)];
    % 
    %     % Only prior
    %    
    %     passing_bins_OnlyPrior    = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    %     array_to_plot_prior = [];
    % 
    %     for ii = 1:length(passing_bins_OnlyPrior)
    %         if passing_bins_OnlyPrior(ii, 6) && ...
    %            (ismember(passing_bins_OnlyPrior(ii, 2:3), max_chans, 'rows'))
    %             array_to_plot_prior = [array_to_plot_prior; passing_bins_OnlyPrior(ii,:)];
    %         end
    %     end
    % 
    %     bins_to_plot_OnlyPrior = array_to_plot_prior(:, 4);
    % 
    %     % Only pretone
    %     TotalTimeBins_OnlyPretone       = length(OnlyPretone_Time_Bins.(freq_band).alltimebins);
    %     PassingTimeBins_OnlyPretone     = length(OnlyPretone_Time_Bins.(freq_band).passingtimebins);
    %     PretoneOnly_Xcorr_Percent_Pass  = (PassingTimeBins_OnlyPretone/TotalTimeBins_OnlyPretone) * 100;
    % 
    %     passing_bins_OnlyPretone = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    %     array_to_plot_pretone = [];
    % 
    %     for ii = 1:length(passing_bins_OnlyPretone)
    %         if passing_bins_OnlyPretone(ii, 6) && ...
    %            (ismember(passing_bins_OnlyPretone(ii, 2:3), max_chans, 'rows'))
    %             array_to_plot_pretone = [array_to_plot_pretone; passing_bins_OnlyPretone(ii,:)];
    %         end
    %     end
    % 
    %     bins_to_plot_OnlyPretone = array_to_plot_pretone(:, 4);
    % 
    %     % Calculate histogram counts without plotting
    %     counts_prior = histcounts(bins_to_plot_OnlyPrior, edges);
    %     counts_pretone = histcounts(bins_to_plot_OnlyPretone, edges);
    % 
    %     % Find the maximum count value between the two histograms
    %     max_count = max([max(counts_prior), max(counts_pretone)]);
    % 
    %    % Calculate the total number of counts across both conditions
    %     total_counts_all = sum(counts_prior) + sum(counts_pretone);
    % 
    %     % Normalize both histograms to the total number of counts across conditions
    %     norm_counts_prior = (counts_prior / total_counts_all) * 100;
    %     norm_counts_pretone = (counts_pretone / total_counts_all) * 100;
    % 
    %     % Plot PriorOnly as an outline in red
    %     plot(edges(1:end-1), norm_counts_prior, 'r-', 'LineWidth', 2, 'DisplayName', 'PriorOnly');
    % 
    %     % Plot PretoneOnly as an outline in blue
    %     plot(edges(1:end-1), norm_counts_pretone, 'b-', 'LineWidth', 2, 'DisplayName', 'PretoneOnly');
    % 
    % 
    %     % Formatting
    %     xlim([-200 200]);
    %     xlabel('Lag (ms)', 'FontSize', 14);
    %     ylabel('Percentage of Trials','FontSize', 14);
    %     %title([Animal ' Absolute XCorr Peaks in ' freq_band]);
    %     legend('LED Only Trials', 'Pretone Only Trials'); 
    % end

    figure (2);
  
    for i = 1:length(freq_band_list)
        freq_band = freq_band_list{i};

        subplot(1, length(freq_band_list), i); hold on

        fName = sprintf('pvalue_max_channels_%s_%s_Coh.mat', Animal, freq_band);
        load(fullfile(maxdir, fName));
        max_results_Coh = max_results;
        max_results_Coh = max_results_Coh(max_results_Coh(:, 1) == 1, :);

        max_chans = [max_results_Coh(:,2) max_results_Coh(:,3)];

        % Only prior

        all_bins_OnlyPrior      = OnlyPrior_Time_Bins.(freq_band).alltimebins;
        passing_bins_OnlyPrior    = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
        Max_PassedBins_OnlyPrior = [];
        Max_AllBins_OnlyPrior    = [];

        for ii = 1:length(passing_bins_OnlyPrior)
            if passing_bins_OnlyPrior(ii, 6) && ...
            (ismember(passing_bins_OnlyPrior(ii, 2:3), max_chans, 'rows'))
            Max_PassedBins_OnlyPrior = [Max_PassedBins_OnlyPrior; passing_bins_OnlyPrior(ii,:)];
            end
        end
         
        for rr = 1:length(all_bins_OnlyPrior)
            if all_bins_OnlyPrior(rr, 6) && ...
            (ismember(all_bins_OnlyPrior(rr, 2:3), max_chans, 'rows'))
            Max_AllBins_OnlyPrior = [Max_AllBins_OnlyPrior; all_bins_OnlyPrior(rr,:)];
            end
        end

         Num_Max_PassedBins_OnlyPrior     = length(Max_PassedBins_OnlyPrior);
         Num_Max_AllBins_OnlyPrior        = length(Max_AllBins_OnlyPrior);
         Percentage_Passed_Bins_PriorOnly = (Num_Max_PassedBins_OnlyPrior/Num_Max_AllBins_OnlyPrior) * 100;

        % Only pretone
    
        all_bins_OnlyPretone     = OnlyPretone_Time_Bins.(freq_band).alltimebins;
        passing_bins_OnlyPretone    = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
        Max_PassedBins_OnlyPretone = [];
        Max_AllBins_OnlyPretone    = [];

        for ii = 1:length(passing_bins_OnlyPretone)
            if passing_bins_OnlyPretone(ii, 6) && ...
            (ismember(passing_bins_OnlyPretone(ii, 2:3), max_chans, 'rows'))
            Max_PassedBins_OnlyPretone = [Max_PassedBins_OnlyPretone; passing_bins_OnlyPretone(ii,:)];
            end
        end
         
        for rr = 1:length(all_bins_OnlyPretone)
            if all_bins_OnlyPretone(rr, 6) && ...
            (ismember(all_bins_OnlyPretone(rr, 2:3), max_chans, 'rows'))
             Max_AllBins_OnlyPretone = [Max_AllBins_OnlyPretone; all_bins_OnlyPretone(rr,:)];
            end
        end

         Num_Max_PassedBins_OnlyPretone    = length(Max_PassedBins_OnlyPretone);
         Num_Max_AllBins_OnlyPretone       = length(Max_AllBins_OnlyPretone);
         Percentage_Passed_Bins_PretoneOnly = (Num_Max_PassedBins_OnlyPretone/Num_Max_AllBins_OnlyPretone) * 100;

% Calculate counts for each condition
Num_Max_NotPassedBins_OnlyPrior = Num_Max_AllBins_OnlyPrior - Num_Max_PassedBins_OnlyPrior;
Num_Max_NotPassedBins_OnlyPretone = Num_Max_AllBins_OnlyPretone - Num_Max_PassedBins_OnlyPretone;

% Create the observed counts
observed_counts = [
    Num_Max_PassedBins_OnlyPrior, Num_Max_NotPassedBins_OnlyPrior; 
    Num_Max_PassedBins_OnlyPretone, Num_Max_NotPassedBins_OnlyPretone
];

% Calculate expected counts
total_passed = Num_Max_PassedBins_OnlyPrior + Num_Max_PassedBins_OnlyPretone;
total_not_passed = Num_Max_NotPassedBins_OnlyPrior + Num_Max_NotPassedBins_OnlyPretone;
total = total_passed + total_not_passed;

expected_counts = [
    (total_passed / total) * (Num_Max_AllBins_OnlyPrior), (total_not_passed / total) * (Num_Max_AllBins_OnlyPrior);
    (total_passed / total) * (Num_Max_AllBins_OnlyPretone), (total_not_passed / total) * (Num_Max_AllBins_OnlyPretone)
];

% Calculate chi-squared statistic
chi2_stat = sum((observed_counts(:) - expected_counts(:)).^2 ./ expected_counts(:));

% Degrees of freedom (df) for a 2x2 contingency table
df = (size(observed_counts, 1) - 1) * (size(observed_counts, 2) - 1);

% Calculate p-value from the chi-squared distribution
p_value = 1 - chi2cdf(chi2_stat, df); 



        % Data for the bar graph
data = [Percentage_Passed_Bins_PriorOnly, Percentage_Passed_Bins_PretoneOnly];
labels = {'LED', 'Pretone'};

% Create the bar graph with specific figure dimensions
figure('Units', 'inches', 'Position', [0, 0, 7.5, 6]);

% Create the bar graph and manually set colors
b = bar(data);
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1]; % Set LED Trials (first bar) to blue
b.CData(2,:) = [1 0 0]; % Set Pretone Trials (second bar) to red

% Set x-tick labels and font size
set(gca, 'XTickLabel', labels, 'FontSize', 24, 'FontName', 'Arial'); 

% Y-axis label and limits
ylabel('Percent Sig. CrossCorr', 'FontSize', 24, 'FontName', 'Arial');
ylim([0 100]);
yticks(0:10:100);

% Adjust Y-axis and X-axis font sizes
ax = gca;
ax.YAxis.FontSize = 24;
ax.YAxis.FontName = 'Arial'; % Set Y-axis font to Arial
ax.XAxis.FontSize = 24;
ax.XAxis.FontName = 'Arial'; % Set X-axis font to Arial

% Check p-value and add text to the top right corner of the plot
if p_value <= 0.05
    textStr = 'p-value < 0.05';
else
    textStr = 'p-value > 0.05';
end

% Position text in the top right corner (northeast)
xPos = 3; % Adjust this based on your x-axis range
yPos = 95; % A suitable y-value for positioning in the plot
text('String', textStr, ...
     'Position', [xPos, yPos], 'FontSize', 18, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontName', 'Arial');

% Save the figure at high resolution
save_dir = 'C:\Users\Corey Roach\Desktop\2024_SFN_Poster'; % Replace with your desired directory
save_file_name = fullfile(save_dir, 'CrossCorr_Percentages.png'); % Full path and file name

% Save as PNG at 300 DPI
print(save_file_name, '-dpng', '-r300'); % '-dpng' for PNG format, '-r300' for 300 DPI resolution

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