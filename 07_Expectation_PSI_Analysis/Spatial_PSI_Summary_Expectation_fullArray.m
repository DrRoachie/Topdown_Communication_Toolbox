%% Step 3: Organize PSI values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'beta';
Behavior          = 'correct';
Congruency        = 'congruent';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius'};                                                   % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Connectivity_FullArray_noshare_testTone_beta_Cassius';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 1:24, 'UniformOutput', false);
OnlyPrior_PSI_grid     = cell(24, 24); 
OnlyPretone_PSI_grid   = cell(24, 24);

%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));
            PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_','OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));

            if ~isempty(PSI_OnlyPrior_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(PSI_OnlyPrior_file)

        PSI_OnlyPrior_filename   = PSI_OnlyPrior_file(k).name;
        PSI_OnlyPretone_filename = PSI_OnlyPretone_file(k).name;

        temp_OnlyPrior_dir   = fullfile(session_path, PSI_OnlyPrior_filename);
        temp_OnlyPretone_dir = fullfile(session_path, PSI_OnlyPretone_filename);

        % Load the .mat file safely
         load(temp_OnlyPrior_dir,'PSI_OnlyPrior' );
         load(temp_OnlyPretone_dir,'PSI_OnlyPretone');

       OnlyPrior_Matrix     = PSI_OnlyPrior.psispctrm; 
       OnlyPrior_labels     = PSI_OnlyPrior.labelcmb;         
       OnlyPretone_Matrix   = PSI_OnlyPretone.psispctrm; 
       OnlyPretone_labels   = PSI_OnlyPretone.labelcmb;       

      end    

           % Populate OnlyPrior_PSI_grid with data 

           for rowIdx_1 = 1:size(OnlyPrior_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(OnlyPrior_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    OnlyPrior_PSI_grid{pfcIdx, acIdx} = OnlyPrior_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    OnlyPrior_PSI_grid{pfcIdx, acIdx} = [OnlyPrior_PSI_grid{pfcIdx, acIdx}; OnlyPrior_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate OnlyPretone_PSI_grid with data 

           for rowIdx_2 = 1:size(OnlyPretone_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = OnlyPretone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = OnlyPretone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(OnlyPretone_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    OnlyPretone_PSI_grid{pfcIdx, acIdx} = OnlyPretone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    OnlyPretone_PSI_grid{pfcIdx, acIdx} = [OnlyPretone_PSI_grid{pfcIdx, acIdx}; OnlyPretone_Matrix(rowIdx_2, :)];
                    end                 
                end
           end

         end
     end
end

%% Average selected bands to output a single PSI value per channel pair (OnlyPrior)


if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
OnlyPrior_averagedGrid = nan(24, 24); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:24
    for acIdx = 1:24
        % Check if the grid position contains data
        if ~isempty(OnlyPrior_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPrior_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPrior_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'alpha') == 1

% Define the range of columns to average
columnRange = 9:14;

% Initialize a new grid to store the averaged values
OnlyPrior_averagedGrid = nan(24, 24); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:24
    for acIdx = 1:24
        % Check if the grid position contains data
        if ~isempty(OnlyPrior_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPrior_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPrior_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
OnlyPrior_averagedGrid = nan(24, 24); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:24
    for acIdx = 1:24
        % Check if the grid position contains data
        if ~isempty(OnlyPrior_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPrior_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPrior_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end


%% Average selected bands to output a single PSI value per channel pair (OnlyPretone)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
OnlyPretone_averagedGrid = nan(24, 24); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:24
    for acIdx = 1:24
        % Check if the grid position contains data
        if ~isempty(OnlyPretone_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPretone_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPretone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;

% Initialize a new grid to store the averaged values
OnlyPretone_averagedGrid = nan(24, 24); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:24
    for acIdx = 1:24
        % Check if the grid position contains data
        if ~isempty(OnlyPretone_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPretone_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPretone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
OnlyPretone_averagedGrid = nan(24, 24); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:24
    for acIdx = 1:24
        % Check if the grid position contains data
        if ~isempty(OnlyPretone_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPretone_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPretone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

%

%% Plot average PSI 

% Replace NaNs with zero
OnlyPrior_averagedGrid(isnan(OnlyPrior_averagedGrid)) = 0;
OnlyPretone_averagedGrid(isnan(OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
% globalMin = -0.1;
% globalMax = 0.1;

% Clamp values to the global color range
% averagedGrid_1_clamped = max(globalMin, min(globalMax, OnlyPrior_averagedGrid));
% averagedGrid_2_clamped = max(globalMin, min(globalMax, OnlyPretone_averagedGrid));

averagedGrid_1_clamped = OnlyPrior_averagedGrid;
averagedGrid_2_clamped = OnlyPretone_averagedGrid;

% Plot the first heatmap
figure(1);
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
%clim([globalMin, globalMax]);
clim([-.05, .05]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPrior-',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 1:24, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 1:24, 'UniformOutput', false));
axis square;

% Plot the second heatmap
subplot(1, 2, 2);
imagesc(averagedGrid_2_clamped);
colormap(jet(256)); % 256 instead of default 64
%clim([globalMin, globalMax]);
clim([-.05, .05]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPretone-',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 1:24, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 1:24, 'UniformOutput', false));
axis square;

% %% Statistics
% 
% % Initialize an empty array to store all trials
% OnlyPrior_all_trials = [];
% 
% % Iterate through each cell in preCueOnset_PSI_grid
% for i = 1:numel(OnlyPrior_PSI_grid)
%     % Concatenate the trials from the current cell
%     OnlyPrior_all_trials = [OnlyPrior_all_trials; OnlyPrior_PSI_grid{i}];
% end
% 
% % Initialize an empty array to store all trials
% OnlyPretone_all_trials = [];
% 
% % Iterate through each cell in testToneOnset_PSI_grid
% for b = 1:numel(OnlyPretone_PSI_grid)
%     % Concatenate the trials from the current cell
%     OnlyPretone_all_trials = [OnlyPretone_all_trials; OnlyPretone_PSI_grid{b}];
% end
% 
% 
% %% Band-wise collect all values for each channel pair, and collapse into column
% 
% if strcmp(Frequency_Band,'theta') == 1
% 
% % Extract columns 5 to 8 and reshape into a single column
% OnlyPrior_theta = OnlyPrior_all_trials(:, 5:8);
% OnlyPrior_theta = OnlyPrior_theta(:);
% 
% OnlyPretone_theta = OnlyPretone_all_trials(:, 5:8);
% OnlyPretone_theta = OnlyPretone_theta(:);
% 
% OnlyPrior_PFC_2_AC = OnlyPrior_theta(OnlyPrior_theta > 0);
% OnlyPrior_AC_2_PFC = OnlyPrior_theta(OnlyPrior_theta < 0);
% 
% OnlyPretone_PFC_2_AC = OnlyPretone_theta(OnlyPretone_theta > 0);
% OnlyPretone_AC_2_PFC = OnlyPretone_theta(OnlyPretone_theta < 0);
% 
% end
% 
% 
% if strcmp(Frequency_Band,'alpha') == 1
% 
% % Extract columns 5 to 8 and reshape into a single column
% OnlyPrior_alpha = preCue_all_trials(:, 9:14);
% OnlyPrior_alpha = OnlyPrior_alpha(:);
% 
% OnlyPretone_alpha = OnlyPretone_all_trials(:, 9:14);
% OnlyPretone_alpha = OnlyPretone_alpha(:);
% 
% OnlyPrior_PFC_2_AC = OnlyPrior_alpha(OnlyPrior_alpha > 0);
% OnlyPrior_AC_2_PFC = OnlyPrior_alpha(OnlyPrior_alpha < 0);
% 
% OnlyPretone_PFC_2_AC = OnlyPretone_alpha(OnlyPretone_alpha > 0);
% OnlyPretone_AC_2_PFC = OnlyPretone_alpha(OnlyPretone_alpha < 0);
% 
% end
% 
% if strcmp(Frequency_Band,'beta') == 1
% 
% % Extract columns 5 to 8 and reshape into a single column
% OnlyPrior_alpha = preCue_all_trials(:, 15:30);
% OnlyPrior_alpha = OnlyPrior_alpha(:);
% 
% OnlyPretone_alpha = OnlyPretone_all_trials(:, 15:30);
% OnlyPretone_alpha = OnlyPretone_alpha(:);
% 
% OnlyPrior_PFC_2_AC = OnlyPrior_alpha(OnlyPrior_alpha > 0);
% OnlyPrior_AC_2_PFC = OnlyPrior_alpha(OnlyPrior_alpha < 0);
% 
% OnlyPretone_PFC_2_AC = OnlyPretone_alpha(OnlyPretone_alpha > 0);
% OnlyPretone_AC_2_PFC = OnlyPretone_alpha(OnlyPretone_alpha < 0);
% 
% end
% 
% %%
% 
% % Combine the two arrays
% OnlyPrior_combined_array   = [OnlyPrior_PFC_2_AC; OnlyPrior_AC_2_PFC];
% OnlyPrior_combined_array   = abs(OnlyPrior_combined_array);
% 
% OnlyPretone_combined_array = [OnlyPretone_PFC_2_AC; OnlyPretone_AC_2_PFC];
% OnlyPretone_combined_array = abs(OnlyPretone_combined_array);
% 
% %%
% 
% % save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Belt_PSI_Values.mat'), 'OnlyPrior_combined_array', 'OnlyPretone_combined_array');
% 
% %%
% % Perform the Mann-Whitney U test (Wilcoxon rank-sum test) twice
% [p, h, stats] = ranksum(OnlyPrior_combined_array, OnlyPretone_combined_array);  
% 
% % Bonferroni correction: adjust p-value threshold for two tests
% alpha = 0.05;
% adjusted_alpha = alpha / 4;
% 
% % Display adjusted p-values and significance results
% % fprintf('Adjusted p-value for test 1: %.4f\n', PFC_2_AC_p);
% % fprintf('Adjusted p-value for test 2: %.4f\n', AC_2_PFC_p);
% % fprintf('Test 1 significant (Bonferroni): %d\n', PFC_2_AC_p < adjusted_alpha);
% % fprintf('Test 2 significant (Bonferroni): %d\n', AC_2_PFC_p < adjusted_alpha);
% 
% % Calculate means and 95% CI for both data sets
% mean1 = mean(OnlyPrior_combined_array);
% mean2 = mean(OnlyPretone_combined_array);
% stderr1 = std(OnlyPrior_combined_array) / sqrt(length(OnlyPrior_combined_array));  % Standard error of the mean
% stderr2 = std(OnlyPretone_combined_array) / sqrt(length(OnlyPretone_combined_array));  % Standard error of the mean
% CI1 = 1.96 * stderr1;  % 95% confidence interval for data1
% CI2 = 1.96 * stderr2;  % 95% confidence interval for data2
% 
% % Plot grouped bars
% hold on 
% bar(1, mean1, 'FaceColor', [0.4, 0.8, 1.0], 'EdgeColor', 'k', 'BarWidth', .75);
% bar(2, mean2, 'FaceColor', [1.0, 0.7, 0.3], 'EdgeColor', 'k', 'BarWidth', .75);
% 
% % Add error bars for 95% CI
% errorbar(1, mean1, CI1, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% errorbar(2, mean2, CI2, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% 
% % Customize the plot
% set(gca, 'XTick', [1 2], 'XTickLabel', {'LED', 'Pretone'});  % Set x-axis labels to 'LED' and 'Pretone'
% ylim([0 .018])
% ylabel('Mean PSI Value');
% %title('Mean PSI with 95% Confidence Interval');
% 
% % Adjust the axis to reduce space between bars
% set(gca, 'XLim', [0.5 2.5]);  % Adjust the x-axis limits to reduce space between bars
% 
% hold off;

% %% plotting histogram 
% 
% %[~, bins] = hist([OnlyPrior_PFC_2_AC; abs(OnlyPrior_AC_2_PFC); OnlyPretone_PFC_2_AC; abs(OnlyPretone_AC_2_PFC)], 100);
% load('C:\Users\Corey Roach\Documents\00_DATA\bins.mat')
% 
% % Compute histograms for each dataset
% posCounts_1 = hist(OnlyPrior_PFC_2_AC, bins);
% negCounts_1 = hist(abs(OnlyPrior_AC_2_PFC), bins);
% posCounts_2 = hist(OnlyPretone_PFC_2_AC, bins);
% negCounts_2 = hist(abs(OnlyPretone_AC_2_PFC), bins);
% 
% % Normalize counts to percentages
% posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
% negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
% posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
% negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;
% 
% % Plot histograms
% 
% figure(2);
% 
% hold on;
% % Condition 1: Red (Positive) & Soft Red (Negative)
% bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
% bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
% 
% 
% % Labels and Legend
% xlabel('PSI Value');
% ylabel('Proportion of Sites');
% % ylim([-65 65])
% % xlim([-.15 .15])
% axis([-.025 .025 -65 65]) 
% leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
% set(leg, 'FontSize', 12); % Reduce legend font size
% 
% % Improve aesthetics
% set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
% box on;
% hold off;

% figure(3)
% 
% hold on 
% bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
% bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
% 
% % Labels and Legend
% xlabel('PSI Value', 'FontSize', 18);
% ylabel('Proportion of Sites', 'FontSize', 18);
% axis([-.025 .025 -65 65]) 
% 
% % Set legend and adjust its font size
% leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
% set(leg, 'FontSize', 12); % Reduce legend font size
% 
% % Improve aesthetics
% set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
% box on;
% hold off;
% 
% %'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
% %FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges;

% %% Statistics
% 
% % Assuming preCue and testTone are both 20x20 cell arrays
% % and each cell contains a struct with numeric fields
% 
% 
% subtracted_PSI_grid = cell(size(OnlyPrior_PSI_grid)); % Preallocate
% 
% for i = 1:size(OnlyPrior_PSI_grid, 1)
%     for j = 1:size(OnlyPrior_PSI_grid, 2)
%         subtracted_PSI_grid{i, j} = OnlyPretone_PSI_grid{i, j} - OnlyPrior_PSI_grid{i, j};
%     end
% end
% 
% 
% subtracted_all_trials = [];
% 
% % Iterate through each cell in preCueOnset_PSI_grid
% for i = 1:numel(OnlyPrior_PSI_grid)
%     % Concatenate the trials from the current cell
%    subtracted_all_trials = [subtracted_all_trials; subtracted_PSI_grid{i}];
% end
% 
% % Extract columns 5 to 8 and reshape into a single column
% diff_theta = subtracted_all_trials(:, 5:8);
% diff_theta = diff_theta(:);
% 
% diff_PFC_2_AC = OnlyPrior_theta(OnlyPrior_theta > 0);
% diff_AC_2_PFC = OnlyPrior_theta(OnlyPrior_theta < 0);

%% Jamie Kuminski Core/Belt Analysis

% Core_OnlyPrior_averagedGrid   = OnlyPrior_averagedGrid;
% Core_OnlyPretone_averagedGrid = OnlyPretone_averagedGrid;
% save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Core_Spatial_PSI_Values.mat'), 'Core_OnlyPrior_averagedGrid', 'Core_OnlyPretone_averagedGrid');

% Belt_OnlyPrior_averagedGrid   = OnlyPrior_averagedGrid;
% Belt_OnlyPretone_averagedGrid = OnlyPretone_averagedGrid;
% save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Belt_Spatial_PSI_Values.mat'), 'Belt_OnlyPrior_averagedGrid', 'Belt_OnlyPretone_averagedGrid');

% Subtract the session-wise means from each other and plot as a histogram.

% % Subtract the conditions
% Belt_Subtracted = Belt_OnlyPrior_averagedGrid - Belt_OnlyPretone_averagedGrid;
% Core_Subtracted = Core_OnlyPrior_averagedGrid - Core_OnlyPretone_averagedGrid; 
% 
% % Reshape the matrices into column vectors
% Belt_Subtracted_Column = Belt_Subtracted(:);
% Core_Subtracted_Column = Core_Subtracted(:);
% 
% % Subtract the conditions
% Belt_Subtracted = Belt_OnlyPrior_averagedGrid - Belt_OnlyPretone_averagedGrid;
% Core_Subtracted = Core_OnlyPrior_averagedGrid - Core_OnlyPretone_averagedGrid; 
% 
% % Reshape the matrices into column vectors
% Belt_Subtracted_Column = Belt_Subtracted(:);
% Core_Subtracted_Column = Core_Subtracted(:);
% 
% % Set the same number of bins for both histograms
% numBins = 20;
% 
% % Calculate the edges for both histograms using the same bin range
% edges = linspace(min([Belt_Subtracted_Column; Core_Subtracted_Column]), max([Belt_Subtracted_Column; Core_Subtracted_Column]), numBins+1);
% 
% % Calculate the histogram counts for each condition using the same bin edges
% [counts_belt, ~] = histcounts(Belt_Subtracted_Column, edges);
% [counts_core, ~] = histcounts(Core_Subtracted_Column, edges);
% 
% % Plot the histograms with transparency
% hold on;
% 
% % Plot the first histogram for Belt_Subtracted with a specific color and transparency
% h1 = bar(edges(1:end-1), counts_belt, 'FaceColor', [0.4, 1.0, 0.4], 'EdgeColor', 'k', 'BarWidth', 1);
% set(h1, 'FaceAlpha', 0.5);  % Set transparency for Belt_Subtracted
% 
% % Plot the second histogram for Core_Subtracted with a different color and transparency
% h2 = bar(edges(1:end-1), counts_core, 'FaceColor', [1.0, 0.4, 0.4], 'EdgeColor', 'k', 'BarWidth', 1);
% set(h2, 'FaceAlpha', 0.5);  % Set transparency for Core_Subtracted
% 
% hold off;
% 
% % Customize the plot
% set(gca, 'XTick', linspace(min(edges), max(edges), numBins), 'XTickLabel', round(linspace(min(edges), max(edges), numBins), 2));  % Set x-axis labels
% ylabel('Frequency');
% xlabel('Subtracted mean PSI values');% 
% % Optionally, add a legend to differentiate the histograms
% legend({'Belt Subtracted', 'Core Subtracted'});



%%

% % Perform K-S test for normality on Core_Subtracted and Belt_Subtracted
% % Null Hypothesis: Data follows a normal distribution
% 
% % Test on Core_Subtracted (shifted to zero mean)
% [h_core, p_core] = kstest(Core_Subtracted_Column - mean(Core_Subtracted_Column));  % Center around zero
% fprintf('K-S Test for Core Subtracted: h = %d, p = %.4f\n', h_core, p_core);
% 
% % Test on Belt_Subtracted (shifted to zero mean)
% [h_belt, p_belt] = kstest(Belt_Subtracted_Column - mean(Belt_Subtracted_Column));  % Center around zero
% fprintf('K-S Test for Belt Subtracted: h = %d, p = %.4f\n', h_belt, p_belt);
% 
% % Perform K-S test to compare the distributions of Core and Belt Subtracted
% % Null Hypothesis: The two samples are from the same distribution
% 
% % Perform K-S test between the two distributions
% [h_ks, p_ks] = kstest2(Core_Subtracted_Column, Belt_Subtracted_Column);  % Compare distributions
% fprintf('K-S Test between Core and Belt Subtracted: h = %d, p = %.4f\n', h_ks, p_ks);
% 
% 
% % Skewness and Kurtosis for Core_Subtracted
% skew_core = skewness(Core_Subtracted_Column);
% kurtosis_core = kurtosis(Core_Subtracted_Column);
% 
% % Skewness and Kurtosis for Belt_Subtracted
% skew_belt = skewness(Belt_Subtracted_Column);
% kurtosis_belt = kurtosis(Belt_Subtracted_Column);
% 
% % Display the results
% fprintf('Core Subtracted - Skewness: %.4f, Kurtosis: %.4f\n', skew_core, kurtosis_core);
% fprintf('Belt Subtracted - Skewness: %.4f, Kurtosis: %.4f\n', skew_belt, kurtosis_belt);

