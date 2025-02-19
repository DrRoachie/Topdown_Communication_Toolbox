%% Define data

Frequency_Band  = 'theta';                       % Options: 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'PSI';                         % Options: 'Granger' or 'Coherence'
animals         = {'MrCassius'};          % Options: 'MrCassius' and/or 'MrM'
Epoch           = 'preCueOnset';
Behavior        = 'Correct';                     % Options: Correct or Wrong

% Select PSI data directory
% as of 1/14/25 the most important input is the directory because it
% determines what two conditions are stored in the file 

rootdir  = 'C:\Users\Corey Roach\Documents\00_DATA\PretoneOnly_PreCue_Correct_Congruency_PSI';
%savedir = 'C:\Users\Corey Roach\Documents\00_DATA\PretoneOnly_TestTone_HighSNR_Correct_Congruency_PSI';                                        % path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
sessions = dir(fullfile(rootdir, '19*'));

% designations are used for outfile naming convention only
SNR_designation = 'Both';                     % Options: HighSNR or LowSNR
Condition_1_designation = 'Congruent_Trials';    % Options: Can be almost any 
Condition_2_designation = 'Incongruent_Trials';

%%

% Define the PFC and AC channel labels and generate temporary grids

channels = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Condition_1_PSI_grid = cell(20, 20); 
Condition_2_PSI_grid = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

% Preallocate cell arrays for storing results
num_sessions = length(sessions);
monte_Results = cell(num_sessions, 1);
Wilcox_Results = cell(num_sessions, 1);

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        
        if exist(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct'), 'dir')
            if strcmp(Statistic, 'PSI')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', ['*' Frequency_Band '*PSI_data.mat']));
            end
        else
        continue;
        end
        
        if ~isempty(files)
            fprintf('Processing session %s: ', num2str(RecDate));
        end

            for k = 1:length(files) % runs if there are any processed channel pairs in given frequency band
        
            file_name = files(k).name;

            % pull this sessions, PSI data 
            if strcmp(Statistic, 'PSI')
            channel_files = dir(fullfile(rootdir, RecDate, Animal, 'testTone', 'Correct', [extractBefore(file_name, 'PSI') '*']));
            end

            % 
            datadir = fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', file_name);
            load(datadir);
           
            dataMatrix_1 = PSI_Condition_1.psispctrm; 
            labels_1 = PSI_Condition_1.labelcmb;         
            dataMatrix_2 = PSI_Condition_2.psispctrm; 
            labels_2 = PSI_Condition_2.labelcmb;         
    
           % Populate Condition_1_PSI_grid with data
        
           for rowIdx_1 = 1:size(dataMatrix_1, 1)

            % Extract channel pair from label
            pfcLabel = labels_1{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = labels_1{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
            
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
        
            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx = find(strcmp(channels, strrep(acChannel, 'AC_', '')));
        
                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Condition_1_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Condition_1_PSI_grid{pfcIdx, acIdx} = dataMatrix_1(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Condition_1_PSI_grid{pfcIdx, acIdx} = [Condition_1_PSI_grid{pfcIdx, acIdx}; dataMatrix_1(rowIdx_1, :)];
                    end

                end
           end

           % Populate Condition_2_PSI_grid with data 

           for rowIdx_2 = 1:size(dataMatrix_2, 1)

            % Extract channel pair from label
            pfcLabel = labels_2{rowIdx_2, 1}; % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = labels_2{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
            
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
        
            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));
        
                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Condition_2_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Condition_2_PSI_grid{pfcIdx, acIdx} = dataMatrix_2(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Condition_2_PSI_grid{pfcIdx, acIdx} = [Condition_2_PSI_grid{pfcIdx, acIdx}; dataMatrix_2(rowIdx_2, :)];
                    end                 
                end
           end
          
            
           [Wilcox_results, monte_results] = PSIScratchPad_v2(Frequency_Band, PSI_Condition_1, PSI_Condition_2);

            % Store results for this session
            monte_Results{i} = monte_results;
            Wilcox_Results{i} = Wilcox_results;


          end
      end
  end


%% 

num_sessions = length(Wilcox_Results);

% Preallocate column arrays
PFC_to_AC_Monte = [];
AC_to_PFC_Monte = [];

for i = 1:num_sessions
    % Check if the current cell in monte_Results is empty
    if isempty(monte_Results{i})
        continue; % Skip this iteration if the cell is empty
    end
    
    % Check if the required fields exist before accessing them
    if isfield(monte_Results{i}, 'pfc_to_ac') && isfield(monte_Results{i}, 'ac_to_pfc')
        % Extract PFC_to_AC and AC_to_PFC values safely
        PFC_to_AC_Monte = [PFC_to_AC_Monte; monte_Results{i}.pfc_to_ac(:)];
        AC_to_PFC_Monte = [AC_to_PFC_Monte; monte_Results{i}.ac_to_pfc(:)];
    end
end


%% percentage

% Count the number of values ≤ 0.05
PFC_to_AC_num_significant = sum(PFC_to_AC_Monte <= 0.05);
PFC_to_AC_total_values = numel(PFC_to_AC_Monte);
PFC_to_AC_percent_significant = (PFC_to_AC_num_significant / PFC_to_AC_total_values) * 100;

AC_to_PFC_num_significant = sum(AC_to_PFC_Monte <= 0.05);
AC_to_PFC_total_values = numel(AC_to_PFC_Monte);
AC_to_PFC_percent_significant = (AC_to_PFC_num_significant / AC_to_PFC_total_values) * 100;

%% Chi-Squared Test for Proportions 

observed = [PFC_to_AC_num_significant, PFC_to_AC_total_values - PFC_to_AC_num_significant; AC_to_PFC_num_significant, AC_to_PFC_total_values - AC_to_PFC_num_significant];

% Perform chi-square test
[h, p, stats] = fishertest(observed);


%% Average selected bands to output a single PSI value per channel pair (Condition 1)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
averagedGrid_1 = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Condition_1_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = Condition_1_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            averagedGrid_1(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end
% 
% if strcmp(Frequency_Band,'alpha') == 1
% % Define the range of columns to average
% columnRange = 8:14;
% 
% % Initialize a new grid to store the averaged values
% averagedGrid_1 = nan(20, 20); % Use NaN to handle empty grid positions
% 
% % Loop through each position in the grid
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         % Check if the grid position contains data
%         if ~isempty(Condition_1_PSI_grid{pfcIdx, acIdx})
%             % Extract the data matrix for the current position
%             dataposition = Condition_1_PSI_grid{pfcIdx, acIdx};
% 
%             % Compute the average of the selected columns
%             columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');
% 
%             % Store the average in the new grid
%             averagedGrid_1(pfcIdx, acIdx) = columnAverage;
%         end
%     end
% end
% end
% 
% if strcmp(Frequency_Band,'beta') == 1
% % Define the range of columns to average
% columnRange = 15:30;
% 
% % Initialize a new grid to store the averaged values
% averagedGrid_1 = nan(20, 20); % Use NaN to handle empty grid positions
% 
% % Loop through each position in the grid
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         % Check if the grid position contains data
%         if ~isempty(Condition_1_PSI_grid{pfcIdx, acIdx})
%             % Extract the data matrix for the current position
%             dataposition = Condition_1_PSI_grid{pfcIdx, acIdx};
% 
%             % Compute the average of the selected columns
%             columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');
% 
%             % Store the average in the new grid
%             averagedGrid_1(pfcIdx, acIdx) = columnAverage;
%         end
%     end
% end
% end
% 
%% Average selected bands to output a single PSI value per channel pair (Condition 2)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
averagedGrid_2 = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Condition_2_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = Condition_2_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            averagedGrid_2(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

% if strcmp(Frequency_Band,'alpha') == 1
% % Define the range of columns to average
% columnRange = 8:14;
% 
% % Initialize a new grid to store the averaged values
% averagedGrid_2 = nan(20, 20); % Use NaN to handle empty grid positions
% 
% % Loop through each position in the grid
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         % Check if the grid position contains data
%         if ~isempty(Condition_2_PSI_grid{pfcIdx, acIdx})
%             % Extract the data matrix for the current position
%             dataposition = Condition_2_PSI_grid{pfcIdx, acIdx};
% 
%             % Compute the average of the selected columns
%             columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');
% 
%             % Store the average in the new grid
%             averagedGrid_2(pfcIdx, acIdx) = columnAverage;
%         end
%     end
% end
% end
% 
% if strcmp(Frequency_Band,'beta') == 1
% % Define the range of columns to average
% columnRange = 15:30;
% 
% % Initialize a new grid to store the averaged values
% averagedGrid_2 = nan(20, 20); % Use NaN to handle empty grid positions
% 
% % Loop through each position in the grid
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         % Check if the grid position contains data
%         if ~isempty(Condition_2_PSI_grid{pfcIdx, acIdx})
%             % Extract the data matrix for the current position
%             dataposition = Condition_2_PSI_grid{pfcIdx, acIdx};
% 
%             % Compute the average of the selected columns
%             columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');
% 
%             % Store the average in the new grid
%             averagedGrid_2(pfcIdx, acIdx) = columnAverage;
%         end
%     end
% end
% end

% %%
% 
% % Determine the global color range
% % globalMin = min([averagedGrid_1(:); averagedGrid_2(:)]);
% % globalMax = max([averagedGrid_1(:); averagedGrid_2(:)]);
% 
% % Determine the global color range
% globalMin = .001;
% globalMax = .15;
% 
% % Plot the first heatmap
% figure;
% subplot(1, 2, 1); % Create a subplot for side-by-side comparison
% imagesc(averagedGrid_1);
% colormap(jet); % Apply the same colormap
% caxis([globalMin, globalMax]); % Set consistent color limits
% % caxis([-.1, .1]); % Set consistent color limits
% 
% % Plot the first heatmap
% figure;
% subplot(1, 2, 1); % Create a subplot for side-by-side comparison
% imagesc(averagedGrid_1);
% colormap(jet); % Apply the same colormap
% caxis([globalMin, globalMax]); % Set consistent color limits
% % caxis([-.1, .1]); % Set consistent color limits
% 
% colorbar; % Add color bar
% 
% % Add axis labels and title
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title([ 'Congruent  ', '', Statistic, '-', Frequency_Band, '-', Animal]);
% 
% % Adjust axis ticks to reflect channel numbers
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% 
% % Improve visual appearance
% axis square; % Make the heatmap square
% 
% % Plot the second heatmap
% subplot(1, 2, 2); % Second subplot
% imagesc(averagedGrid_2);
% colormap(jet); % Apply the same colormap
% caxis([globalMin, globalMax]); % Set consistent color limits
% % caxis([-.1, .1]); % Set consistent color limits
% colorbar; % Add color bar
% 
% % Add axis labels and title
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title([ 'Incongruent  ', '', Statistic, '-', Frequency_Band, '-', Animal]);
% 
% % Adjust axis ticks to reflect channel numbers
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% 
% % Improve visual appearance
% axis square; % Make the heatmap square

% Determine the global color range
globalMin = -0.03;
globalMax = 0.03;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, averagedGrid_1));
averagedGrid_2_clamped = max(globalMin, min(globalMax, averagedGrid_2));

% Plot the first heatmap
figure;
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
%colormap(jet)
caxis([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['Congruent  ', Statistic, '-', Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

% Plot the second heatmap
subplot(1, 2, 2);
imagesc(averagedGrid_2_clamped);
colormap(jet(256)); % 256 instead of default 64
caxis([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['Incongruent  ', Statistic, '-', Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

%%
% Assume averagedGrid_1 is your 20x20 array

% Create matrices for positive and negative values
averagedGrid_1_pos = zeros(size(averagedGrid_1));
averagedGrid_1_neg = zeros(size(averagedGrid_1));

% Assign positive values
averagedGrid_1_pos(averagedGrid_1 > 0) = averagedGrid_1(averagedGrid_1 > 0);
% Assign negative values
averagedGrid_1_neg(averagedGrid_1 < 0) = averagedGrid_1(averagedGrid_1 < 0);

% Assume averagedGrid_2 is your 20x20 array

% Create matrices for positive and negative values
averagedGrid_2_pos = zeros(size(averagedGrid_2));
averagedGrid_2_neg = zeros(size(averagedGrid_2));

% Assign positive values
averagedGrid_2_pos(averagedGrid_2 > 0) = averagedGrid_2(averagedGrid_2 > 0);
% Assign negative values
averagedGrid_2_neg(averagedGrid_2 < 0) = averagedGrid_2(averagedGrid_2 < 0);


%%

Summary = PosNegHistogram(averagedGrid_1_pos, averagedGrid_1_neg, averagedGrid_2_pos, averagedGrid_2_neg, 20);



