%% Step 3: Organize PSI values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals    = {'MrCassius'};                     % Options: 'MrCassius' and/or 'MrM'
rootdir    = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis';                                 
sessions   = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
preCueOnset_PSI_grid   = cell(20, 20); 
testToneOnset_PSI_grid = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            preCue_file   = dir(fullfile(session_path, ['*', 'preCueOnset_' Frequency_Band '*PSI_data.mat']));
            testTone_file = dir(fullfile(session_path, ['*', 'testToneOnset_' Frequency_Band '*PSI_data.mat']));

            if ~isempty(preCue_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(preCue_file)

        preCue_filename   = preCue_file(k).name;
        testTone_filename = testTone_file(k).name;

        temp_preCue_dir   = fullfile(session_path, preCue_filename);
        temp_testTone_dir = fullfile(session_path, testTone_filename);

        % Load the .mat file safely
         load(temp_preCue_dir,'PSI_preCueOnset' );
         load(temp_testTone_dir,'PSI_testToneOnset');

        preCue_Matrix     = PSI_preCueOnset.psispctrm; 
        preCue_labels     = PSI_preCueOnset.labelcmb;         
        testTone_Matrix   = PSI_testToneOnset.psispctrm; 
        testTone_labels   = PSI_testToneOnset.labelcmb;       

      end    

           % Populate preCueOnset_PSI_grid with data 

           for rowIdx_1 = 1:size(preCue_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = preCue_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = preCue_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    preCueOnset_PSI_grid{pfcIdx, acIdx} = preCue_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    preCueOnset_PSI_grid{pfcIdx, acIdx} = [preCueOnset_PSI_grid{pfcIdx, acIdx}; preCue_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate testToneOnset_PSI_grid with data 

           for rowIdx_2 = 1:size(testTone_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = testTone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = testTone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    testToneOnset_PSI_grid{pfcIdx, acIdx} = testTone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    testToneOnset_PSI_grid{pfcIdx, acIdx} = [testToneOnset_PSI_grid{pfcIdx, acIdx}; testTone_Matrix(rowIdx_2, :)];
                    end                 
                end
           end

         end
     end
end

%% Average selected bands to output a single PSI value per channel pair (preCueOnset)


if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
preCue_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = preCueOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            preCue_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;

% Initialize a new grid to store the averaged values
preCue_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = preCueOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            preCue_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
preCue_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = preCueOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            preCue_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end


%% Average selected bands to output a single PSI value per channel pair (testToneOnset)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
testTone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = testToneOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            testTone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;

% Initialize a new grid to store the averaged values
testTone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = testToneOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            testTone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
testTone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = testToneOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            testTone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end

%% Plot average PSI 

% Replace NaNs with zero
preCue_averagedGrid(isnan(preCue_averagedGrid)) = 0;
testTone_averagedGrid(isnan(testTone_averagedGrid)) = 0;

% Determine the global color range
globalMin = -0.02;
globalMax = 0.02;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, preCue_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, testTone_averagedGrid));

% Plot the first heatmap
figure(1);
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
caxis([globalMin, globalMax]);

% Colorbar formatting
cb1 = colorbar;
ylabel(cb1, 'PSI Value', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
cb1.Label.Rotation = 270;

% Axes labels
xlabel('AC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('PFC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');

% Tick marks and labels (simplified 1:20)
xticks(1:20);
yticks(1:20);
xticklabels(1:20);
yticklabels(1:20);

% Set font for axes ticks
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

axis square;

% Plot the second heatmap
subplot(1, 2, 2);
imagesc(averagedGrid_2_clamped);
colormap(jet(256));
caxis([globalMin, globalMax]);

% Colorbar formatting
cb2 = colorbar;
ylabel(cb2, 'PSI Value', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
cb2.Label.Rotation = 270;

% Axes labels
xlabel('AC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
ylabel('PFC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');

% Tick marks and labels (simplified 1:20)
xticks(1:20);
yticks(1:20);
xticklabels(1:20);
yticklabels(1:20);

% Set font for axes ticks
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

axis square;


%% Statistics


% Initialize an empty array to store all trials
preCue_all_trials = [];

% Iterate through each cell in preCueOnset_PSI_grid
for i = 1:numel(preCueOnset_PSI_grid)
    % Concatenate the trials from the current cell
    preCue_all_trials = [preCue_all_trials; preCueOnset_PSI_grid{i}];
end

% Initialize an empty array to store all trials
testTone_all_trials = [];

% Iterate through each cell in testToneOnset_PSI_grid
for b = 1:numel(testToneOnset_PSI_grid)
    % Concatenate the trials from the current cell
    testTone_all_trials = [testTone_all_trials; testToneOnset_PSI_grid{i}];
end



%% Band-wise collect all values for each channel pair, and collapse into column

if strcmp(Frequency_Band,'theta') == 1

% Extract columns 5 to 8 and reshape into a single column
preCue_theta = preCue_all_trials(:, 5:8);
preCue_theta = preCue_theta(:);

testTone_theta = testTone_all_trials(:, 5:8);
testTone_theta = testTone_theta(:);

preCue_PFC_2_AC = preCue_theta(preCue_theta > 0);
preCue_AC_2_PFC = preCue_theta(preCue_theta < 0);

testTone_PFC_2_AC = testTone_theta(testTone_theta > 0);
testTone_AC_2_PFC = testTone_theta(testTone_theta < 0);

end

%% plotting histogram 

%[~, bins] = hist([preCue_PFC_2_AC; abs(preCue_AC_2_PFC); testTone_PFC_2_AC; abs(testTone_AC_2_PFC)], 100);
load('C:\Users\Corey Roach\Documents\00_DATA\bins.mat')

% Compute histograms for each dataset
posCounts_1 = hist(preCue_PFC_2_AC, bins);
negCounts_1 = hist(abs(preCue_AC_2_PFC), bins);
posCounts_2 = hist(testTone_PFC_2_AC, bins);
negCounts_2 = hist(abs(testTone_AC_2_PFC), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure(2);

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% xlim([-.15 .15])
axis([-.03 .03 -100 100]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

figure(3)

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.03 .03 -100 100]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

% 
% %% Statistics
% 
% % Assuming preCue and testTone are both 20x20 cell arrays
% % and each cell contains a struct with numeric fields
% 
% 
% subtracted_PSI_grid = cell(size(preCueOnset_PSI_grid)); % Preallocate
% 
% for i = 1:size(preCueOnset_PSI_grid, 1)
%     for j = 1:size(preCueOnset_PSI_grid, 2)
%         subtracted_PSI_grid{i, j} = testToneOnset_PSI_grid{i, j} - preCueOnset_PSI_grid{i, j};
%     end
% end
% 
% 
% subtracted_all_trials = [];
% 
% % Iterate through each cell in preCueOnset_PSI_grid
% for i = 1:numel(preCueOnset_PSI_grid)
%     % Concatenate the trials from the current cell
%    subtracted_all_trials = [subtracted_all_trials; subtracted_PSI_grid{i}];
% end
% 
% % Extract columns 5 to 8 and reshape into a single column
% diff_theta = subtracted_all_trials(:, 5:8);
% diff_theta = diff_theta(:);
% 
% diff_PFC_2_AC = preCue_theta(preCue_theta > 0);
% diff_AC_2_PFC = preCue_theta(preCue_theta < 0);



% %% Plot average PSI 
% 
% % Replace NaNs with zero
% preCue_averagedGrid(isnan(preCue_averagedGrid)) = 0;
% testTone_averagedGrid(isnan(testTone_averagedGrid)) = 0;
% 
% % Determine the global color range
% globalMin = -0.015;
% globalMax = 0.015;
% 
% % Clamp values to the global color range
% averagedGrid_1_clamped = max(globalMin, min(globalMax, preCue_averagedGrid));
% averagedGrid_2_clamped = max(globalMin, min(globalMax, testTone_averagedGrid));
% 
% % Plot the first heatmap
% figure(1);
% subplot(1, 2, 1);
% imagesc(averagedGrid_1_clamped);
% colormap(jet(256)); % 256 instead of default 64
% caxis([globalMin, globalMax]);
% 
% % Colorbar formatting
% cb1 = colorbar;
% ylabel(cb1, 'PSI Value', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% cb1.Label.Rotation = 270;
% 
% % Axes labels
% xlabel('AC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% ylabel('PFC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% 
% % Tick marks and labels (simplified 1:20)
% xticks(1:20);
% yticks(1:20);
% xticklabels(1:20);
% yticklabels(1:20);
% 
% % Set font for axes ticks
% set(gca, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
% 
% axis square;
% 
% % Plot the second heatmap
% subplot(1, 2, 2);
% imagesc(averagedGrid_2_clamped);
% colormap(jet(256));
% caxis([globalMin, globalMax]);
% 
% % Colorbar formatting
% cb2 = colorbar;
% ylabel(cb2, 'PSI Value', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% cb2.Label.Rotation = 270;
% 
% % Axes labels
% xlabel('AC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% ylabel('PFC Channels', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% 
% % Tick marks and labels (simplified 1:20)
% xticks(1:20);
% yticks(1:20);
% xticklabels(1:20);
% yticklabels(1:20);
% 
% % Set font for axes ticks
% set(gca, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
% 
% axis square;
