%% Define data

Frequency_Band  = 'theta';               % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'PSI';                 % 'Granger' or 'Coherence'
animals         = {'MrCassius'};               % 'MrCassius' and/or 'MrM'
Epoch           = 'testToneOnset';

rootdir  = 'C:\Users\Corey Roach\Documents\00_DATA\PretoneOnly_TestTone_HighSNR_Correct_Congruency_PSI';
sessions = dir(fullfile(rootdir, '19*'));

if strcmp(Statistic, 'Coherence')
    ChanPair_Array_PSI    = zeros(20,20);   % used to store mean PSI collapsed across all sessions
end

%%

% Define the PFC and AC channel labels and generate temporary grids

channels = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Condition_1_PSI_grid = cell(20, 20); 
Condition_2_PSI_grid = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

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
           
            dataMatrix_1 = PSI_Condition_1.psispctrm; % Adjust based on your file structure
            labels_1 = PSI_Condition_1.labelcmb;         % Adjust based on your file structure
            dataMatrix_2 = PSI_Condition_2.psispctrm; % Adjust based on your file structure
            labels_2 = PSI_Condition_2.labelcmb;         % Adjust based on your file structure
    
           % Populate Condition_1_PSI_grid with data
        
           for rowIdx_1 = 1:size(dataMatrix_1, 1)

            % Extract channel pair from label
            pfcLabel = labels_1{rowIdx_1, 1}; % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel = labels_1{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
            
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
        
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
            acLabel = labels_2{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
            
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
        
            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx = find(strcmp(channels, strrep(acChannel, 'AC_', '')));
        
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
          end
      end
  end

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

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 8:14;

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

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

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

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 8:14;

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

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

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
%%

% Plot the heatmap
figure;
imagesc(averagedGrid_1);

% Customize the color map and add color bar
colormap(jet); % Use the 'jet' colormap (adjust as preferred)
colorbar; % Add a color bar to indicate the scale
caxis([min(averagedGrid_1(:)), max(averagedGrid_1(:))]); % Scale color axis to data range


% Add axis labels and title
xlabel('AC Channels');
ylabel('PFC Channels');
title([ 'Congruent  ', '', Statistic, '-', Frequency_Band, '-', Animal]);

% Adjust axis ticks to reflect channel numbers
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));

% Improve visual appearance
axis square; % Make the heatmap square

%%
% Plot the heatmap
figure;
imagesc(averagedGrid_2);

% Customize the color map and add color bar
colormap(jet); % Use the 'jet' colormap (adjust as preferred)
colorbar; % Add a color bar to indicate the scale
caxis([min(averagedGrid_2(:)), max(averagedGrid_2(:))]); % Scale color axis to data range

% Add axis labels and title
xlabel('AC Channels');
ylabel('PFC Channels');
title([ 'Incongruent  ', '', Statistic, '-', Frequency_Band, '-', Animal]);

% Adjust axis ticks to reflect channel numbers
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));

% Improve visual appearance
axis square; % Make the heatmap square

% %%
% 
% % Determine the global color range
% globalMin = min([averagedGrid_1(:); averagedGrid_2(:)]);
% globalMax = max([averagedGrid_1(:); averagedGrid_2(:)]);
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

%%

% Define a function to replace 0 or negative values with NaN
replaceNegativesWithNaN = @(matrix) matrix .* (matrix > 0) + NaN .* (matrix <= 0);

% Initialize a new cell array for the modified data
PriorOnly_Congruent_NaN = cell(size(PriorOnly_Congruent));

% Loop through each cell in the original 20x20 cell array
for i = 1:size(PriorOnly_Congruent, 1) % Loop over rows
    for j = 1:size(PriorOnly_Congruent, 2) % Loop over columns
        % Apply the function and store the result in the new cell array
        PriorOnly_Congruent_NaN{i, j} = replaceNegativesWithNaN(PriorOnly_Congruent{i, j});
    end
end


