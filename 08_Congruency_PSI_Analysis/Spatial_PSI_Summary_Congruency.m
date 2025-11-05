%% Congruency Analysis
% Updated: CR 10/1/25
% This code consolidates the PSI values from two separate directories
% outputed from the 'Expectation_Connectivity_Test*.m' where one contains
% the PSI values for congruent trials and the other contains the values for
% incongruent trials. Ultimately, we end up with a 2x2x2 conditions
% [Expectation x Congruency x Direction]

% OnlyPrior~Congruent 
% OnlyPrior~Incongruent
% OnlyPretone ~ Congruent
% OnlyPretone ~ Incongruent

%% Organize PSI values by channel pair. [CONGRUENT]

Frequency_Band    = 'theta';
Behavior          = 'correct';
animals           = {'MrCassius'; 'MrM'};                                                   % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_07_10_PSI_Expectation_Analysis_testTone_Congruent';      % This is the congruent                            
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

Congruent_channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Congruent_OnlyPrior_PSI_grid     = cell(20, 20); 
Congruent_OnlyPretone_PSI_grid   = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair [CONGRUENT]

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            Congruent_PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', 'congruent', '_', Frequency_Band,'_PSI.mat']));
            Congruent_PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_','OnlyPretone_', Behavior,'_', 'congruent', '_', Frequency_Band,'_PSI.mat']));

            if ~isempty(Congruent_PSI_OnlyPrior_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(Congruent_PSI_OnlyPrior_file)

        Congruent_PSI_OnlyPrior_filename   = Congruent_PSI_OnlyPrior_file(k).name;
        Congruent_PSI_OnlyPretone_filename = Congruent_PSI_OnlyPretone_file(k).name;

        temp_Congruent_OnlyPrior_dir   = fullfile(session_path, Congruent_PSI_OnlyPrior_filename);
        temp_Congruent_OnlyPretone_dir = fullfile(session_path, Congruent_PSI_OnlyPretone_filename);

        % Load the .mat file safely
         load(temp_Congruent_OnlyPrior_dir,'PSI_OnlyPrior' );
         load(temp_Congruent_OnlyPretone_dir,'PSI_OnlyPretone');

       Congruent_OnlyPrior_Matrix     = PSI_OnlyPrior.psispctrm; 
       Congruent_OnlyPrior_labels     = PSI_OnlyPrior.labelcmb;         
       Congruent_OnlyPretone_Matrix   = PSI_OnlyPretone.psispctrm; 
       Congruent_OnlyPretone_labels   = PSI_OnlyPretone.labelcmb;       

      end    

           % Populate Congruent_OnlyPrior_PSI_grid with data 

           for rowIdx_1 = 1:size(Congruent_OnlyPrior_Matrix, 1)

            % Extract channel pair from label
            Congruent_pfcLabel = Congruent_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            Congruent_acLabel  = Congruent_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            Congruent_pfcChannel = regexp(Congruent_pfcLabel, 'PFC_ch\d+', 'match', 'once');    % Extract PFC channel
            Congruent_acChannel  = regexp(Congruent_acLabel,  'AC_ch\d+',  'match', 'once');    % Extract AC channel

            % Find indices in the grid
            Congruent_pfcIdx = find(strcmp(Congruent_channels,  strrep(Congruent_pfcChannel, 'PFC_', '')));
            Congruent_acIdx  = find(strcmp(Congruent_channels,  strrep(Congruent_acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(Congruent_pfcIdx) && ~isempty(Congruent_acIdx)
                    if isempty(Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx})
                    % Initialize the cell with the first row of data
                    Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx} = Congruent_OnlyPrior_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx} = [Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx}; Congruent_OnlyPrior_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate OnlyPretone_PSI_grid with data 

           for rowIdx_2 = 1:size(Congruent_OnlyPretone_Matrix, 1)

            % Extract channel pair from label
            Congruent_pfcLabel = Congruent_OnlyPretone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            Congruent_acLabel  = Congruent_OnlyPretone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            Congruent_pfcChannel = regexp(Congruent_pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            Congruent_acChannel  = regexp(Congruent_acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            Congruent_pfcIdx = find(strcmp(Congruent_channels, strrep(Congruent_pfcChannel, 'PFC_', '')));
            Congruent_acIdx  = find(strcmp(Congruent_channels, strrep(Congruent_acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(Congruent_pfcIdx) && ~isempty(Congruent_acIdx)
                    if isempty(Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx})
                    % Initialize the cell with the first row of data
                    Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx} = Congruent_OnlyPretone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx} = [Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx}; Congruent_OnlyPretone_Matrix(rowIdx_2, :)];
                    end                 
                end
           end

         end
     end
end

%% Average selected bands to output a single PSI value per channel pair (Congruent_OnlyPrior)


if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
Congruent_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Congruent_pfcIdx = 1:20
    for Congruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Congruent_OnlyPrior_averagedGrid(Congruent_pfcIdx, Congruent_acIdx) = columnAverage;
        end
    end
end
end



if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
Congruent_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Congruent_pfcIdx = 1:20
    for Congruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Congruent_OnlyPrior_PSI_grid{Congruent_pfcIdx, Congruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Congruent_OnlyPrior_averagedGrid(Congruent_pfcIdx, Congruent_acIdx) = columnAverage;
        end
    end
end
end


%% Average selected bands to output a single PSI value per channel pair (Congruent_OnlyPretone)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
Congruent_OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Congruent_pfcIdx = 1:20
    for Congruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Congruent_OnlyPretone_averagedGrid(Congruent_pfcIdx, Congruent_acIdx) = columnAverage;
        end
    end
end
end


if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
Congruent_OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Congruent_pfcIdx = 1:20
    for Congruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Congruent_OnlyPretone_PSI_grid{Congruent_pfcIdx, Congruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Congruent_OnlyPretone_averagedGrid(Congruent_pfcIdx, Congruent_acIdx) = columnAverage;
        end
    end
end
end

%% Incongruent_Analysis

% Step 3: Organize PSI values by channel pair. 

rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_24_PSI_Expectation_Analysis_testTone_Incongruent';      % This is the incongruent conditions                          
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

Incongruent_channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Incongruent_OnlyPrior_PSI_grid     = cell(20, 20); 
Incongruent_OnlyPretone_PSI_grid   = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair [INCONGRUENT]

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            Incongruent_PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', 'incongruent', '_', Frequency_Band,'_PSI.mat']));
            Incongruent_PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_','OnlyPretone_', Behavior,'_', 'incongruent', '_', Frequency_Band,'_PSI.mat']));

            if ~isempty(Incongruent_PSI_OnlyPrior_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(Incongruent_PSI_OnlyPrior_file)

        Incongruent_PSI_OnlyPrior_filename   = Incongruent_PSI_OnlyPrior_file(k).name;
        Incongruent_PSI_OnlyPretone_filename = Incongruent_PSI_OnlyPretone_file(k).name;

        temp_Incongruent_OnlyPrior_dir   = fullfile(session_path, Incongruent_PSI_OnlyPrior_filename);
        temp_Incongruent_OnlyPretone_dir = fullfile(session_path, Incongruent_PSI_OnlyPretone_filename);

        % Load the .mat file safely
         load(temp_Incongruent_OnlyPrior_dir,'PSI_OnlyPrior' );
         load(temp_Incongruent_OnlyPretone_dir,'PSI_OnlyPretone');

       Incongruent_OnlyPrior_Matrix     = PSI_OnlyPrior.psispctrm; 
       Incongruent_OnlyPrior_labels     = PSI_OnlyPrior.labelcmb;         
       Incongruent_OnlyPretone_Matrix   = PSI_OnlyPretone.psispctrm; 
       Incongruent_OnlyPretone_labels   = PSI_OnlyPretone.labelcmb;       

      end    

           % Populate Congruent_OnlyPrior_PSI_grid with data 

           for rowIdx_1 = 1:size(Incongruent_OnlyPrior_Matrix, 1)

            % Extract channel pair from label
            Incongruent_pfcLabel = Incongruent_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            Incongruent_acLabel  = Incongruent_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            Incongruent_pfcChannel = regexp(Incongruent_pfcLabel, 'PFC_ch\d+', 'match', 'once');    % Extract PFC channel
            Incongruent_acChannel  = regexp(Incongruent_acLabel,  'AC_ch\d+',  'match', 'once');    % Extract AC channel

            % Find indices in the grid
            Incongruent_pfcIdx = find(strcmp(Incongruent_channels,  strrep(Incongruent_pfcChannel, 'PFC_', '')));
            Incongruent_acIdx  = find(strcmp(Incongruent_channels,  strrep(Incongruent_acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(Incongruent_pfcIdx) && ~isempty(Incongruent_acIdx)
                    if isempty(Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx})
                    % Initialize the cell with the first row of data
                    Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx} = Incongruent_OnlyPrior_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx} = [Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx}; Incongruent_OnlyPrior_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate OnlyPretone_PSI_grid with data 

           for rowIdx_2 = 1:size(Incongruent_OnlyPretone_Matrix, 1)

            % Extract channel pair from label
            Incongruent_pfcLabel = Incongruent_OnlyPretone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            Incongruent_acLabel  = Incongruent_OnlyPretone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            Incongruent_pfcChannel = regexp(Incongruent_pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            Incongruent_acChannel  = regexp(Incongruent_acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            Incongruent_pfcIdx = find(strcmp(Incongruent_channels, strrep(Incongruent_pfcChannel, 'PFC_', '')));
            Incongruent_acIdx  = find(strcmp(Incongruent_channels, strrep(Incongruent_acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(Incongruent_pfcIdx) && ~isempty(Incongruent_acIdx)
                    if isempty(Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx})
                    % Initialize the cell with the first row of data
                    Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx} = Incongruent_OnlyPretone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx} = [Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx}; Incongruent_OnlyPretone_Matrix(rowIdx_2, :)];
                    end                 
                end
           end

         end
     end
end

%% Average selected bands to output a single PSI value per channel pair (Incongruent_OnlyPrior)


if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
Incongruent_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Incongruent_pfcIdx = 1:20
    for Incongruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Incongruent_OnlyPrior_averagedGrid(Incongruent_pfcIdx, Incongruent_acIdx) = columnAverage;
        end
    end
end
end



if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
Incongruent_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Incongruent_pfcIdx = 1:20
    for Incongruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Incongruent_OnlyPrior_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Incongruent_OnlyPrior_averagedGrid(Incongruent_pfcIdx, Incongruent_acIdx) = columnAverage;
        end
    end
end
end


%% Average selected bands to output a single PSI value per channel pair (Incongruent_OnlyPretone)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
Incongruent_OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Incongruent_pfcIdx = 1:20
    for Incongruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Incongruent_OnlyPretone_averagedGrid(Incongruent_pfcIdx, Incongruent_acIdx) = columnAverage;
        end
    end
end
end


if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;

% Initialize a new grid to store the averaged values
Incongruent_OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for Incongruent_pfcIdx = 1:20
    for Incongruent_acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx})
            % Extract the data matrix for the current position
            dataposition = Incongruent_OnlyPretone_PSI_grid{Incongruent_pfcIdx, Incongruent_acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            Incongruent_OnlyPretone_averagedGrid(Incongruent_pfcIdx, Incongruent_acIdx) = columnAverage;
        end
    end
end
end

%% Plot average PSI [Congruent]

% Replace NaNs with zero
Congruent_OnlyPrior_averagedGrid(isnan(Congruent_OnlyPrior_averagedGrid)) = 0;
Congruent_OnlyPretone_averagedGrid(isnan(Congruent_OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
globalMin = -0.02;
globalMax = 0.02;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Congruent_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Congruent_OnlyPretone_averagedGrid));

% averagedGrid_1_clamped = Congruent_OnlyPrior_averagedGrid;
% averagedGrid_2_clamped = Congruent_OnlyPretone_averagedGrid;

% Plot the first heatmap
figure(1);
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
%clim([-.08, .08]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPrior-congruent',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

% Plot the second heatmap
subplot(1, 2, 2);
imagesc(averagedGrid_2_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
%clim([-.08, .08]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPretone-congruent',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

%% Plot average PSI [INCONGRUENT]

% Replace NaNs with zero
Incongruent_OnlyPrior_averagedGrid(isnan(Incongruent_OnlyPrior_averagedGrid)) = 0;
Incongruent_OnlyPretone_averagedGrid(isnan(Incongruent_OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
globalMin = -0.02;
globalMax = 0.02;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Incongruent_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Incongruent_OnlyPretone_averagedGrid));

% averagedGrid_1_clamped = Congruent_OnlyPrior_averagedGrid;
% averagedGrid_2_clamped = Congruent_OnlyPretone_averagedGrid;

% Plot the first heatmap
figure(2);
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
%clim([-.08, .08]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPrior-incongruent-', Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

% Plot the second heatmap
subplot(1, 2, 2);
imagesc(averagedGrid_2_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
%clim([-.08, .08]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPretone-incongruent-',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


% %% Plot with Means + 95% CI, Legend, and High-Resolution Export [CONGRUENT]
% 
% Congruent_abs_OnlyPrior   = abs(Congruent_OnlyPrior_averagedGrid);
% Congruent_abs_OnlyPretone = abs(Congruent_OnlyPretone_averagedGrid);
% 
% % Mask out zeros (placeholder values)
% mask_nonzero_OnlyPrior   = Congruent_OnlyPrior_averagedGrid ~= 0;
% mask_nonzero_OnlyPretone = Congruent_OnlyPretone_averagedGrid ~= 0;
% 
% % Extract directional values
% Congruent_PFC_2_AC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Congruent_OnlyPrior_averagedGrid > 0);
% Congruent_AC_2_PFC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Congruent_OnlyPrior_averagedGrid < 0);
% Congruent_PFC_2_AC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Congruent_OnlyPretone_averagedGrid > 0);
% Congruent_AC_2_PFC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Congruent_OnlyPretone_averagedGrid < 0);
% 
% figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 11, 8.5]); hold on;
% 
% group_data = {Congruent_PFC_2_AC_OnlyPrior, Congruent_AC_2_PFC_OnlyPrior, Congruent_PFC_2_AC_OnlyPretone, Congruent_AC_2_PFC_OnlyPretone};
% means = cellfun(@mean, group_data);
% SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
% CIs   = SEMs * 1.96;  % 95% CI
% 
% bar_colors = [222 30 38; 71 133 197; 222 30 38; 71 133 197]/255;
% 
% % Bar plot
% for i = 1:4
%     b = bar(i, means(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
%     hold on;
% end
% 
% % Add CI error bars
% errorbar(1:4, means, CIs, 'k.', 'LineWidth', 1.5);
% 
% % Format axes
% set(gca, 'XTick', [1.5 3.5], ...
%          'XTickLabel', {'OnlyPrior','OnlyPretone'}, ...
%          'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold', ...
%          'YGrid', 'off');
% ylabel('Mean |PSI Value|', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% 
% % Format y-axis ticks (no scientific notation)
% ax = gca;
% ax.YAxis.Exponent = 0;
% ax.YAxis.TickLabelFormat = '%.4f';
% 
% % Add legend in top-left corner (no box)
% legend({'PFC→AC', 'AC→PFC'}, ...
%        'Location', 'northwest', ...
%        'FontSize', 14, ...
%        'FontName', 'Arial', ...
%        'Box', 'off');
% 
% 
% box off
% axis tight
% 
% % Optional: Export high-DPI image
% exportgraphics(gcf, 'PSI_BarPlot_Labeled_HighRes.png', 'Resolution', 600);
% 
% %% Plot with Means + 95% CI, Legend, and High-Resolution Export [Incongruent]
% 
% Congruent_abs_OnlyPrior   = abs(Incongruent_OnlyPrior_averagedGrid);
% Congruent_abs_OnlyPretone = abs(Incongruent_OnlyPretone_averagedGrid);
% 
% % Mask out zeros (placeholder values)
% mask_nonzero_OnlyPrior   = Incongruent_OnlyPrior_averagedGrid ~= 0;
% mask_nonzero_OnlyPretone = Incongruent_OnlyPretone_averagedGrid ~= 0;
% 
% % Extract directional values
% Congruent_PFC_2_AC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Incongruent_OnlyPrior_averagedGrid > 0);
% Congruent_AC_2_PFC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Incongruent_OnlyPrior_averagedGrid < 0);
% Congruent_PFC_2_AC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Incongruent_OnlyPretone_averagedGrid > 0);
% Congruent_AC_2_PFC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Incongruent_OnlyPretone_averagedGrid < 0);
% 
% figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 11, 8.5]); hold on;
% 
% group_data = {Congruent_PFC_2_AC_OnlyPrior, Congruent_AC_2_PFC_OnlyPrior, Congruent_PFC_2_AC_OnlyPretone, Congruent_AC_2_PFC_OnlyPretone};
% means = cellfun(@mean, group_data);
% SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
% CIs   = SEMs * 1.96;  % 95% CI
% 
% bar_colors = [222 30 38; 71 133 197; 222 30 38; 71 133 197]/255;
% 
% % Bar plot
% for i = 1:4
%     b = bar(i, means(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
%     hold on;
% end
% 
% % Add CI error bars
% errorbar(1:4, means, CIs, 'k.', 'LineWidth', 1.5);
% 
% % Format axes
% set(gca, 'XTick', [1.5 3.5], ...
%          'XTickLabel', {'OnlyPrior','OnlyPretone'}, ...
%          'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold', ...
%          'YGrid', 'off');
% ylabel('Mean |PSI Value|', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
% 
% % Format y-axis ticks (no scientific notation)
% ax = gca;
% ax.YAxis.Exponent = 0;
% ax.YAxis.TickLabelFormat = '%.4f';
% 
% % Add legend in top-left corner (no box)
% legend({'PFC→AC', 'AC→PFC'}, ...
%        'Location', 'northwest', ...
%        'FontSize', 14, ...
%        'FontName', 'Arial', ...
%        'Box', 'off');
% 
% 
% box off
% axis tight
% 
% % Optional: Export high-DPI image
% exportgraphics(gcf, 'PSI_BarPlot_Labeled_HighRes.png', 'Resolution', 600);

%% Grouped (touching) bars with aligned cluster labels + legend

% --- helper to pull |PSI| by direction, zero-masked ---
extract_dir = @(M) deal( ...
    abs(M(M~=0 & M > 0)), ... % PFC→AC
    abs(M(M~=0 & M < 0)));    % AC→PFC

% Congruent sets
[cOP_p2a, cOP_a2p] = extract_dir(Congruent_OnlyPrior_averagedGrid);
[cOT_p2a, cOT_a2p] = extract_dir(Congruent_OnlyPretone_averagedGrid);

% Incongruent sets
[iOP_p2a, iOP_a2p] = extract_dir(Incongruent_OnlyPrior_averagedGrid);
[iOT_p2a, iOT_a2p] = extract_dir(Incongruent_OnlyPretone_averagedGrid);

% Means + 95% CIs (omit NaNs)
mfun   = @(x) mean(x,'omitnan');
sfun   = @(x) std(x,'omitnan');
nfun   = @(x) sum(~isnan(x));
semfun = @(x) sfun(x) / max(1, sqrt(nfun(x)));
ci95   = @(x) 1.96 * semfun(x);

% Rows = groups on x; Cols = [Congruent, Incongruent]
means_mat = [ ...
    mfun(cOP_p2a), mfun(iOP_p2a);   % OnlyPrior · PFC→AC
    mfun(cOP_a2p), mfun(iOP_a2p);   % OnlyPrior · AC→PFC
    mfun(cOT_p2a), mfun(iOT_p2a);   % OnlyPretone · PFC→AC
    mfun(cOT_a2p), mfun(iOT_a2p)];  % OnlyPretone · AC→PFC

ci_mat = [ ...
    ci95(cOP_p2a), ci95(iOP_p2a);
    ci95(cOP_a2p), ci95(iOP_a2p);
    ci95(cOT_p2a), ci95(iOT_p2a);
    ci95(cOT_a2p), ci95(iOT_a2p)];

% --- plot ---
figure('Color','w','Units','inches','Position',[1,1,11,8.5]); hold on;
set(gca,'FontSize',16,'FontName','Arial','FontWeight','bold');

% Two touching bars per group
bh = bar(means_mat, 'grouped', 'BarWidth', 1); % 4×2 -> 4 clusters, 2 bars each

% Colors emphasize congruency
col_cong = [76 153 0]/255;     % Congruent
col_incg = [230 85 13]/255;    % Incongruent
set(bh(1), 'FaceColor', col_cong, 'EdgeColor','none');
set(bh(2), 'FaceColor', col_incg, 'EdgeColor','none');

% Error bars (R2019b+ uses XEndPoints/YEndPoints)
x1 = bh(1).XEndPoints;  y1 = bh(1).YEndPoints;  e1 = ci_mat(:,1);
x2 = bh(2).XEndPoints;  y2 = bh(2).YEndPoints;  e2 = ci_mat(:,2);
errorbar(x1, y1, e1, 'k.', 'LineWidth', 1.5);
errorbar(x2, y2, e2, 'k.', 'LineWidth', 1.5);

% ----- PERFECTLY ALIGNED X-AXIS LABELS -----
% Compute cluster centers from the two series endpoints
groupCenters = mean([x1(:) x2(:)], 2);

% Apply two-line labels centered on each cluster
ax = gca;
ax.XTick = groupCenters.';
ax.XTickLabel = { ...
    'OnlyPrior\newlinePFC→AC', ...
    'OnlyPrior\newlineAC→PFC', ...
    'OnlyPretone\newlinePFC→AC', ...
    'OnlyPretone\newlineAC→PFC' };
ax.TickLabelInterpreter = 'tex'; % enable \newline

% Y axis formatting
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%.4f';
ylabel('Mean |PSI Value|', 'FontSize',16, 'FontName','Arial', 'FontWeight','bold');

% Subtle separator between conditions (after cluster 2)
yl = ylim;
sepX = mean([groupCenters(2), groupCenters(3)]);
plot([sepX sepX], yl, '-', 'Color', [0 0 0 0.12], 'LineWidth', 2)

% Condition headers above halves
% text(mean(groupCenters(1:2)), yl(2) - 0.04*(yl(2)-yl(1)), 'OnlyPrior', ...
%      'HorizontalAlignment','center','FontName','Arial','FontSize',16,'FontWeight','bold');
% text(mean(groupCenters(3:4)), yl(2) - 0.04*(yl(2)-yl(1)), 'OnlyPretone', ...
%      'HorizontalAlignment','center','FontName','Arial','FontSize',16,'FontWeight','bold');

% ----- LEGEND (congruent vs incongruent) -----
legend(bh, {'Congruent','Incongruent'}, 'Location','northwest', ...
       'FontSize',14, 'FontName','Arial', 'Box','off');

box off; axis tight;

% Export (optional)
exportgraphics(gcf, 'PSI_BarPlot_Grouped_Aligned_CongVsIncong.png', 'Resolution', 600);



% %% Step 10: Nonparametric 2×2 ANOVA on Ranked PSI Values
% 
% PSI_values       = [Congruent_PFC_2_AC_OnlyPrior(:); Congruent_AC_2_PFC_OnlyPrior(:); Congruent_PFC_2_AC_OnlyPretone(:); Congruent_AC_2_PFC_OnlyPretone(:)];
% Expectation      = [repmat({'OnlyPrior'}, numel(Congruent_PFC_2_AC_OnlyPrior) + numel(Congruent_AC_2_PFC_OnlyPrior), 1);
%                     repmat({'OnlyPretone'}, numel(Congruent_PFC_2_AC_OnlyPretone) + numel(Congruent_AC_2_PFC_OnlyPretone), 1)];
% Direction        = [repmat({'PFC→AC'}, numel(Congruent_PFC_2_AC_OnlyPrior), 1);
%                     repmat({'AC→PFC'}, numel(Congruent_AC_2_PFC_OnlyPrior), 1);
%                     repmat({'PFC→AC'}, numel(Congruent_PFC_2_AC_OnlyPretone), 1);
%                     repmat({'AC→PFC'}, numel(Congruent_AC_2_PFC_OnlyPretone), 1)];
% 
% ranked_PSI = tiedrank(PSI_values);
% 
% [p, tbl, stats, terms]   = anovan(ranked_PSI, ...
%                                 {Expectation, Direction}, ...
%                                 'model', 2, ...
%                                 'varnames', {'Expectationion','Direction'}, ...
%                                 'display', 'on');
% 
% anova_table = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
% disp(anova_table);
% 
% %% Step 11: Follow-up Pairwise Rank-Sum Tests
% 
% % ----- Claim 1: OnlyPrior > OnlyPretone (all directions combined) -----
% OnlyPrior_all   = [Congruent_PFC_2_AC_OnlyPrior(:); Congruent_AC_2_PFC_OnlyPrior(:)];
% OnlyPretone_all = [Congruent_PFC_2_AC_OnlyPretone(:); Congruent_AC_2_PFC_OnlyPretone(:)];
% [p_expectation, ~, stats_expectation] = ranksum(OnlyPrior_all, OnlyPretone_all);
% fprintf('\nClaim 1 — OnlyPrior vs OnlyPretone:\n');
% fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
%     p_expectation, stats_expectation.zval, numel(OnlyPrior_all), numel(OnlyPretone_all));
% 
% % ----- Claim 2: Directional difference within OnlyPretone -----
% [p_OnlyPretone_dir, ~, stats_OnlyPretone_dir] = ranksum(Congruent_PFC_2_AC_OnlyPretone, Congruent_AC_2_PFC_OnlyPretone);
% fprintf('\nClaim 2 — OnlyPretone: PFC→AC vs AC→PFC:\n');
% fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
%     p_OnlyPretone_dir, stats_OnlyPretone_dir.zval, ...
%     numel(Congruent_PFC_2_AC_OnlyPretone), numel(Congruent_AC_2_PFC_OnlyPretone));
% 
% % ----- Claim 3: Directional difference within PreCue -----
% [p_OnlyPrior_dir, ~, stats_OnlyPrior_dir] = ranksum(Congruent_PFC_2_AC_OnlyPrior, Congruent_AC_2_PFC_OnlyPrior);
% fprintf('\nClaim 3 — OnlyPrior: PFC→AC vs AC→PFC:\n');
% fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
%     p_OnlyPrior_dir, stats_OnlyPrior_dir.zval, ...
%     numel(Congruent_PFC_2_AC_OnlyPrior), numel(Congruent_AC_2_PFC_OnlyPrior));


%% Congruent vs Incongruent (independent samples) per (Condition × Direction)
% Inputs needed in workspace:
%   Congruent_OnlyPrior_averagedGrid
%   Congruent_OnlyPretone_averagedGrid
%   Incongruent_OnlyPrior_averagedGrid
%   Incongruent_OnlyPretone_averagedGrid
%
% Convention: positive PSI -> PFC→AC; negative PSI -> AC→PFC.
% Zeros are placeholders and are excluded.

% --- helpers ---
extract_dirs = @(M) deal(abs(M(M>0)), abs(M(M<0)));  % returns [PFC→AC, AC→PFC]
mclean = @(v) v(~isnan(v) & v~=0);

% Split by cell
[cOP_p2a, cOP_a2p] = extract_dirs(Congruent_OnlyPrior_averagedGrid);
[iOP_p2a, iOP_a2p] = extract_dirs(Incongruent_OnlyPrior_averagedGrid);

[cOT_p2a, cOT_a2p] = extract_dirs(Congruent_OnlyPretone_averagedGrid);
[iOT_p2a, iOT_a2p] = extract_dirs(Incongruent_OnlyPretone_averagedGrid);

labels = { ...
    'OnlyPrior · PFC→AC', ...
    'OnlyPrior · AC→PFC', ...
    'OnlyPretone · PFC→AC', ...
    'OnlyPretone · AC→PFC'};

X = { cOP_p2a, cOP_a2p, cOT_p2a, cOT_a2p };  % Congruent
Y = { iOP_p2a, iOP_a2p, iOT_p2a, iOT_a2p };  % Incongruent

% Run rank-sum tests
m = numel(X);
p_raw = nan(m,1);
W     = nan(m,1);
Z     = nan(m,1);
nC    = nan(m,1);
nI    = nan(m,1);

for k = 1:m
    x = mclean(X{k});
    y = mclean(Y{k});
    nC(k) = numel(x); 
    nI(k) = numel(y);

    if isempty(x) || isempty(y)
        % leave NaNs; nothing to test
        continue
    end

    [p,~,stats] = ranksum(x,y);     % independent-samples nonparametric "t-test"
    p_raw(k) = p;
    W(k)     = stats.ranksum;       % sum of ranks for x (congruent)
    Z(k)     = stats.zval;          % normal-approx z
end

% Holm–Bonferroni (FWER) across the 4 tests
[ps, idx] = sort(p_raw);
adj = (m - (1:m)' + 1) .* ps;       % step-down
adj = min(1, adj);
for i = 2:m
    adj(i) = max(adj(i), adj(i-1)); % enforce monotonicity
end
p_holm = nan(m,1);
p_holm(idx) = adj;

% Tidy table
T = table(labels(:), nC, nI, W, Z, p_raw, p_holm, ...
    'VariableNames', {'Cell','N_Cong','N_Incong','RankSum_W','Z','p_raw','p_Holm'});

disp(T)

% Optional: save
% writetable(T, 'Congruent_vs_Incongruent_Holm.csv');

%% Paired nonparametric tests (Wilcoxon signed-rank):
% Congruent vs Incongruent within each (Condition × Direction)
% Inputs needed in workspace:
%   Congruent_OnlyPrior_averagedGrid
%   Congruent_OnlyPretone_averagedGrid
%   Incongruent_OnlyPrior_averagedGrid
%   Incongruent_OnlyPretone_averagedGrid
%
% Pairing rule (per direction):
%   PFC→AC: use positions where BOTH matrices > 0
%   AC→PFC: use positions where BOTH matrices < 0
% Zeros are placeholders and excluded.

% --- helper: build paired |PSI| vectors for a direction ---
% dirFlag = +1 for PFC→AC (both >0), -1 for AC→PFC (both <0)
pair_dir = @(MC, MI, dirFlag) deal( ...
    abs(MC( (dirFlag>0) & MC>0 & MI>0 | (dirFlag<0) & MC<0 & MI<0 )), ...
    abs(MI( (dirFlag>0) & MC>0 & MI>0 | (dirFlag<0) & MC<0 & MI<0 )) );

% Clean function (drop NaN/zero just in case)
mclean = @(v) v(~isnan(v) & v~=0);

% ---- Build pairs for each (Condition × Direction) ----
% OnlyPrior
[x1, y1] = pair_dir(Congruent_OnlyPrior_averagedGrid,   Incongruent_OnlyPrior_averagedGrid,   +1); % PFC→AC
[x2, y2] = pair_dir(Congruent_OnlyPrior_averagedGrid,   Incongruent_OnlyPrior_averagedGrid,   -1); % AC→PFC

% OnlyPretone
[x3, y3] = pair_dir(Congruent_OnlyPretone_averagedGrid, Incongruent_OnlyPretone_averagedGrid, +1); % PFC→AC
[x4, y4] = pair_dir(Congruent_OnlyPretone_averagedGrid, Incongruent_OnlyPretone_averagedGrid, -1); % AC→PFC

labels = { ...
    'OnlyPrior · PFC→AC', ...
    'OnlyPrior · AC→PFC', ...
    'OnlyPretone · PFC→AC', ...
    'OnlyPretone · AC→PFC'};

X = {x1, x2, x3, x4};   % Congruent |PSI|
Y = {y1, y2, y3, y4};   % Incongruent |PSI|

% ---- Run Wilcoxon signed-rank tests (paired) ----
m = numel(X);
Npairs = nan(m,1);
W      = nan(m,1);   % signed-rank statistic (sum of positive signed ranks)
Z      = nan(m,1);   % normal approx z
p_raw  = nan(m,1);

for k = 1:m
    x = mclean(X{k});
    y = mclean(Y{k});

    % keep only indices where both present (pairing already enforced by masks)
    n = min(numel(x), numel(y));
    if n==0
        continue
    end
    % ensure equal length (should already be equal if built from same mask sizes,
    % but guard against any mismatch due to earlier cleaning)
    x = x(1:n);
    y = y(1:n);

    Npairs(k) = n;

    % Signed-rank on paired differences of |PSI| (Congruent - Incongruent)
    d = x - y;
    if all(d==0 | isnan(d))
        % nothing to test
        continue
    end

    [p,~,stats] = signrank(x, y);   % two-sided by default
    p_raw(k) = p;
    W(k)     = stats.signedrank;    % sum of ranks of positive differences
    if isfield(stats,'zval')
        Z(k) = stats.zval;          % available in newer MATLAB
    end
end

% ---- Holm–Bonferroni correction across the 4 tests ----
[ps, idx] = sort(p_raw);
adj = (m - (1:m)' + 1) .* ps;
adj = min(1, adj);
for i = 2:m
    adj(i) = max(adj(i), adj(i-1));
end
p_holm = nan(m,1);
p_holm(idx) = adj;

% ---- Results table ----
T = table(labels(:), Npairs, W, Z, p_raw, p_holm, ...
    'VariableNames', {'Cell','N_Pairs','SignedRank_W','Z','p_raw','p_Holm'});

disp(T)

% Optional: save
% writetable(T, 'Congruent_vs_Incongruent_SignedRank_Holm.csv');

%%

Cat = Congruent_OnlyPretone_PSI_grid(:);
