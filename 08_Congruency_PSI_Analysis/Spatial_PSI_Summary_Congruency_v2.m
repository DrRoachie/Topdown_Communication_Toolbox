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

%% Step 1: Book-keeping for CONGRUENT analysis

Frequency_Band    = 'theta';
Behavior          = 'correct';
animals           = {'MrCassius'; 'MrM'};                                                   % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_07_10_PSI_Expectation_Analysis_testTone_Congruent';      % This is the congruent                            
sessions          = dir(fullfile(rootdir, '19*'));


%% Step 2: Loop through all sessions and sort PSI values by channel-pair [CONGRUENT]

% Define the PFC and AC channel labels and generate temporary grids

Congruent_channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Congruent_OnlyPrior_PSI_grid     = cell(20, 20); 
Congruent_OnlyPretone_PSI_grid   = cell(20, 20);

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

%% Step 3: Average selected bands to output a single PSI value per channel pair (Congruent_OnlyPrior)

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


% Average selected bands to output a single PSI value per channel pair (Congruent_OnlyPretone)

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



%% Step 5: Book-keeping for Incongruent_Analysis

% Step 3: Organize PSI values by channel pair. 

rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_24_PSI_Expectation_Analysis_testTone_Incongruent';      % This is the incongruent conditions                          
sessions          = dir(fullfile(rootdir, '19*'));


%% Step 6: Loop through all sessions and sort PSI values by channel-pair [INCONGRUENT]

% Define the PFC and AC channel labels and generate temporary grids

Incongruent_channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Incongruent_OnlyPrior_PSI_grid     = cell(20, 20); 
Incongruent_OnlyPretone_PSI_grid   = cell(20, 20);

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

%% Step 7: Average selected bands to output a single PSI value per channel pair (Incongruent_OnlyPrior)

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

% Average selected bands to output a single PSI value per channel pair (Incongruent_OnlyPretone)

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



%% Step 8: Plot average PSI [Congruent]

% Replace NaNs with zero
Congruent_OnlyPrior_averagedGrid(isnan(Congruent_OnlyPrior_averagedGrid)) = 0;
Congruent_OnlyPretone_averagedGrid(isnan(Congruent_OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
globalMin = -0.02;
globalMax = 0.02;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Congruent_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Congruent_OnlyPretone_averagedGrid));

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
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPretone-congruent',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

%% Step 9 Plot average PSI [INCONGRUENT]

% Replace NaNs with zero
Incongruent_OnlyPrior_averagedGrid(isnan(Incongruent_OnlyPrior_averagedGrid)) = 0;
Incongruent_OnlyPretone_averagedGrid(isnan(Incongruent_OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
globalMin = -0.02;
globalMax = 0.02;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Incongruent_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Incongruent_OnlyPretone_averagedGrid));


% Plot the first heatmap
figure(2);
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
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
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['OnlyPretone-incongruent-',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%% Step 10: Extract Directional Data for Summary Plot with Means + 95% CI, Legend, and High-Resolution Export [CONGRUENT]

Congruent_abs_OnlyPrior   = abs(Congruent_OnlyPrior_averagedGrid);
Congruent_abs_OnlyPretone = abs(Congruent_OnlyPretone_averagedGrid);

% Mask out zeros (placeholder values)
mask_nonzero_OnlyPrior   = Congruent_OnlyPrior_averagedGrid ~= 0;
mask_nonzero_OnlyPretone = Congruent_OnlyPretone_averagedGrid ~= 0;

% Extract directional values
Congruent_PFC_2_AC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Congruent_OnlyPrior_averagedGrid > 0);
Congruent_AC_2_PFC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Congruent_OnlyPrior_averagedGrid < 0);
Congruent_PFC_2_AC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Congruent_OnlyPretone_averagedGrid > 0);
Congruent_AC_2_PFC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Congruent_OnlyPretone_averagedGrid < 0);

group_data = {Congruent_PFC_2_AC_OnlyPrior, Congruent_AC_2_PFC_OnlyPrior, Congruent_PFC_2_AC_OnlyPretone, Congruent_AC_2_PFC_OnlyPretone};
means = cellfun(@mean, group_data);
SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
CIs   = SEMs * 1.96;  % 95% CI



%% Step 11: Extract Directional Data for Summary Plot with Means + 95% CI, Legend, and High-Resolution Export [Incongruent]

Congruent_abs_OnlyPrior   = abs(Incongruent_OnlyPrior_averagedGrid);
Congruent_abs_OnlyPretone = abs(Incongruent_OnlyPretone_averagedGrid);

% Mask out zeros (placeholder values)
mask_nonzero_OnlyPrior   = Incongruent_OnlyPrior_averagedGrid ~= 0;
mask_nonzero_OnlyPretone = Incongruent_OnlyPretone_averagedGrid ~= 0;

% Extract directional values
Congruent_PFC_2_AC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Incongruent_OnlyPrior_averagedGrid > 0);
Congruent_AC_2_PFC_OnlyPrior      = Congruent_abs_OnlyPrior(mask_nonzero_OnlyPrior & Incongruent_OnlyPrior_averagedGrid < 0);
Congruent_PFC_2_AC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Incongruent_OnlyPretone_averagedGrid > 0);
Congruent_AC_2_PFC_OnlyPretone    = Congruent_abs_OnlyPretone(mask_nonzero_OnlyPretone & Incongruent_OnlyPretone_averagedGrid < 0);

group_data = {Congruent_PFC_2_AC_OnlyPrior, Congruent_AC_2_PFC_OnlyPrior, Congruent_PFC_2_AC_OnlyPretone, Congruent_AC_2_PFC_OnlyPretone};
means = cellfun(@mean, group_data);
SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
CIs   = SEMs * 1.96;  % 95% CI


%% Step 12: Grouped (touching) bars with aligned cluster labels + legend

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

% % Export (optional)
% exportgraphics(gcf, 'PSI_BarPlot_Grouped_Aligned_CongVsIncong.png', 'Resolution', 600);


%% Step 13: Independent nonparametric 2×2 (SRH): Congruency × Condition + 4 pairwise tests
% Factors:
%   - Congruency: Congruent vs Incongruent
%   - Condition:  OnlyPrior vs OnlyPretone
% DV: |PSI| magnitudes, observations treated as independent

% ---------- Assemble data for SRH (collapse directions into magnitude) ----------
OP_Cong = [cOP_p2a(:); cOP_a2p(:)];   % OnlyPrior · Congruent (both directions)
OP_Incg = [iOP_p2a(:); iOP_a2p(:)];   % OnlyPrior · Incongruent
OT_Cong = [cOT_p2a(:); cOT_a2p(:)];   % OnlyPretone · Congruent
OT_Incg = [iOT_p2a(:); iOT_a2p(:)];   % OnlyPretone · Incongruent

% Remove NaNs (if any)
OP_Cong = OP_Cong(~isnan(OP_Cong));
OP_Incg = OP_Incg(~isnan(OP_Incg));
OT_Cong = OT_Cong(~isnan(OT_Cong));
OT_Incg = OT_Incg(~isnan(OT_Incg));

Y = [OP_Cong; OP_Incg; OT_Cong; OT_Incg];

Congruency = [ ...
    repmat({'Congruent'},  numel(OP_Cong), 1); ...
    repmat({'Incongruent'},numel(OP_Incg), 1); ...
    repmat({'Congruent'},  numel(OT_Cong), 1); ...
    repmat({'Incongruent'},numel(OT_Incg), 1)];

Condition = [ ...
    repmat({'OnlyPrior'},   numel(OP_Cong), 1); ...
    repmat({'OnlyPrior'},   numel(OP_Incg), 1); ...
    repmat({'OnlyPretone'}, numel(OT_Cong), 1); ...
    repmat({'OnlyPretone'}, numel(OT_Incg), 1)];

% ---------- Scheirer–Ray–Hare (independent 2×2 on ranks) ----------
SRH = srh_two_way(Y, Congruency, Condition); % returns a struct with table-like fields

fprintf('\n=== Scheirer–Ray–Hare (2×2, independent) on |PSI| ===\n');
fprintf('Factor                H        df        p\n');
fprintf('---------------------------------------------\n');
fprintf('%-20s %8.4f   %2d   %10.6f\n','Congruency', SRH.HA, SRH.dfA, SRH.pA);
fprintf('%-20s %8.4f   %2d   %10.6f\n','Condition',  SRH.HB, SRH.dfB, SRH.pB);
fprintf('%-20s %8.4f   %2d   %10.6f\n','A×B',        SRH.HAB,SRH.dfAB,SRH.pAB);
fprintf('N = %d (tie-correction C = %.6f)\n\n', SRH.N, SRH.tie_correction);

% ---------- 4 pairwise Wilcoxon rank-sum (Congruent vs Incongruent) with Holm ----------
pairs = { ...
    'OnlyPrior · PFC→AC',   cOP_p2a, iOP_p2a; ...
    'OnlyPrior · AC→PFC',   cOP_a2p, iOP_a2p; ...
    'OnlyPretone · PFC→AC', cOT_p2a, iOT_p2a; ...
    'OnlyPretone · AC→PFC', cOT_a2p, iOT_a2p };

nPairs = size(pairs,1);
rows(nPairs,1) = struct('Cell',[],'N_Cong',[],'N_Incong',[],'RankSum_W',[],'Z',[],'p_raw',[],'p_Holm',[]);
pvals = nan(nPairs,1);

for r = 1:nPairs
    label = pairs{r,1};
    xCong = pairs{r,2}; xCong = xCong(~isnan(xCong));
    xIncg = pairs{r,3}; xIncg = xIncg(~isnan(xIncg));

    % Wilcoxon rank-sum (Mann–Whitney), independent groups
    [p,~,stats] = ranksum(xCong, xIncg); % stats.ranksum (W) and stats.zval (Z)
    pvals(r) = p;

    rows(r).Cell      = label;
    rows(r).N_Cong    = numel(xCong);
    rows(r).N_Incong  = numel(xIncg);
    rows(r).RankSum_W = stats.ranksum;
    rows(r).Z         = stats.zval;
    rows(r).p_raw     = p;
end

% Holm–Bonferroni (step-down)
p_holm = holm_bonferroni(pvals);
for r = 1:nPairs
    rows(r).p_Holm = p_holm(r);
end

% Display table
PairwiseTbl = struct2table(rows);
disp('=== Pairwise: Congruent vs Incongruent within each Condition–Direction ===');
disp(PairwiseTbl(:,{'Cell','N_Cong','N_Incong','RankSum_W','Z','p_raw','p_Holm'}));

% ---------- Helper functions (SRH + Holm) ----------
function S = srh_two_way(Y, A, B)
    % A, B are cellstr/categorical with 2 levels each (generalizable)
    % Implements Scheirer–Ray–Hare with tie correction.
    Y = Y(:);
    A = categorical(A(:));
    B = categorical(B(:));

    % Ranks with average ties
    [R, tie_correction] = rank_with_ties(Y);
    N = numel(Y);

    levelsA = categories(A); a = numel(levelsA);
    levelsB = categories(B); b = numel(levelsB);

    % Rank sums by A, by B, and by A×B cells
    Ra = zeros(a,1); na = zeros(a,1);
    for i = 1:a
        idx = A==levelsA{i};
        Ra(i) = nansum(R(idx));
        na(i) = nnz(idx);
    end

    Rb = zeros(b,1); nb = zeros(b,1);
    for j = 1:b
        idx = B==levelsB{j};
        Rb(j) = nansum(R(idx));
        nb(j) = nnz(idx);
    end

    Rab = zeros(a,b); nab = zeros(a,b);
    for i = 1:a
        for j = 1:b
            idx = (A==levelsA{i}) & (B==levelsB{j});
            Rab(i,j) = nansum(R(idx));
            nab(i,j) = nnz(idx);
        end
    end

    C = tie_correction;                   % tie correction factor
    K = 12 / (N*(N+1));
    grand = 3*(N+1);

    HA  = (K * nansum( (Ra.^2) ./ max(1,na) ))/C  - grand;
    HB  = (K * nansum( (Rb.^2) ./ max(1,nb) ))/C  - grand;
    HAB = (K * nansum( (Rab.^2) ./ max(1,nab), 'all'))/C - HA - HB - grand;

    dfA  = a-1;
    dfB  = b-1;
    dfAB = (a-1)*(b-1);

    pA  = 1 - chi2cdf(HA,  dfA);
    pB  = 1 - chi2cdf(HB,  dfB);
    pAB = 1 - chi2cdf(HAB, dfAB);

    S = struct('HA',HA,'HB',HB,'HAB',HAB, ...
               'dfA',dfA,'dfB',dfB,'dfAB',dfAB, ...
               'pA',pA,'pB',pB,'pAB',pAB, ...
               'N',N,'tie_correction',C, ...
               'levelsA',{levelsA},'levelsB',{levelsB});
end

function [R, C] = rank_with_ties(x)
    % Return average ranks and tie-correction C (Kruskal–Wallis style)
    [sx, idx] = sort(x);
    ranks = (1:numel(x))';
    R = nan(size(x));
    % average tied ranks
    k = 1;
    while k <= numel(x)
        t = find(sx(k:end) ~= sx(k), 1, 'first');
        if isempty(t), t = numel(x)-k+1; end
        span = k:(k+t-1);
        meanRank = mean(ranks(span));
        R(idx(span)) = meanRank;
        k = k + t;
    end
    % tie correction
    % sum over tie groups of (t^3 - t)
    tie_counts = [];
    k = 1;
    while k <= numel(x)
        t = find(sx(k:end) ~= sx(k), 1, 'first');
        if isempty(t), t = numel(x)-k+1; end
        if t > 1, tie_counts(end+1) = t; %#ok<AGROW>
        end
        k = k + t;
    end
    N = numel(x);
    if isempty(tie_counts)
        C = 1;
    else
        C = 1 - sum(tie_counts.^3 - tie_counts) / (N^3 - N);
    end
end

function p_adj = holm_bonferroni(p)
    % Step-down Holm adjustment
    m = numel(p);
    [ps, ord] = sort(p(:));
    p_adj = nan(m,1);
    running_max = 0;
    for k = 1:m
        adj = (m - k + 1) * ps(k);
        adj = min(1, adj);
        running_max = max(running_max, adj);
        p_adj(k) = running_max;
    end
    % map back
    p_adj(ord) = p_adj;
end


