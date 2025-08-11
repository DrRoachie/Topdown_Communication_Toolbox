%% Step 3: Organize PSI values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
Behavior          = 'correct';
Congruency        = 'congruent';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius'};                                                     % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\Expectation_Analysis';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
OnlyPrior_PSI_grid   = cell(20, 20); 
OnlyPretone_PSI_grid = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',   Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));
            PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));

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
         load(temp_OnlyPrior_dir,  'PSI_OnlyPrior' );
         load(temp_OnlyPretone_dir,'PSI_OnlyPretone');

       OnlyPrior_Matrix     = PSI_OnlyPrior.psispctrm; 
       OnlyPrior_labels     = PSI_OnlyPrior.labelcmb;         
       OnlyPretone_Matrix   = PSI_OnlyPretone.psispctrm; 
       OnlyPretone_labels   = PSI_OnlyPretone.labelcmb;       

       % === Request 1: Load SpecBalanced_theta ===
        SpecBalanced_rootdir = 'C:\Users\Corey Roach\Documents\00_DATA\SpecBalanced_Eval_Overlap_Condition';
        SpecBalanced_filename = fullfile(SpecBalanced_rootdir, RecDate, Animal, ...
            sprintf('OverlappingLabelsConditions_%s_%s.mat', RecDate, Animal));
        if exist(SpecBalanced_filename, 'file')
            load(SpecBalanced_filename, 'OverlappingLabelsConditions');
        else
            warning('Missing OverlappingLabelsConditions file: %s — skipping session.', SpecBalanced_filename);
            continue;
        end
        
        % === Request 2 (finalized): Generate all possible PFC-AC label pairs as nx2 cell array ===
        
        % OnlyPrior combinations
        [PFC_prior, AC_prior] = ndgrid(OverlappingLabelsConditions.OnlyPrior.PFC(:), OverlappingLabelsConditions.OnlyPrior.AC(:));
        OnlyPrior_pairs       = [PFC_prior(:), AC_prior(:)];
        
        % OnlyPretone combinations
        [PFC_pretone, AC_pretone] = ndgrid(OverlappingLabelsConditions.OnlyPretone.PFC(:), OverlappingLabelsConditions.OnlyPretone.AC(:));
        OnlyPretone_pairs         = [PFC_pretone(:), AC_pretone(:)];
        

        % % === Request 3: Clean PSI labels and match with SpecBalanced channel pairs ===
        % 
        % Step 1: Clean PSI labels by replacing 'D#_' with '*'
        % clean_OnlyPrior_labels = cellfun(@(x) regexprep(x, '^D\d_', '*'), ...
        %                                  PSI_OnlyPrior.labelcmb, 'UniformOutput', false);
        % clean_OnlyPretone_labels = cellfun(@(x) regexprep(x, '^D\d_', '*'), ...
        %                                    PSI_OnlyPretone.labelcmb, 'UniformOutput', false);
        % 
        % Step 2: Reshape into n x 2 cell arrays
        % clean_OnlyPrior_labels   = reshape(clean_OnlyPrior_labels, [], 2);
        % clean_OnlyPretone_labels = reshape(clean_OnlyPretone_labels, [], 2);
        % 
        % 
        % % Step 3: Create logical indices for valid rows using row-wise check
        % Shared_OnlyPrior_idx = false(size(clean_OnlyPrior_labels,1),1);
        % for i = 1:size(clean_OnlyPrior_labels,1)
        %     Shared_OnlyPrior_idx(i) = any( ...
        %         cellfun(@(a,b) strcmp(clean_OnlyPrior_labels{i,1}, a) && strcmp(clean_OnlyPrior_labels{i,2}, b), ...
        %                 OnlyPrior_pairs(:,1), OnlyPrior_pairs(:,2)) );
        % end
        % 
        % Shared_OnlyPretone_idx = false(size(clean_OnlyPretone_labels,1),1);
        % for i = 1:size(clean_OnlyPretone_labels,1)
        %     Shared_OnlyPretone_idx(i) = any( ...
        %         cellfun(@(a,b) strcmp(clean_OnlyPretone_labels{i,1}, a) && strcmp(clean_OnlyPretone_labels{i,2}, b), ...
        %                 OnlyPretone_pairs(:,1), OnlyPretone_pairs(:,2)) );
        % end
        % 
        % % Step 4: Subset PSI matrices and original labelcmb
        % OnlyPrior_Matrix = OnlyPrior_Matrix(Shared_OnlyPrior_idx, :);
        % OnlyPrior_labels = PSI_OnlyPrior.labelcmb(Shared_OnlyPrior_idx, :);
        % 
        % OnlyPretone_Matrix = OnlyPretone_Matrix(Shared_OnlyPretone_idx, :);
        % OnlyPretone_labels = PSI_OnlyPretone.labelcmb(Shared_OnlyPretone_idx, :);

        % === Request 3: Directly find shared channel pairs ===
        
        % Convert both to keys like '*PFC_ch03-*AC_ch06' for row-wise comparison
        OnlyPrior_label_keys   = strcat(OnlyPrior_labels(:,1), '-', OnlyPrior_labels(:,2));
        OnlyPretone_label_keys = strcat(OnlyPretone_labels(:,1), '-', OnlyPretone_labels(:,2));
        
        OnlyPrior_pair_keys    = strcat(OnlyPrior_pairs(:,1), '-', OnlyPrior_pairs(:,2));
        OnlyPretone_pair_keys  = strcat(OnlyPretone_pairs(:,1), '-', OnlyPretone_pairs(:,2));
        
        % Find intersection
        Shared_OnlyPrior_idx   = ismember(OnlyPrior_label_keys,   OnlyPrior_pair_keys);
        Shared_OnlyPretone_idx = ismember(OnlyPretone_label_keys, OnlyPretone_pair_keys);
        
        % Subset PSI data and labelcmbs
        OnlyPrior_Matrix = OnlyPrior_Matrix(Shared_OnlyPrior_idx, :);
        OnlyPrior_labels = OnlyPrior_labels(Shared_OnlyPrior_idx, :);
        
        OnlyPretone_Matrix = OnlyPretone_Matrix(Shared_OnlyPretone_idx, :);
        OnlyPretone_labels = OnlyPretone_labels(Shared_OnlyPretone_idx, :);


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
OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
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
OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
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
OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
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
OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
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
OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
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
OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
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
globalMin = -0.03;
globalMax = 0.03;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, OnlyPretone_averagedGrid));

% Plot the first heatmap
figure(1);
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
%colormap(jet)
clim([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['SpecBalanced-OnlyPrior-',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
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
title(['SpecBalanced-OnlyPretone-',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;
