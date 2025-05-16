%% Step 3: Organize Coherence values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
Behavior          = 'correct';
Congruency        = 'congruent';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrM'};                                                     % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\Expectation_Analysis\01_Belt';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Mean_OnlyPrior_Coh_grid   = cell(20, 20); 
Mean_OnlyPretone_Coh_grid = cell(20, 20);
Z_OnlyPrior_Coh_grid   = cell(20, 20); 
Z_OnlyPretone_Coh_grid = cell(20, 20);
ZDiff_Coh_grid     = cell(20, 20); 


%% Loop through all sessions and sort Coherence values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            Coherence_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', Congruency, '_', Frequency_Band,'_Coh.mat']));
            Coherence_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_','OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_Coh.mat']));

            if ~isempty(Coherence_OnlyPrior_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(Coherence_OnlyPrior_file)

        Coh_OnlyPrior_filename   = Coherence_OnlyPrior_file(k).name;
        Coh_OnlyPretone_filename = Coherence_OnlyPretone_file(k).name;

        temp_OnlyPrior_dir   = fullfile(session_path, Coh_OnlyPrior_filename);
        temp_OnlyPretone_dir = fullfile(session_path, Coh_OnlyPretone_filename);

        % Load the .mat file safely
         load(temp_OnlyPrior_dir,'Coh_OnlyPrior' );
         load(temp_OnlyPretone_dir,'Coh_OnlyPretone');

       Mean_OnlyPrior_Matrix   = Coh_OnlyPrior.cohspctrm; 
       Mean_OnlyPretone_Matrix = Coh_OnlyPretone.cohspctrm; 
       Z_OnlyPrior_Matrix      = Coh_OnlyPrior.cohspctrm; 
       Z_OnlyPrior_labels      = Coh_OnlyPrior.labelcmb;         
       Z_OnlyPretone_Matrix    = Coh_OnlyPretone.cohspctrm; 
       Z_OnlyPretone_labels    = Coh_OnlyPretone.labelcmb;  

        % z-score the data 

        % Step 1: Pool all values from both matrices
      
        pooledData = [Z_OnlyPrior_Matrix(:); Z_OnlyPretone_Matrix(:)];

        % Step 2: Compute global mean and std
        globalMean = mean(pooledData);
        globalStd  = std(pooledData);


        % Step 3: Z-score both matrices
        Z_OnlyPrior_Matrix   = (Z_OnlyPrior_Matrix   - globalMean) / globalStd;
        Z_OnlyPretone_Matrix = (Z_OnlyPretone_Matrix - globalMean) / globalStd;
        zDiff_Matrix         = Z_OnlyPrior_Matrix - Z_OnlyPretone_Matrix;

      end    

           % Populate Z_OnlyPrior_Coherence_grid with data 

           for rowIdx_1 = 1:size(Z_OnlyPrior_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Z_OnlyPrior_Coh_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = Mean_OnlyPrior_Matrix(rowIdx_1, :);
                    Z_OnlyPrior_Coh_grid{pfcIdx, acIdx} = Z_OnlyPrior_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx}; Mean_OnlyPrior_Matrix(rowIdx_1, :)];
                    Z_OnlyPrior_Coh_grid{pfcIdx, acIdx} = [Z_OnlyPrior_Coh_grid{pfcIdx, acIdx}; Z_OnlyPrior_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate OnlyPretone_PSI_grid with data 

           for rowIdx_2 = 1:size(Z_OnlyPretone_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_OnlyPretone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_OnlyPretone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Z_OnlyPretone_Coh_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx} = Mean_OnlyPretone_Matrix(rowIdx_2, :);
                    Z_OnlyPretone_Coh_grid{pfcIdx, acIdx} = Z_OnlyPretone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx}; Mean_OnlyPretone_Matrix(rowIdx_2, :)];
                    Z_OnlyPretone_Coh_grid{pfcIdx, acIdx} = [Z_OnlyPretone_Coh_grid{pfcIdx, acIdx}; Z_OnlyPretone_Matrix(rowIdx_2, :)];
                    end                 
                end
           end
           
           % Populate ZDiff_Coherence_PSI_grid with data 

           for rowIdx_1 = 1:size(zDiff_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(ZDiff_Coh_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    ZDiff_Coh_grid{pfcIdx, acIdx} = zDiff_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    ZDiff_Coh_grid{pfcIdx, acIdx} = [ZDiff_Coh_grid{pfcIdx, acIdx}; zDiff_Matrix(rowIdx_1, :)];
                    end

                end
           end
            

         end
     end
end

%% Average selected bands to output a single Coherence value per channel pair (OnlyPrior)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;
end

% Initialize a new grid to store the averaged values
Z_OnlyPrior_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions
Mean_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions
std_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_OnlyPrior_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition    = Z_OnlyPrior_Coh_grid{pfcIdx, acIdx};
            Mean_dataposition = Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx};
            
            % Compute the average of the selected columns
            Z_columnAverage    = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage = mean(Mean_dataposition(:, columnRange), 'all', 'omitnan');
            std_columnAverage = std(Mean_dataposition(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_OnlyPrior_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage;
            Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage;
            std_OnlyPrior_averagedGrid(pfcIdx, acIdx)  = std_columnAverage;

        end
    end
end

%% Average selected bands to output a single Coherence value per channel pair (OnlyPretone)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;
end

% Initialize a new grid to store the averaged values
Z_OnlyPretone_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions
Mean_OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions
std_OnlyPretone_averagedGrid  = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_OnlyPretone_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition    = Z_OnlyPretone_Coh_grid{pfcIdx, acIdx};
            Mean_dataposition = Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage    = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage = mean(Mean_dataposition(:, columnRange), 'all', 'omitnan');
            std_columnAverage = std(Mean_dataposition(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_OnlyPretone_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage;
            Mean_OnlyPretone_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage;
            std_OnlyPretone_averagedGrid(pfcIdx, acIdx)  = std_columnAverage;
        end
    end
end

%% Average selected bands to output a single Coherence value per channel pair (ZDiff)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;
end

% Initialize a new grid to store the averaged values
ZDiff_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(ZDiff_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition = ZDiff_Coh_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            ZDiff_averagedGrid(pfcIdx, acIdx) = Z_columnAverage;
        end
    end
end


%% Plot the Z-scored Coherhence  

% Replace NaNs with zero
Z_OnlyPrior_averagedGrid(isnan(Z_OnlyPrior_averagedGrid)) = 0;
Z_OnlyPretone_averagedGrid(isnan(Z_OnlyPretone_averagedGrid)) = 0;
ZDiff_averagedGrid(isnan(Z_OnlyPrior_averagedGrid)) = 0;

% Determine the global color range
globalMin = -7;
globalMax = 7;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Z_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Z_OnlyPretone_averagedGrid));
averagedGrid_3_clamped = max(globalMin, min(globalMax, ZDiff_averagedGrid));

% Plot the first heatmap
figure();
subplot(1, 3, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['Zscore-OnlyPrior-',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

% Plot the second heatmap
subplot(1, 3, 2);
imagesc(averagedGrid_2_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['Zscore-OnlyPretone-',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

% Plot the third heatmap
subplot(1, 3, 3);
imagesc(averagedGrid_3_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['ZDiff',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%% Plot the Mean Coherhence  

% Replace NaNs with zero
Mean_OnlyPrior_averagedGrid(isnan(Mean_OnlyPrior_averagedGrid))     = 0;
Mean_OnlyPretone_averagedGrid(isnan(Mean_OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
globalMin = 0;
globalMax = .35;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Mean_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Mean_OnlyPretone_averagedGrid));

% Plot the first heatmap
figure();
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['Mean-OnlyPrior-',  Frequency_Band, '-', Animal]);
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
title(['Mean-OnlyPretone-',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

%% Plot the std Coherhence  

% Replace NaNs with zero
std_OnlyPrior_averagedGrid(isnan(std_OnlyPrior_averagedGrid))     = 0;
std_OnlyPretone_averagedGrid(isnan(std_OnlyPretone_averagedGrid)) = 0;

% Determine the global color range
globalMin = 0;
globalMax = .35;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, std_OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, std_OnlyPretone_averagedGrid));

% Plot the first heatmap
figure();
subplot(1, 2, 1);
imagesc(averagedGrid_1_clamped);
colormap(jet(256)); % 256 instead of default 64
clim([globalMin, globalMax]);
colorbar;
xlabel('AC Channels');
ylabel('PFC Channels');
title(['std-OnlyPrior-',  Frequency_Band, '-', Animal]);
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
title(['std-OnlyPretone-',  Frequency_Band, '-', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

%% plot the average full coherence spectra across sessions with std on them 


% Initialize an empty array to store all trials
Coherence_OnlyPrior_all_trials = [];

% Iterate through each cell in preCueOnset_PSI_grid
for i = 1:numel(Mean_OnlyPrior_Coh_grid)
    % Concatenate the trials from the current cell
    Coherence_OnlyPrior_all_trials = [Coherence_OnlyPrior_all_trials; Mean_OnlyPrior_Coh_grid{i}];
end

% Initialize an empty array to store all trials
Coherence_OnlyPretone_all_trials = [];

% Iterate through each cell in testToneOnset_PSI_grid
for b = 1:numel(Mean_OnlyPretone_Coh_grid)
    % Concatenate the trials from the current cell
    Coherence_OnlyPretone_all_trials = [Coherence_OnlyPretone_all_trials; Mean_OnlyPretone_Coh_grid{b}];
end

% Compute mean across trials
mean_data_Prior = mean(Coherence_OnlyPrior_all_trials, 1);
mean_data_Pretone = mean(Coherence_OnlyPretone_all_trials, 1);

% Calculate standard deviation (for error bars)
std_Prior = std(Coherence_OnlyPrior_all_trials, 0, 1);
std_Pretone = std(Coherence_OnlyPretone_all_trials, 0, 1);

% Plotting
freq = 1:101;  % Assuming 101 frequencies

figure;
hold on;

% Plot mean data with error bars using standard deviation
errorbar(freq, mean_data_Prior, std_Prior, 'b', 'LineWidth', 1.5);
errorbar(freq, mean_data_Pretone, std_Pretone, 'r', 'LineWidth', 1.5);

% Customize plot labels and legend
xlabel('Frequency');
ylabel('Coherence');
title('Coherence Comparison: Prior vs Pretone');
legend('Prior', 'Pretone', 'Location', 'best');

% Set x-axis limits and plot aesthetics
xlim([4 100]);  % Adjust x-axis limits to 4:100
set(gca, 'FontSize', 12);
grid on;
box on;
hold off;

%% plot the average full coherence spectra across sessions with 95% confidence interval around them.


% Initialize an empty array to store all trials
Coherence_OnlyPrior_all_trials = [];

% Iterate through each cell in preCueOnset_PSI_grid
for i = 1:numel(Mean_OnlyPrior_Coh_grid)
    % Concatenate the trials from the current cell
    Coherence_OnlyPrior_all_trials = [Coherence_OnlyPrior_all_trials; Mean_OnlyPrior_Coh_grid{i}];
end

% Initialize an empty array to store all trials
Coherence_OnlyPretone_all_trials = [];

% Iterate through each cell in testToneOnset_PSI_grid
for b = 1:numel(Mean_OnlyPretone_Coh_grid)
    % Concatenate the trials from the current cell
    Coherence_OnlyPretone_all_trials = [Coherence_OnlyPretone_all_trials; Mean_OnlyPretone_Coh_grid{b}];
end

% Compute mean and standard error across trials
mean_data_Prior = mean(Coherence_OnlyPrior_all_trials, 1);
mean_data_Pretone = mean(Coherence_OnlyPretone_all_trials, 1);

% Calculate standard error
SE_Prior = std(Coherence_OnlyPrior_all_trials, 0, 1) / sqrt(size(Coherence_OnlyPrior_all_trials, 1));
SE_Pretone = std(Coherence_OnlyPretone_all_trials, 0, 1) / sqrt(size(Coherence_OnlyPretone_all_trials, 1));

% Calculate 95% confidence intervals
CI_Prior = tinv(0.975, size(Coherence_OnlyPrior_all_trials, 1) - 1) * SE_Prior;
CI_Pretone = tinv(0.975, size(Coherence_OnlyPretone_all_trials, 1) - 1) * SE_Pretone;

% Plotting
freq = 1:101;  % Assuming 101 frequencies

figure;
hold on;

% Plot mean data with error bars
errorbar(freq, mean_data_Prior, CI_Prior, 'b', 'LineWidth', 1.5);
errorbar(freq, mean_data_Pretone, CI_Pretone, 'r', 'LineWidth', 1.5);

% Customize plot labels and legend
xlabel('Frequency');
ylabel('Coherence');
title('Coherence Comparison: Prior vs Pretone');
legend('Prior', 'Pretone', 'Location', 'best');

% Set plot aesthetics
set(gca, 'FontSize', 12);
grid on;
box on;
hold off;


% %% Band-wise collect all values for each channel pair, and collapse into column
% 
% if strcmp(Frequency_Band,'theta') == 1
% 
% % Extract columns 5 to 8 and reshape into a single column
% OnlyPrior_theta = Coherence_OnlyPrior_all_trials(:, 5:8);
% OnlyPrior_theta = OnlyPrior_theta(:);
% 
% OnlyPretone_theta = Coherence_OnlyPretone_all_trials(:, 5:8);
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
% OnlyPretone_alpha = Coherence_OnlyPretone_all_trials(:, 9:14);
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
% OnlyPretone_alpha = Coherence_OnlyPretone_all_trials(:, 15:30);
% OnlyPretone_alpha = OnlyPretone_alpha(:);
% 
% OnlyPrior_PFC_2_AC = OnlyPrior_alpha(OnlyPrior_alpha > 0);
% OnlyPrior_AC_2_PFC = OnlyPrior_alpha(OnlyPrior_alpha < 0);
% 
% OnlyPretone_PFC_2_AC = OnlyPretone_alpha(OnlyPretone_alpha > 0);
% OnlyPretone_AC_2_PFC = OnlyPretone_alpha(OnlyPretone_alpha < 0);
% 
% end

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

