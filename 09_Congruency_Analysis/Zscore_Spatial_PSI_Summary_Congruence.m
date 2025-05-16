%% Step 3: Organize Coherence values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
Behavior          = 'correct';
% Congruency        = 'congruent';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius', 'MrM'};                                                     % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\Congruency_Analysis';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Make str. list out of 'animals' for figure titles
% Join the animal names with commas or another separator
animalStr = strjoin(animals, ', ');


% Define the PFC and AC channel labels and generate temporary grids

channels                  = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Mean_Congruent_PSI_grid   = cell(20, 20); 
Mean_Incongruent_PSI_grid = cell(20, 20);
Z_Congruent_PSI_grid      = cell(20, 20); 
Z_Incongruent_PSI_grid    = cell(20, 20);
ZDiff_PSI_grid            = cell(20, 20); 


%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            PSI_Congruent_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_',  Behavior,'_', 'congruent',   '_', Frequency_Band,'_PSI.mat']));
            PSI_Incongruent_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_',  Behavior,'_', 'incongruent', '_', Frequency_Band,'_PSI.mat']));

            if ~isempty(PSI_Congruent_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(PSI_Congruent_file)

        PSI_Congruent_filename   = PSI_Congruent_file(k).name;
        PSI_Incongruent_filename = PSI_Incongruent_file(k).name;

        temp_Congruent_dir   = fullfile(session_path, PSI_Congruent_filename);
        temp_Incongruent_dir = fullfile(session_path, PSI_Incongruent_filename);

        % Load the .mat file safely
        load(temp_Congruent_dir,'PSI_congruent' );
        load(temp_Incongruent_dir,'PSI_incongruent');

        Mean_Congruent_Matrix   = PSI_congruent.psispctrm; 
        Mean_Incongruent_Matrix = PSI_incongruent.psispctrm; 
        Z_Congruent_Matrix      = PSI_congruent.psispctrm; 
        Z_Congruent_labels      = PSI_congruent.labelcmb;         
        Z_Incongruent_Matrix    = PSI_incongruent.psispctrm; 
        Z_Incongruent_labels    = PSI_incongruent.labelcmb;  

        % z-score the data 

        % Step 1: Pool all values from both matrices
      
        pooledData = [Z_Congruent_Matrix(:); Z_Incongruent_Matrix(:)];

        % Step 2: Compute global mean and std
        globalMean = mean(pooledData);
        globalStd  = std(pooledData);


        % Step 3: Z-score both matrices
        Z_Congruent_Matrix   = (Z_Congruent_Matrix   - globalMean) / globalStd;
        Z_Incongruent_Matrix = (Z_Incongruent_Matrix - globalMean) / globalStd;
        zDiff_Matrix         = Z_Congruent_Matrix - Z_Incongruent_Matrix;

      end    

           % Populate Z_Congruent_Coherence_grid with data 

           for rowIdx_1 = 1:size(Z_Congruent_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_Congruent_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_Congruent_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Z_Congruent_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_Congruent_PSI_grid{pfcIdx, acIdx} = Mean_Congruent_Matrix(rowIdx_1, :);
                    Z_Congruent_PSI_grid{pfcIdx, acIdx} = Z_Congruent_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Mean_Congruent_PSI_grid{pfcIdx, acIdx} = [Mean_Congruent_PSI_grid{pfcIdx, acIdx}; Mean_Congruent_Matrix(rowIdx_1, :)];
                    Z_Congruent_PSI_grid{pfcIdx, acIdx} = [Z_Congruent_PSI_grid{pfcIdx, acIdx}; Z_Congruent_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate Incongruent_PSI_grid with data 

           for rowIdx_2 = 1:size(Z_Incongruent_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_Incongruent_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_Incongruent_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Z_Incongruent_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_Incongruent_PSI_grid{pfcIdx, acIdx} = Mean_Incongruent_Matrix(rowIdx_2, :);
                    Z_Incongruent_PSI_grid{pfcIdx, acIdx} = Z_Incongruent_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Mean_Incongruent_PSI_grid{pfcIdx, acIdx} = [Mean_Incongruent_PSI_grid{pfcIdx, acIdx}; Mean_Incongruent_Matrix(rowIdx_2, :)];
                    Z_Incongruent_PSI_grid{pfcIdx, acIdx} = [Z_Incongruent_PSI_grid{pfcIdx, acIdx}; Z_Incongruent_Matrix(rowIdx_2, :)];
                    end                 
                end
           end
           
           % Populate ZDiff_Coherence_PSI_grid with data 

           for rowIdx_1 = 1:size(zDiff_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_Congruent_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_Congruent_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(ZDiff_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    ZDiff_PSI_grid{pfcIdx, acIdx} = zDiff_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    ZDiff_PSI_grid{pfcIdx, acIdx} = [ZDiff_PSI_grid{pfcIdx, acIdx}; zDiff_Matrix(rowIdx_1, :)];
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
Z_Congruent_averagedGrid     = nan(20, 20); % Use NaN to handle empty grid positions
Mean_Congruent_averagedGrid  = nan(20, 20); % Use NaN to handle empty grid positions
std_Congruent_averagedGrid   = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_Congruent_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition    = Z_Congruent_PSI_grid{pfcIdx, acIdx};
            Mean_dataposition = Mean_Congruent_PSI_grid{pfcIdx, acIdx};
            
            % Compute the average of the selected columns
            Z_columnAverage    = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage = mean(Mean_dataposition(:, columnRange), 'all', 'omitnan');
            std_columnAverage  = std(Mean_dataposition(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_Congruent_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage;
            Mean_Congruent_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage;
            std_Congruent_averagedGrid(pfcIdx, acIdx)  = std_columnAverage;

        end
    end
end

%% Average selected bands to output a single PSI value per channel pair (OnlyPretone)

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
Z_Incongruent_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions
Mean_Incongruent_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions
std_Incongruent_averagedGrid  = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_Incongruent_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition    = Z_Incongruent_PSI_grid{pfcIdx, acIdx};
            Mean_dataposition = Mean_Incongruent_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage    = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage = mean(Mean_dataposition(:, columnRange), 'all', 'omitnan');
            std_columnAverage = std(Mean_dataposition(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_Incongruent_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage;
            Mean_Incongruent_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage;
            std_Incongruent_averagedGrid(pfcIdx, acIdx)  = std_columnAverage;
        end
    end
end

%% Average selected bands to output a single PSI value per channel pair (ZDiff)

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
        if ~isempty(ZDiff_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition = ZDiff_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            ZDiff_averagedGrid(pfcIdx, acIdx) = Z_columnAverage;
        end
    end
end


%% Plot the Z-scored PSI  

% Replace NaNs with zero
Z_Congruent_averagedGrid(isnan(Z_Congruent_averagedGrid)) = 0;
Z_Incongruent_averagedGrid(isnan(Z_Incongruent_averagedGrid)) = 0;
ZDiff_averagedGrid(isnan(Z_Congruent_averagedGrid)) = 0;

% Determine the global color range
globalMin = -0.4;
globalMax = 0.4;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, Z_Congruent_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, Z_Incongruent_averagedGrid));
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
title(['Zscore-Congruent-',  Frequency_Band, '-', animalStr]);
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
title(['Zscore-Incongruent-',  Frequency_Band, '-', animalStr]);
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
title(['ZDiff',  Frequency_Band, '-', animalStr]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%% plot the average full PSI spectra across sessions with 95% confidence interval around them.

% Initialize an empty array to store all trials
PSI_Congruent_all_trials = [];

% Iterate through each cell in preCueOnset_PSI_grid
for i = 1:numel(Z_Congruent_PSI_grid)
    % Concatenate the trials from the current cell
    PSI_Congruent_all_trials = [PSI_Congruent_all_trials; Z_Congruent_PSI_grid{i}];
end

% Initialize an empty array to store all trials
PSI_Incongruent_all_trials = [];

% Iterate through each cell in testToneOnset_PSI_grid
for b = 1:numel(Z_Incongruent_PSI_grid)
    % Concatenate the trials from the current cell
    PSI_Incongruent_all_trials = [PSI_Incongruent_all_trials; Z_Incongruent_PSI_grid{b}];
end

% 

% Compute mean and standard error across trials
mean_data_Congruent = mean(PSI_Congruent_all_trials, 1);
mean_data_Incongruent = mean(PSI_Incongruent_all_trials, 1);

% Apply absolute values to the means
mean_data_Congruent = abs(mean_data_Congruent);
mean_data_Incongruent = abs(mean_data_Incongruent);

% Calculate standard error
SE_Congruent = std(PSI_Congruent_all_trials, 0, 1) / sqrt(size(PSI_Congruent_all_trials, 1));
SE_Incongruent = std(PSI_Incongruent_all_trials, 0, 1) / sqrt(size(PSI_Incongruent_all_trials, 1));

% Calculate 95% confidence intervals
CI_Congruent = tinv(0.975, size(PSI_Congruent_all_trials, 1) - 1) * SE_Congruent;
CI_Incongruent = tinv(0.975, size(PSI_Incongruent_all_trials, 1) - 1) * SE_Incongruent;

% Plotting
freq = 1:101;  % Assuming 101 frequencies

figure;
hold on;

% Plot mean data with error bars (absolute values)
errorbar(freq, mean_data_Congruent, CI_Congruent, 'b', 'LineWidth', 1.5);
errorbar(freq, mean_data_Incongruent, CI_Incongruent, 'r', 'LineWidth', 1.5);

% Customize plot labels and legend
xlabel('Frequency');
ylabel('|PSI|');
xlim([4 35])
title('PSI Comparison: Congruent vs Incongruent');
legend('Congruent', 'Incongruent', 'Location', 'best');

% Set plot aesthetics
set(gca, 'FontSize', 12);
grid on;
box on;
hold off;

%%

% Count true zeros in OnlyPrior_theta
num_zeros_Prior = sum(PSI_Congruent_all_trials == 0);

% Count true zeros in OnlyPretone_theta
num_zeros_Pretone = sum(PSI_Incongruent_all_trials == 0);

%% Band-wise collect all values for each channel pair, and collapse into column

if strcmp(Frequency_Band,'theta') == 1

    % Extract columns 5 to 8 and reshape into a single column
    Congruent_theta   = PSI_Congruent_all_trials(:, 5:8);
    Congruent_theta   = Congruent_theta(:);
    Incongruent_theta = PSI_Incongruent_all_trials(:, 5:8);
    Incongruent_theta = Incongruent_theta(:);

    % % Split by direction
    Congruent_PFC_2_AC   = Congruent_theta(Congruent_theta > 0);
    Incongruent_PFC_2_AC = Incongruent_theta(Incongruent_theta > 0);
    Congruent_AC_2_PFC   = Congruent_theta(Congruent_theta < 0);
    Incongruent_AC_2_PFC = Incongruent_theta(Incongruent_theta < 0);

    % Define bin edges from -6 to +6 (centered)
    binEdges = linspace(-10, 10, 50);

    % Plot all histograms on the same figure
    figure;
    hold on;

    histogram(Congruent_theta, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], ...
              'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Blue
    histogram(Incongruent_theta, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], ...
              'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Light red

    xlabel('Phase Slope Index (PSI)');
    ylabel('Count');
    xlim([-10 10])
    legend({'Congruent', 'Incongruent'});
    title('Z-Scored PSI Values');
    box on;
    set(gca, 'FontSize', 14);
    hold off;

end

%% plotting histogram 

%[~, bins] = hist([OnlyPrior_PFC_2_AC; abs(OnlyPrior_AC_2_PFC); OnlyPretone_PFC_2_AC; abs(OnlyPretone_AC_2_PFC)], 100);
edges = linspace(0, 10, 10);                % 100 bins from 0 to 0.1
bins  = edges(1:end-1) + diff(edges)/2;       % Bin centers


% Compute histograms for each dataset
Congruent_PFC_2_AC_counts   = hist(Congruent_PFC_2_AC, bins);
Congruent_AC_2_PFC_counts   = hist(abs(Congruent_AC_2_PFC), bins);
Incongruent_PFC_2_AC_counts = hist(Incongruent_PFC_2_AC, bins);
Incongruent_AC_2_PFC_counts = hist(abs(Incongruent_AC_2_PFC), bins);

% Normalize counts to percentages
Congruent_PFC_2_AC_counts   = Congruent_PFC_2_AC_counts / sum(Congruent_PFC_2_AC_counts) * 100;
Congruent_AC_2_PFC_counts   = Congruent_AC_2_PFC_counts / sum(Congruent_AC_2_PFC_counts) * 100;
Incongruent_PFC_2_AC_counts = Incongruent_PFC_2_AC_counts / sum(Incongruent_PFC_2_AC_counts) * 100;
Incongruent_AC_2_PFC_counts = Incongruent_AC_2_PFC_counts / sum(Incongruent_AC_2_PFC_counts) * 100;

% Plot histograms

figure();

hold on;

% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   Congruent_PFC_2_AC_counts,   'FaceColor', [1, 0.6, 0.6],   'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins, -Congruent_AC_2_PFC_counts,   'FaceColor', [0.6, 0.6, 1],   'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
bar(bins,   Incongruent_PFC_2_AC_counts, 'FaceColor', [0.8, 0.1, 0.1], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins, -Incongruent_AC_2_PFC_counts, 'FaceColor', [0.2, 0.4, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
xlim([-10 10])
%axis([-.025 .025 -65 65]) 
% Adjust legend to include all conditions
leg = legend({'Congruent PFC→AC', 'Congruent AC→PFC', 'Incongruent PFC→AC', 'Incongruent AC→PFC'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size
% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;



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

% %% Jamie Kuminski Core/Belt Analysis
% 
% Core_OnlyPrior_averagedGrid   = OnlyPrior_averagedGrid;
% Core_OnlyPretone_averagedGrid = OnlyPretone_averagedGrid;
% save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Core_Spatial_PSI_Values.mat'), 'Core_OnlyPrior_averagedGrid', 'Core_OnlyPretone_averagedGrid');
% 
% Belt_OnlyPrior_averagedGrid   = OnlyPrior_averagedGrid;
% Belt_OnlyPretone_averagedGrid = OnlyPretone_averagedGrid;
% save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Belt_Spatial_PSI_Values.mat'), 'Belt_OnlyPrior_averagedGrid', 'Belt_OnlyPretone_averagedGrid');
% 
% Subtract the session-wise means from each other and plot as a histogram.
% 
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
% 
% 
% 
% %
% 
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

