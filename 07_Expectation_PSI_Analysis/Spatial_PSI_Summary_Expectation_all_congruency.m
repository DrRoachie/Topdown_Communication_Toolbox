%% Expectation Analysis
% Updated: CR 9/30/25
% Overall, this code compares the observed phase slope index values during a target tone between two expectation conditions (OnlyPrior and OnlyPretone). 
% I tierank the values and use a standard anova logic. After the omnibus test, i run pairwise wilcoxin tests to make claims about my data.

%% Step 3: Organize PSI values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
Behavior          = 'correct';
%Congruency        = [];
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius'};                                                   % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_24_PSI_Expectation_Analysis_All_Congruency';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
OnlyPrior_PSI_grid     = cell(20, 20); 
OnlyPretone_PSI_grid   = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', 'all-congruency', '_', Frequency_Band,'_PSI.mat']));
            PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_','OnlyPretone_', Behavior,'_', 'all-congruency', '_', Frequency_Band,'_PSI.mat']));

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
globalMin = -0.02;
globalMax = 0.02;

% Clamp values to the global color range
averagedGrid_1_clamped = max(globalMin, min(globalMax, OnlyPrior_averagedGrid));
averagedGrid_2_clamped = max(globalMin, min(globalMax, OnlyPretone_averagedGrid));

% averagedGrid_1_clamped = OnlyPrior_averagedGrid;
% averagedGrid_2_clamped = OnlyPretone_averagedGrid;

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
title(['OnlyPrior-',  Frequency_Band, '-', Animal]);
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
title(['OnlyPretone-',  Frequency_Band, '-', Animal]);
xticks(1:24);
yticks(1:24);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%% Plot with Means + 95% CI, Legend, and High-Resolution Export

abs_OnlyPrior   = abs(OnlyPrior_averagedGrid);
abs_OnlyPretone = abs(OnlyPretone_averagedGrid);

% Mask out zeros (placeholder values)
mask_nonzero_OnlyPrior   = OnlyPrior_averagedGrid ~= 0;
mask_nonzero_OnlyPretone = OnlyPretone_averagedGrid ~= 0;

% Extract directional values
PFC_2_AC_OnlyPrior     = abs_OnlyPrior(mask_nonzero_OnlyPrior & OnlyPrior_averagedGrid > 0);
AC_2_PFC_OnlyPrior     = abs_OnlyPrior(mask_nonzero_OnlyPrior & OnlyPrior_averagedGrid < 0);
PFC_2_AC_OnlyPretone    = abs_OnlyPretone(mask_nonzero_OnlyPretone & OnlyPretone_averagedGrid > 0);
AC_2_PFC_OnlyPretone    = abs_OnlyPretone(mask_nonzero_OnlyPretone & OnlyPretone_averagedGrid < 0);

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 11, 8.5]); hold on;

group_data = {PFC_2_AC_OnlyPrior, AC_2_PFC_OnlyPrior, PFC_2_AC_OnlyPretone, AC_2_PFC_OnlyPretone};
means = cellfun(@mean, group_data);
SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
CIs   = SEMs * 1.96;  % 95% CI

bar_colors = [222 30 38; 71 133 197; 222 30 38; 71 133 197]/255;

% Bar plot
for i = 1:4
    b = bar(i, means(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
    hold on;
end

% Add CI error bars
errorbar(1:4, means, CIs, 'k.', 'LineWidth', 1.5);

% Format axes
set(gca, 'XTick', [1.5 3.5], ...
         'XTickLabel', {'OnlyPrior','OnlyPretone'}, ...
         'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold', ...
         'YGrid', 'off');
ylabel('Mean |PSI Value|', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');

% Format y-axis ticks (no scientific notation)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%.4f';

% Add legend in top-left corner (no box)
legend({'PFC→AC', 'AC→PFC'}, ...
       'Location', 'northwest', ...
       'FontSize', 14, ...
       'FontName', 'Arial', ...
       'Box', 'off');


box off
axis tight

% Optional: Export high-DPI image
exportgraphics(gcf, 'PSI_BarPlot_Labeled_HighRes.png', 'Resolution', 600);


%% Step 10: Nonparametric 2×2 ANOVA on Ranked PSI Values

PSI_values       = [PFC_2_AC_OnlyPrior(:); AC_2_PFC_OnlyPrior(:); PFC_2_AC_OnlyPretone(:); AC_2_PFC_OnlyPretone(:)];
Expectation      = [repmat({'OnlyPrior'}, numel(PFC_2_AC_OnlyPrior) + numel(AC_2_PFC_OnlyPrior), 1);
                    repmat({'OnlyPretone'}, numel(PFC_2_AC_OnlyPretone) + numel(AC_2_PFC_OnlyPretone), 1)];
Direction        = [repmat({'PFC→AC'}, numel(PFC_2_AC_OnlyPrior), 1);
                    repmat({'AC→PFC'}, numel(AC_2_PFC_OnlyPrior), 1);
                    repmat({'PFC→AC'}, numel(PFC_2_AC_OnlyPretone), 1);
                    repmat({'AC→PFC'}, numel(AC_2_PFC_OnlyPretone), 1)];

ranked_PSI = tiedrank(PSI_values);

[p, tbl, stats, terms]   = anovan(ranked_PSI, ...
                                {Expectation, Direction}, ...
                                'model', 2, ...
                                'varnames', {'Expectationion','Direction'}, ...
                                'display', 'on');

anova_table = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
disp(anova_table);

%% Step 11: Follow-up Pairwise Rank-Sum Tests

% ----- Claim 1: OnlyPrior > OnlyPretone (all directions combined) -----
OnlyPrior_all   = [PFC_2_AC_OnlyPrior(:); AC_2_PFC_OnlyPrior(:)];
OnlyPretone_all = [PFC_2_AC_OnlyPretone(:); AC_2_PFC_OnlyPretone(:)];
[p_expectation, ~, stats_expectation] = ranksum(OnlyPrior_all, OnlyPretone_all);
fprintf('\nClaim 1 — OnlyPrior vs OnlyPretone:\n');
fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
    p_expectation, stats_expectation.zval, numel(OnlyPrior_all), numel(OnlyPretone_all));

% ----- Claim 2: Directional difference within OnlyPretone -----
[p_OnlyPretone_dir, ~, stats_OnlyPretone_dir] = ranksum(PFC_2_AC_OnlyPretone, AC_2_PFC_OnlyPretone);
fprintf('\nClaim 2 — OnlyPretone: PFC→AC vs AC→PFC:\n');
fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
    p_OnlyPretone_dir, stats_OnlyPretone_dir.zval, ...
    numel(PFC_2_AC_OnlyPretone), numel(AC_2_PFC_OnlyPretone));

% ----- Claim 3: Directional difference within PreCue -----
[p_OnlyPrior_dir, ~, stats_OnlyPrior_dir] = ranksum(PFC_2_AC_OnlyPrior, AC_2_PFC_OnlyPrior);
fprintf('\nClaim 3 — OnlyPrior: PFC→AC vs AC→PFC:\n');
fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
    p_OnlyPrior_dir, stats_OnlyPrior_dir.zval, ...
    numel(PFC_2_AC_OnlyPrior), numel(AC_2_PFC_OnlyPrior));