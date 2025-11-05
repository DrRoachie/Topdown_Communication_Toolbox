%% Step 3: Organize PSI values by channel pair. 

Condition         = 'OnlyPretone';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'beta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals    = {'MrM'};                     % Options: 'MrCassius' and/or 'MrM'
rootdir    = 'C:\Users\Corey Roach\Documents\00_DATA\2025_19_12_PSI_Epoch_Analysis_OnlyPretone';                                 
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


%% Step 9: Plot with Means + 95% CI, Legend, and High-Resolution Export

abs_preCue = abs(preCue_averagedGrid);
abs_testTone = abs(testTone_averagedGrid);

% Mask out zeros (placeholder values)
mask_nonzero_pre = preCue_averagedGrid ~= 0;
mask_nonzero_test = testTone_averagedGrid ~= 0;

% Extract directional values
PFC_2_AC_pre     = abs_preCue(mask_nonzero_pre & preCue_averagedGrid > 0);
AC_2_PFC_pre     = abs_preCue(mask_nonzero_pre & preCue_averagedGrid < 0);
PFC_2_AC_test    = abs_testTone(mask_nonzero_test & testTone_averagedGrid > 0);
AC_2_PFC_test    = abs_testTone(mask_nonzero_test & testTone_averagedGrid < 0);

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 11, 8.5]); hold on;

group_data = {PFC_2_AC_pre, AC_2_PFC_pre, PFC_2_AC_test, AC_2_PFC_test};
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
         'XTickLabel', {'LED','Target Tone'}, ...
         'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold', ...
         'YGrid', 'off');
ylabel('Mean |PSI Value|', 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');

% Format y-axis ticks (no scientific notation)
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%.3f';

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

PSI_values = [PFC_2_AC_pre(:); AC_2_PFC_pre(:); PFC_2_AC_test(:); AC_2_PFC_test(:)];
Epoch      = [repmat({'PreCue'}, numel(PFC_2_AC_pre) + numel(AC_2_PFC_pre), 1);
              repmat({'TestTone'}, numel(PFC_2_AC_test) + numel(AC_2_PFC_test), 1)];
Direction  = [repmat({'PFC→AC'}, numel(PFC_2_AC_pre), 1);
              repmat({'AC→PFC'}, numel(AC_2_PFC_pre), 1);
              repmat({'PFC→AC'}, numel(PFC_2_AC_test), 1);
              repmat({'AC→PFC'}, numel(AC_2_PFC_test), 1)];

ranked_PSI = tiedrank(PSI_values);

[p, tbl, stats, terms] = anovan(ranked_PSI, ...
                                {Epoch, Direction}, ...
                                'model', 2, ...
                                'varnames', {'Epoch','Direction'}, ...
                                'display', 'on');

anova_table = cell2table(tbl(2:end,:), 'VariableNames', tbl(1,:));
disp(anova_table);

%% Step 11: Follow-up Pairwise Rank-Sum Tests

% ----- Claim 1: PreCue > TestTone (all directions combined) -----
PreCue_all   = [PFC_2_AC_pre(:); AC_2_PFC_pre(:)];
TestTone_all = [PFC_2_AC_test(:); AC_2_PFC_test(:)];
[p_epoch, ~, stats_epoch] = ranksum(PreCue_all, TestTone_all);
fprintf('\nClaim 1 — PreCue vs TestTone:\n');
fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
    p_epoch, stats_epoch.zval, numel(PreCue_all), numel(TestTone_all));

% ----- Claim 2: Directional difference within TestTone -----
[p_testtone_dir, ~, stats_testtone_dir] = ranksum(PFC_2_AC_test, AC_2_PFC_test);
fprintf('\nClaim 2 — TestTone: PFC→AC vs AC→PFC:\n');
fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
    p_testtone_dir, stats_testtone_dir.zval, ...
    numel(PFC_2_AC_test), numel(AC_2_PFC_test));

% ----- Claim 3: Directional difference within PreCue -----
[p_precue_dir, ~, stats_precue_dir] = ranksum(PFC_2_AC_pre, AC_2_PFC_pre);
fprintf('\nClaim 3 — PreCue: PFC→AC vs AC→PFC:\n');
fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
    p_precue_dir, stats_precue_dir.zval, ...
    numel(PFC_2_AC_pre), numel(AC_2_PFC_pre));



%% All data point statistics

% Initialize an empty array to store all trials
% preCue_all_trials = [];
% 
% % Iterate through each cell in preCueOnset_PSI_grid
% for i = 1:numel(preCueOnset_PSI_grid)
%     % Concatenate the trials from the current cell
%     preCue_all_trials = [preCue_all_trials; preCueOnset_PSI_grid{i}];
% end
% 
% % Initialize an empty array to store all trials
% testTone_all_trials = [];
% 
% % Iterate through each cell in testToneOnset_PSI_grid
% for b = 1:numel(testToneOnset_PSI_grid)
%     % Concatenate the trials from the current cell
%     testTone_all_trials = [testTone_all_trials; testToneOnset_PSI_grid{i}];
% end


% Band-wise collect all values for each channel pair, and collapse into column

% if strcmp(Frequency_Band,'theta') == 1
% 
% % Extract columns 5 to 8 and reshape into a single column
% preCue_theta = preCue_all_trials(:, 5:8);
% preCue_theta = preCue_theta(:);
% 
% testTone_theta = testTone_all_trials(:, 5:8);
% testTone_theta = testTone_theta(:);
% 
% preCue_PFC_2_AC = preCue_theta(preCue_theta > 0);
% preCue_AC_2_PFC = preCue_theta(preCue_theta < 0);
% 
% testTone_PFC_2_AC = testTone_theta(testTone_theta > 0);
% testTone_AC_2_PFC = testTone_theta(testTone_theta < 0);
% 
% end

% %% plotting histogram 
% 
% %[~, bins] = hist([preCue_PFC_2_AC; abs(preCue_AC_2_PFC); testTone_PFC_2_AC; abs(testTone_AC_2_PFC)], 100);
% load('C:\Users\Corey Roach\Documents\00_DATA\bins.mat')
% 
% % Compute histograms for each dataset
% posCounts_1 = hist(preCue_PFC_2_AC, bins);
% negCounts_1 = hist(abs(preCue_AC_2_PFC), bins);
% posCounts_2 = hist(testTone_PFC_2_AC, bins);
% negCounts_2 = hist(abs(testTone_AC_2_PFC), bins);
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
% % xlim([-.15 .15])
% axis([-.03 .03 -100 100]) 
% leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
% set(leg, 'FontSize', 12); % Reduce legend font size
% 
% % Improve aesthetics
% set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
% box on;
% hold off;
% 
% figure(3)
% 
% hold on 
% bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
% bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
% 
% % Labels and Legend
% xlabel('PSI Value', 'FontSize', 18);
% ylabel('Proportion of Sites', 'FontSize', 18);
% axis([-.03 .03 -100 100]) 
% 
% % Set legend and adjust its font size
% leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
% set(leg, 'FontSize', 12); % Reduce legend font size
% 
% % Improve aesthetics
% set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
% box on;
% hold off;

% %% this runs a statistical test on all values significant values n = > 1000 versus the test of the average where n = number of channel pairs {DO NOT DELETE, maybe useful}
% 
% % Assume these four variables are already column vectors of theta-band PSI
% % preCue_PFC_2_AC, preCue_AC_2_PFC, testTone_PFC_2_AC, testTone_AC_2_PFC
% 
% % Take absolute values
% preCue_PFC_2_AC_abs   = abs(preCue_PFC_2_AC);
% preCue_AC_2_PFC_abs   = abs(preCue_AC_2_PFC);
% testTone_PFC_2_AC_abs = abs(testTone_PFC_2_AC);
% testTone_AC_2_PFC_abs = abs(testTone_AC_2_PFC);
% 
% % Combine data into a single vector
% allData = [preCue_PFC_2_AC_abs; preCue_AC_2_PFC_abs; ...
%            testTone_PFC_2_AC_abs; testTone_AC_2_PFC_abs];
% 
% % Create grouping variable
% group = [repmat({'preCue PFC→AC'}, length(preCue_PFC_2_AC_abs), 1);
%          repmat({'preCue AC→PFC'}, length(preCue_AC_2_PFC_abs), 1);
%          repmat({'testTone PFC→AC'}, length(testTone_PFC_2_AC_abs), 1);
%          repmat({'testTone AC→PFC'}, length(testTone_AC_2_PFC_abs), 1)];
% 
% % --- Run one-way ANOVA ---
% [p,tbl,stats] = anova1(allData, group, 'off'); % 'off' suppresses default figure
% disp('One-way ANOVA results:');
% disp(tbl);
% 
% % Post-hoc multiple comparisons (Tukey-Kramer)
% c = multcompare(stats, 'CType','tukey-kramer');
% disp('Post-hoc comparisons (Tukey-Kramer):');
% disp(c);
% 
% % --- Compute means and 95% confidence intervals ---
% dataVars = {preCue_PFC_2_AC_abs, preCue_AC_2_PFC_abs, ...
%             testTone_PFC_2_AC_abs, testTone_AC_2_PFC_abs};
% 
% means = cellfun(@mean, dataVars);
% 
% % Compute 95% CI
% conf95 = @(x) 1.96 * std(x)/sqrt(length(x));
% ci95 = cellfun(conf95, dataVars);
% 
% % --- Plot fat bars with 95% CI, large square figure ---
% figure('Color','w','Units','inches','Position',[1 1 8 8],'PaperPositionMode','auto'); hold on;
% 
% barWidth = 0.8; % fat bars
% 
% % Assign colors: PFC→AC = red, AC→PFC = blue (RGB 0-255 scaled to 0-1)
% barColors = [222 30 38;    % red
%              71 133 197;   % blue
%              222 30 38;    % red
%              71 133 197]/255;
% 
% b = bar(1:4, means, barWidth, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth',1.5);
% for k = 1:4
%     b.CData(k,:) = barColors(k,:);
% end
% 
% % Add 95% CI error bars
% errorbar(1:4, means, ci95, 'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',16);
% 
% % Aesthetics
% set(gca, 'XTick', [1.5 3.5], ... % position labels between bars
%          'XTickLabel', {'LED','Target Tone'}, ...
%          'FontName','Arial', 'FontSize', 16, 'FontWeight','bold');
% 
% ylabel('Absolute Theta PSI', 'FontName','Arial','FontSize',16,'FontWeight','bold');
% 
% ylim([0 max(means+ci95)*1.15]);
% 
% 
% % Reduce whitespace
% ax = gca;
% ax.Position = [0.15 0.15 0.75 0.75]; % [left bottom width height]
% 
% %title('Absolute Theta-band PSI', 'FontSize',18,'FontWeight','bold','FontName','Arial');
% 
% % Make the figure scalable and high-resolution when saving
% set(gcf,'Renderer','painters');

