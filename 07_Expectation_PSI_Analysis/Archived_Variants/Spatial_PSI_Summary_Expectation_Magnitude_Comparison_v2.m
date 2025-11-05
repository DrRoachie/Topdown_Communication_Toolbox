%% Expectation Analysis
% Updated: CR 9/30/25
% Overall, this code compares the observed phase slope index values during a target tone between two expectation conditions (OnlyPrior and OnlyPretone). 
% I tierank the values and use a standard anova logic. After the omnibus test, i run pairwise wilcoxin tests to make claims about my data.

%% Step 1: Book-Keeping 

Frequency_Band    = 'theta';
Behavior          = 'correct';
%Congruency        = [];
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius'};                                                   % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_24_PSI_Expectation_Analysis_All_Congruency';                                 
sessions          = dir(fullfile(rootdir, '19*'));


%% Step 2: Loop through all sessions and sort PSI values by channel-pair 

% Define the PFC and AC channel labels and generate temporary grids
channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
OnlyPrior_PSI_grid     = cell(20, 20); 
OnlyPretone_PSI_grid   = cell(20, 20);

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

%% Step 3: Average selected bands to output a single PSI value per channel pair 

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

%% Step 4: Generate arrays with all significant interactions 

% Initialize an empty array to store all trials
OnlyPrior_all_sig_interactions = [];
OnlyPretone_all_sig_interactions = [];

% Iterate through each cell in OnlyPrior_PSI_grid
for i = 1:numel(OnlyPrior_PSI_grid)
    % Concatenate the trials from the current cell
    OnlyPrior_all_sig_interactions = [OnlyPrior_all_sig_interactions; OnlyPrior_PSI_grid{i}];
end

% Iterate through each cell in OnlyPretone_PSI_grid
for b = 1:numel(OnlyPretone_PSI_grid)
    % Concatenate the trials from the current cell
    OnlyPretone_all_sig_interactions = [OnlyPretone_all_sig_interactions; OnlyPretone_PSI_grid{b}];
end

% Band-wise collect all values for each channel pair, and collapse into column

if strcmp(Frequency_Band,'theta') == 1

% Extract columns 5 to 8 and reshape into a single column
OnlyPrior_all_sig_interactions_theta = OnlyPrior_all_sig_interactions(:, 5:8);
OnlyPrior_all_sig_interactions_theta = OnlyPrior_all_sig_interactions_theta(:);

OnlyPretone_all_sig_interactions_theta = OnlyPretone_all_sig_interactions(:, 5:8);
OnlyPretone_all_sig_interactions_theta = OnlyPretone_all_sig_interactions_theta(:);

OnlyPrior_all_sig_interactions_PFC_2_AC = OnlyPrior_all_sig_interactions_theta(OnlyPrior_all_sig_interactions_theta > 0);
OnlyPrior_all_sig_interactions_AC_2_PFC = abs(OnlyPrior_all_sig_interactions_theta(OnlyPrior_all_sig_interactions_theta < 0));

OnlyPretone_all_sig_interactions_PFC_2_AC = OnlyPretone_all_sig_interactions_theta(OnlyPretone_all_sig_interactions_theta > 0);
OnlyPretone_all_sig_interactions_AC_2_PFC = abs(OnlyPretone_all_sig_interactions_theta(OnlyPretone_all_sig_interactions_theta < 0));

end
%% Step 5: Plot the average PSI value of each channel pair [HEAT MAP]

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


%% Step 6: Plot with Mean of Means + 95% CI [SUMMARY OF HEAT MAP]

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


%% Step 7: Plot Mean of all Sig Interactions 

figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 11, 8.5]); hold on;

group_data = {OnlyPrior_all_sig_interactions_PFC_2_AC, OnlyPrior_all_sig_interactions_AC_2_PFC, OnlyPretone_all_sig_interactions_PFC_2_AC, OnlyPretone_all_sig_interactions_AC_2_PFC};
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

%% Step 8: Test the influence of expectation on mean of mean PSI values extracted from heat maps 

% --- Build per-pair summaries from your existing grids (theta cols 5:8) ---
pair_means_prior   = nan(20,20);
pair_counts_prior  = zeros(20,20);
pair_posfrac_prior = nan(20,20);  % fraction of positive samples per pair

for pfcIdx = 1:20
  for acIdx = 1:20
    X = OnlyPrior_PSI_grid{pfcIdx, acIdx};
    if ~isempty(X)
      v = X(:,5:8); v = v(:);
      pair_means_prior(pfcIdx,acIdx)  = mean(v, 'omitnan');    % signed per-pair mean
      pair_counts_prior(pfcIdx,acIdx) = sum(~isnan(v));         % sample count per pair
      pair_posfrac_prior(pfcIdx,acIdx)= mean(v>0, 'omitnan');   % within-pair sign mix
    end
  end
end

% Flatten to vectors (keep valid pairs only)
m  = pair_means_prior(:);
n  = pair_counts_prior(:);
pf = pair_posfrac_prior(:);
valid = ~isnan(m) & (n>0);
m  = m(valid); n = n(valid); pf = pf(valid);

% --- Unweighted vs pooled (weighted) means over *magnitude* like your bars ---
unweighted_mean_mag = mean(abs(m), 'omitnan');
pooled_mean_mag     = sum(abs(m).*n)/sum(n);

fprintf('Unweighted mean(|per-pair mean|): %.6f\n', unweighted_mean_mag);
fprintf('Pooled (n-weighted) mean(|per-pair mean|): %.6f\n', pooled_mean_mag);
fprintf('Ratio pooled / unweighted: %.3f\n', pooled_mean_mag / unweighted_mean_mag);

% --- Is n correlated with magnitude? (this is the big tell) ---
[r, p] = corr(n, abs(m), 'type','Spearman');  % robust to nonlinearity
fprintf('Spearman corr between n and |per-pair mean|: r=%.3f, p=%.3g\n', r, p);

% --- How much do mixed signs damp per-pair means? ---
% If many pairs have pf ~ 0.5, per-pair signed mean shrinks while pooled-by-sign stays large.
[rpf, ppf] = corr(pf, abs(m), 'type','Spearman');
fprintf('Sign-mix (positive fraction) vs |per-pair mean|: r=%.3f, p=%.3g\n', rpf, ppf);

% --- Sanity check: replicate your Step 7 pooled mean directly from samples ---
all_samples = [];
for pfcIdx = 1:20
  for acIdx = 1:20
    X = OnlyPrior_PSI_grid{pfcIdx, acIdx};
    if ~isempty(X), all_samples = [all_samples; X(:,5:8)]; end %#ok<AGROW>
  end
end
pooled_from_samples = mean(all_samples(:), 'omitnan');  % signed pooled
fprintf('Signed pooled mean (all samples): %.6f\n', pooled_from_samples);


% %% Step 10: Nonparametric 2×2 ANOVA on Ranked PSI Values
% 
% PSI_values       = [PFC_2_AC_OnlyPrior(:); AC_2_PFC_OnlyPrior(:); PFC_2_AC_OnlyPretone(:); AC_2_PFC_OnlyPretone(:)];
% Expectation      = [repmat({'OnlyPrior'}, numel(PFC_2_AC_OnlyPrior) + numel(AC_2_PFC_OnlyPrior), 1);
%                     repmat({'OnlyPretone'}, numel(PFC_2_AC_OnlyPretone) + numel(AC_2_PFC_OnlyPretone), 1)];
% Direction        = [repmat({'PFC→AC'}, numel(PFC_2_AC_OnlyPrior), 1);
%                     repmat({'AC→PFC'}, numel(AC_2_PFC_OnlyPrior), 1);
%                     repmat({'PFC→AC'}, numel(PFC_2_AC_OnlyPretone), 1);
%                     repmat({'AC→PFC'}, numel(AC_2_PFC_OnlyPretone), 1)];
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
% OnlyPrior_all   = [PFC_2_AC_OnlyPrior(:); AC_2_PFC_OnlyPrior(:)];
% OnlyPretone_all = [PFC_2_AC_OnlyPretone(:); AC_2_PFC_OnlyPretone(:)];
% [p_expectation, ~, stats_expectation] = ranksum(OnlyPrior_all, OnlyPretone_all);
% fprintf('\nClaim 1 — OnlyPrior vs OnlyPretone:\n');
% fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
%     p_expectation, stats_expectation.zval, numel(OnlyPrior_all), numel(OnlyPretone_all));
% 
% % ----- Claim 2: Directional difference within OnlyPretone -----
% [p_OnlyPretone_dir, ~, stats_OnlyPretone_dir] = ranksum(PFC_2_AC_OnlyPretone, AC_2_PFC_OnlyPretone);
% fprintf('\nClaim 2 — OnlyPretone: PFC→AC vs AC→PFC:\n');
% fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
%     p_OnlyPretone_dir, stats_OnlyPretone_dir.zval, ...
%     numel(PFC_2_AC_OnlyPretone), numel(AC_2_PFC_OnlyPretone));
% 
% % ----- Claim 3: Directional difference within PreCue -----
% [p_OnlyPrior_dir, ~, stats_OnlyPrior_dir] = ranksum(PFC_2_AC_OnlyPrior, AC_2_PFC_OnlyPrior);
% fprintf('\nClaim 3 — OnlyPrior: PFC→AC vs AC→PFC:\n');
% fprintf('p = %.4f, Z = %.2f, n1 = %d, n2 = %d\n', ...
%     p_OnlyPrior_dir, stats_OnlyPrior_dir.zval, ...
%     numel(PFC_2_AC_OnlyPrior), numel(AC_2_PFC_OnlyPrior));
% 
% %% Step 13: Collective Magnitude 
% 
% x = OnlyPrior_all(:);
% y = OnlyPretone_all(:);
% 
% % Wilcoxon signed-rank test (paired, nonparametric)
% [p, h, stats] = signrank(x, y);   % stats.signedrank = W+, stats.zval = z
% 
% % Number of non-zero pairs actually used
% nz = sum((x - y) ~= 0);
% 
% % Rank-biserial correlation (paired) as effect size
% T = nz*(nz+1)/2;            % total rank sum
% Wplus = stats.signedrank;   % sum of positive ranks
% r_rb = 2*Wplus/T - 1;       % matched-pairs rank-biserial
% 
% % Optional: report medians & IQRs
% mx = mean(x,'omitnan'); my = median(y,'omitnan');
% iqx = iqr(x); iqy = iqr(y);
% 
% fprintf('Wilcoxon signed-rank: W+=%.0f, z=%.3f, p=%.4g, r_rb=%.3f\n', Wplus, stats.zval, p, r_rb);
% fprintf('OnlyPrior: median=%.3g, IQR=%.3g | OnlyPretone: median=%.3g, IQR=%.3g\n', mx, iqx, my, iqy);
% 
% %%
% 
% 
% % Colors
% colPrior   = [0    0.4470 0.7410];   % blue
% colPretone = [0.8500 0.3250 0.0980]; % orange
% 
% % Summary stats (medians & IQR)
% mx  = median(x,'omitnan');  my  = median(y,'omitnan');
% q1x = quantile(x,0.25);     q3x = quantile(x,0.75);
% q1y = quantile(y,0.25);     q3y = quantile(y,0.75);
% 
% % Asymmetric IQR errorbars around the median
% neg = [mx - q1x, my - q1y];
% pos = [q3x - mx, q3y - my];
% 
% % Figure setup
% figure('Color','w'); clf
% set(groot, 'defaultAxesFontName', 'Arial', ...
%            'defaultTextFontName', 'Arial', ...
%            'defaultAxesFontSize', 16, ...
%            'defaultTextFontSize', 16);
% 
% tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
% 
% %% A) Box plots
% ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
% boxchart(ones(size(x)), x, 'BoxFaceColor', colPrior,   'WhiskerLineColor', colPrior,   'BoxFaceAlpha', 0.35);
% boxchart(2*ones(size(y)), y, 'BoxFaceColor', colPretone,'WhiskerLineColor', colPretone,'BoxFaceAlpha', 0.35);
% set(ax1,'XTick',[1 2],'XTickLabel',{'OnlyPrior','OnlyPretone'});
% ylabel(ax1,'Value');
% title(ax1,'Paired summary (box plots)');
% 
% %% B) Bar plot of medians with IQR error bars
% ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
% bar(1, mx, 'FaceColor', colPrior,   'EdgeColor','none'); 
% bar(2, my, 'FaceColor', colPretone, 'EdgeColor','none');
% 
% % Asymmetric IQR around the median
% errorbar([1 2], [mx my], neg, pos, 'k', 'LineStyle','none', 'LineWidth', 1.5, 'CapSize', 10);
% 
% set(ax2,'XTick',[1 2],'XTickLabel',{'OnlyPrior','OnlyPretone'});
% ylabel(ax2,'Median ± IQR');
% title(ax2, sprintf('Wilcoxon signed-rank: W^+=%d, z=%.3f, p=%.3g, r_{rb}=%.3f', ...
%                    Wplus, stats.zval, p, r_rb));
% 
% % Optional: y-limits, grid, and saving
% % ylim(ax1, [min([x;y]) max([x;y])]);
% % ylim(ax2, [min([x;y]) max([x;y])]);
% % grid(ax1,'on'); grid(ax2,'on');
% % exportgraphics(gcf, 'Paired_SummaryPlots.png', 'Resolution', 300);
% 
% %% Re-analysis: comparing all significant interactions
% % This section is another version of the statistical analysis. Rather than taking the average of each channel pair, and then averaging them together (a mean of means) resulting in N = 400,
% % we average across all significant interactions (N = number of significant
% % interactions coming out of the trial-balanced PSI comparison) 
% 
% 
% %%
% 
% % Combine the two arrays
% OnlyPrior_combined_array   = [OnlyPrior_all_sig_interactions_PFC_2_AC; OnlyPrior_all_sig_interactions_AC_2_PFC];
% OnlyPrior_combined_array   = abs(OnlyPrior_combined_array);
% 
% OnlyPretone_combined_array = [OnlyPretone_all_sig_interactions_PFC_2_AC; OnlyPretone_all_sig_interactions_AC_2_PFC];
% OnlyPretone_combined_array = abs(OnlyPretone_combined_array);
% 
% % Wilcoxon signed-rank test (paired, nonparametric)
% [p, h, stats] = signrank(OnlyPrior_combined_array, OnlyPretone_combined_array);   % stats.signedrank = W+, stats.zval = z
% 
% % Number of non-zero pairs actually used
% nz = sum((OnlyPrior_combined_array - OnlyPretone_combined_array) ~= 0);
% 
% % Rank-biserial correlation (paired) as effect size
% T = nz*(nz+1)/2;            % total rank sum
% Wplus = stats.signedrank;   % sum of positive ranks
% r_rb = 2*Wplus/T - 1;       % matched-pairs rank-biserial
% 
% % Optional: report medians & IQRs
% mx = mean(OnlyPrior_combined_array,'omitnan'); my = median(OnlyPretone_combined_array,'omitnan');
% iqx = iqr(OnlyPrior_combined_array); iqy = iqr(OnlyPretone_combined_array);
% 
% fprintf('Wilcoxon signed-rank: W+=%.0f, z=%.3f, p=%.4g, r_rb=%.3f\n', Wplus, stats.zval, p, r_rb);
% fprintf('OnlyPrior: median=%.3g, IQR=%.3g | OnlyPretone: median=%.3g, IQR=%.3g\n', mx, iqx, my, iqy);
% 
% % Colors
% colPrior   = [0    0.4470 0.7410];   % blue
% colPretone = [0.8500 0.3250 0.0980]; % orange
% 
% % Summary stats (medians & IQR)
% mx  = median(x,'omitnan');  my  = median(y,'omitnan');
% q1x = quantile(x,0.25);     q3x = quantile(x,0.75);
% q1y = quantile(y,0.25);     q3y = quantile(y,0.75);
% 
% % Asymmetric IQR errorbars around the median
% neg = [mx - q1x, my - q1y];
% pos = [q3x - mx, q3y - my];
% 
% % Figure setup
% figure('Color','w'); clf
% set(groot, 'defaultAxesFontName', 'Arial', ...
%            'defaultTextFontName', 'Arial', ...
%            'defaultAxesFontSize', 16, ...
%            'defaultTextFontSize', 16);
% 
% tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
% 
% %% A) Box plots
% ax1 = nexttile; hold(ax1,'on'); box(ax1,'on');
% boxchart(ones(size(OnlyPrior_combined_array)), OnlyPrior_combined_array, 'BoxFaceColor', colPrior,   'WhiskerLineColor', colPrior,   'BoxFaceAlpha', 0.35);
% boxchart(2*ones(size(OnlyPretone_combined_array)), OnlyPretone_combined_array, 'BoxFaceColor', colPretone,'WhiskerLineColor', colPretone,'BoxFaceAlpha', 0.35);
% set(ax1,'XTick',[1 2],'XTickLabel',{'OnlyPrior','OnlyPretone'});
% ylabel(ax1,'Value');
% title(ax1,'Paired summary (box plots)');
% 
% %% B) Bar plot of medians with IQR error bars
% ax2 = nexttile; hold(ax2,'on'); box(ax2,'on');
% bar(1, mx, 'FaceColor', colPrior,   'EdgeColor','none'); 
% bar(2, my, 'FaceColor', colPretone, 'EdgeColor','none');
% 
% % Asymmetric IQR around the median
% errorbar([1 2], [mx my], neg, pos, 'k', 'LineStyle','none', 'LineWidth', 1.5, 'CapSize', 10);
% 
% set(ax2,'XTick',[1 2],'XTickLabel',{'OnlyPrior','OnlyPretone'});
% ylabel(ax2,'Median ± IQR');
% title(ax2, sprintf('Wilcoxon signed-rank: W^+=%d, z=%.3f, p=%.3g, r_{rb}=%.3f', ...
%                    Wplus, stats.zval, p, r_rb));
% 
% % Optional: y-limits, grid, and saving
% % ylim(ax1, [min([x;y]) max([x;y])]);
% % ylim(ax2, [min([x;y]) max([x;y])]);
% % grid(ax1,'on'); grid(ax2,'on');
% % exportgraphics(gcf, 'Paired_SummaryPlots.png', 'Resolution', 300);