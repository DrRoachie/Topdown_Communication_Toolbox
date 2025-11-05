% Organize PSI values by channel pair. 

Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'beta';
Behavior          = 'correct';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrM'};                                                   % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\Coh_Correct_Expectation_Visualization';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
OnlyPrior_Coh_grid     = cell(20, 20); 
OnlyPretone_Coh_grid   = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            Coh_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_',  Frequency_Band,'_Coh.mat']));
            Coh_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_','OnlyPretone_', Behavior,'_',  Frequency_Band,'_Coh.mat']));

            if ~isempty(Coh_OnlyPrior_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(Coh_OnlyPrior_file)

        Coh_OnlyPrior_filename   = Coh_OnlyPrior_file(k).name;
        Coh_OnlyPretone_filename = Coh_OnlyPretone_file(k).name;

        temp_OnlyPrior_dir   = fullfile(session_path, Coh_OnlyPrior_filename);
        temp_OnlyPretone_dir = fullfile(session_path, Coh_OnlyPretone_filename);

        % Load the .mat file safely
         load(temp_OnlyPrior_dir,'Coh_OnlyPrior' );
         load(temp_OnlyPretone_dir,'Coh_OnlyPretone');

       OnlyPrior_Matrix     = Coh_OnlyPrior.cohspctrm; 
       OnlyPrior_labels     = Coh_OnlyPrior.labelcmb;         
       OnlyPretone_Matrix   = Coh_OnlyPretone.cohspctrm; 
       OnlyPretone_labels   = Coh_OnlyPretone.labelcmb;       

      end    

           % Populate OnlyPrior_Coh_grid with data 

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
                    if isempty(OnlyPrior_Coh_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    OnlyPrior_Coh_grid{pfcIdx, acIdx} = OnlyPrior_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    OnlyPrior_Coh_grid{pfcIdx, acIdx} = [OnlyPrior_Coh_grid{pfcIdx, acIdx}; OnlyPrior_Matrix(rowIdx_1, :)];
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
                    if isempty(OnlyPretone_Coh_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    OnlyPretone_Coh_grid{pfcIdx, acIdx} = OnlyPretone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    OnlyPretone_Coh_grid{pfcIdx, acIdx} = [OnlyPretone_Coh_grid{pfcIdx, acIdx}; OnlyPretone_Matrix(rowIdx_2, :)];
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
        if ~isempty(OnlyPrior_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPrior_Coh_grid{pfcIdx, acIdx};

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
        if ~isempty(OnlyPrior_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPrior_Coh_grid{pfcIdx, acIdx};

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
        if ~isempty(OnlyPretone_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPretone_Coh_grid{pfcIdx, acIdx};

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
        if ~isempty(OnlyPretone_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = OnlyPretone_Coh_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            OnlyPretone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end



%% Step 5: Collapse across channels and sessions

% --- OnlyPrior ---
OnlyPrior_all = [];  % Initialize empty matrix

for pfcIdx = 1:20
    for acIdx = 1:20
        if ~isempty(OnlyPrior_Coh_grid{pfcIdx, acIdx})
            OnlyPrior_all = [OnlyPrior_all; OnlyPrior_Coh_grid{pfcIdx, acIdx}]; %#ok<AGROW>
        end
    end
end

% Compute descriptive statistics
OnlyPrior_mean  = mean(OnlyPrior_all, 1, 'omitnan');
OnlyPrior_std   = std(OnlyPrior_all, 0, 1, 'omitnan');
n_prior         = size(OnlyPrior_all, 1);
OnlyPrior_CI95  = 1.96 * OnlyPrior_std ./ sqrt(n_prior);  % 95% CI

% --- OnlyPretone ---
OnlyPretone_all = [];  % Initialize empty matrix

for pfcIdx = 1:20
    for acIdx = 1:20
        if ~isempty(OnlyPretone_Coh_grid{pfcIdx, acIdx})
            OnlyPretone_all = [OnlyPretone_all; OnlyPretone_Coh_grid{pfcIdx, acIdx}]; %#ok<AGROW>
        end
    end
end


%% Step 6: Box plot for band-specific coherence



    % --- Select columns corresponding to the band ---
    switch Frequency_Band
        case 'theta'
            freqCols = 4:8;      % columns 4-8 for theta
            yLims = [0 0.4];
        case 'beta'
            freqCols = 15:30;    % columns 15-30 for beta
            yLims = [0 0.3];
    end

    % --- Collapse across selected frequency bins ---
    OnlyPrior_collapsed  = mean(OnlyPrior_all(:, freqCols), 2);   % mean across selected bins per trial
    OnlyPretone_collapsed = mean(OnlyPretone_all(:, freqCols), 2);

    % --- Combine data for boxplot ---
    allData = [OnlyPrior_collapsed; OnlyPretone_collapsed];
    group   = [ones(size(OnlyPrior_collapsed)); 2*ones(size(OnlyPretone_collapsed))];

    % --- Colors for printing ---
    boxColors = [0 0.6 0.6;    % turquoise
                 0.9 0.7 0];  % yellow

    % --- Create figure with custom size (inches) ---
    figure('Units','inches','Position',[1 1 8 8]); % 6x6 inches square
    hold on;

    % Create boxplot
    h = boxplot(allData, group, 'Widths',0.6, 'Colors','k', 'Whisker',1.5, 'Symbol',''); 
    % Remove default outlier symbols for cleaner print

    % Customize box colors
    boxes = findobj(h, 'Tag', 'Box');
    for j = 1:length(boxes)
        patch(get(boxes(j),'XData'), get(boxes(j),'YData'), boxColors(j,:), ...
            'FaceAlpha',0.5, 'EdgeColor','k', 'LineWidth',1.5);
    end

    % Median lines bold
    medians = findobj(h, 'Tag', 'Median');
    set(medians, 'Color', 'k', 'LineWidth', 2);

    % Whiskers bold
    whiskers = findobj(h, 'Tag', 'Whisker');
    set(whiskers, 'Color', 'k', 'LineWidth', 1.5);

    % Labels and formatting
    set(gca, 'XTickLabel', {'Informative','NonInformative'}, 'FontWeight','bold', 'FontSize', 12);
    ylabel('Coherence', 'FontWeight','bold', 'FontSize', 12);
   % title([Frequency_Band, ' Coherence (Collapsed Across Band)'], 'FontWeight','bold', 'FontSize', 18);

    ylim(yLims);    % Set y-axis limits
    xlim([0.5 2.5]); % Remove white space before first box
    set(gca,'Position',[0.18 0.18 0.72 0.72]); % Make figure more square/tight
    box on;

    hold off;

    % --- Save as high-res TIFF for publication ---
    print(gcf, [Frequency_Band '_Coherence_BoxPlot.tif'], '-dtiff', '-r600');




%% Step 7: Statistical comparison (OnlyPrior vs OnlyPretone, theta & beta bands)

if strcmp(Frequency_Band,'theta') == 1

    % Define theta column indices
    thetaCols = 5:8;

    % Collapse across frequencies to get single value per row (channel pair)
    OnlyPrior_theta_avg   = mean(OnlyPrior_all(:, thetaCols), 2, 'omitnan');   
    OnlyPretone_theta_avg = mean(OnlyPretone_all(:, thetaCols), 2, 'omitnan');

    % Ensure matching number of rows
    minRows = min(length(OnlyPrior_theta_avg), length(OnlyPretone_theta_avg));
    OnlyPrior_theta_avg   = OnlyPrior_theta_avg(1:minRows);
    OnlyPretone_theta_avg = OnlyPretone_theta_avg(1:minRows);

    % Wilcoxon signed-rank test
    [pval,~,stats] = signrank(OnlyPrior_theta_avg, OnlyPretone_theta_avg);

    % Compute descriptive statistics (median + 25th/75th percentiles)
    median_prior   = median(OnlyPrior_theta_avg);
    median_pretone = median(OnlyPretone_theta_avg);
    pct25_prior    = prctile(OnlyPrior_theta_avg, 25);
    pct75_prior    = prctile(OnlyPrior_theta_avg, 75);
    pct25_pretone  = prctile(OnlyPretone_theta_avg, 25);
    pct75_pretone  = prctile(OnlyPretone_theta_avg, 75);

    % Display results
    fprintf('Theta band (4-8 Hz) comparison:\n');
    fprintf('OnlyPrior: Median = %.4f, IQR = %.4f - %.4f\n', median_prior, pct25_prior, pct75_prior);
    fprintf('OnlyPretone: Median = %.4f, IQR = %.4f - %.4f\n', median_pretone, pct25_pretone, pct75_pretone);
    fprintf('Wilcoxon signed-rank test: p = %.4f, z = %.2f\n', pval, stats.zval);

end

if strcmp(Frequency_Band,'beta') == 1

    % Define beta column indices
    betaCols = 15:30;

    % Collapse across frequencies
    OnlyPrior_beta_avg   = mean(OnlyPrior_all(:, betaCols), 2, 'omitnan');   
    OnlyPretone_beta_avg = mean(OnlyPretone_all(:, betaCols), 2, 'omitnan');

    % Ensure matching number of rows
    minRows = min(length(OnlyPrior_beta_avg), length(OnlyPretone_beta_avg));
    OnlyPrior_beta_avg   = OnlyPrior_beta_avg(1:minRows);
    OnlyPretone_beta_avg = OnlyPretone_beta_avg(1:minRows);

    % Wilcoxon signed-rank test
    [pval,~,stats] = signrank(OnlyPrior_beta_avg, OnlyPretone_beta_avg);

    % Compute descriptive statistics
    median_prior   = median(OnlyPrior_beta_avg);
    median_pretone = median(OnlyPretone_beta_avg);
    pct25_prior    = prctile(OnlyPrior_beta_avg, 25);
    pct75_prior    = prctile(OnlyPrior_beta_avg, 75);
    pct25_pretone  = prctile(OnlyPretone_beta_avg, 25);
    pct75_pretone  = prctile(OnlyPretone_beta_avg, 75);

    % Display results
    fprintf('Beta band (15-30 Hz) comparison:\n');
    fprintf('OnlyPrior: Median = %.4f, IQR = %.4f - %.4f\n', median_prior, pct25_prior, pct75_prior);
    fprintf('OnlyPretone: Median = %.4f, IQR = %.4f - %.4f\n', median_pretone, pct25_pretone, pct75_pretone);
    fprintf('Wilcoxon signed-rank test: p = %.4f, z = %.2f\n', pval, stats.zval);

end










