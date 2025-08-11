%% SpecBalanced RAW PSI & Coherence Visualization

% This code visualizes connectivity between PFC and AC based on SpecBalanced channel pairs
% using both Phase Slope Index (PSI) and Coherence metrics.
% Visualize:
%    - RAW PSI and Coherence as heatmaps (Prior, Pretone)
%    - Spectra plots across frequency with 95% CI bands for Prior and Pretone (PSI and Coherence)
%    - Histograms of PSI values (signed → no abs) to visualize PFC→AC vs AC→PFC directionality
%    - Histograms of Coherence values (prior vs pretone) 


% Input conditions 

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
Mean_OnlyPrior_PSI_grid   = cell(20, 20); 
Mean_OnlyPretone_PSI_grid = cell(20, 20);


% Adding grids for coherence 
Mean_OnlyPrior_Coh_grid   = cell(20, 20); 
Mean_OnlyPretone_Coh_grid = cell(20, 20);


%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));
            PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));

            Coh_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', Congruency, '_', Frequency_Band,'_Coh.mat']));
            Coh_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_Coh.mat']));

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

       % Load all files
        load(temp_OnlyPrior_dir,'PSI_OnlyPrior' );
        load(temp_OnlyPretone_dir,'PSI_OnlyPretone');
        load(fullfile(session_path, Coh_OnlyPrior_file.name), 'Coh_OnlyPrior');
        load(fullfile(session_path, Coh_OnlyPretone_file.name), 'Coh_OnlyPretone');
        
        Mean_OnlyPrior_Matrix   = PSI_OnlyPrior.psispctrm; 
        Mean_OnlyPretone_Matrix = PSI_OnlyPretone.psispctrm; 
        
        % Define grid labels
        channels   = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
        pfc_labels = arrayfun(@(x) sprintf('PFC_ch%02d', x), 3:22, 'UniformOutput', false);
        ac_labels  = arrayfun(@(x) sprintf('AC_ch%02d',  x), 3:22, 'UniformOutput', false);
        
        % === LOOP THROUGH PSI ONLYPRIOR ===
        for rowIdx_1 = 1:size(Mean_OnlyPrior_Matrix, 1)
            pfcLabel = strrep(PSI_OnlyPrior.labelcmb{rowIdx_1, 1}, '*', '');
            acLabel  = strrep(PSI_OnlyPrior.labelcmb{rowIdx_1, 2}, '*', '');
        
            pfcIdx = find(strcmp(pfc_labels, pfcLabel));
            acIdx  = find(strcmp(ac_labels,  acLabel));
        
            if ~isempty(pfcIdx) && ~isempty(acIdx)
                if isempty(Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx})
                    Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx} = Mean_OnlyPrior_Matrix(rowIdx_1, :);
                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = Coh_OnlyPrior.cohspctrm(rowIdx_1, :);
                else
                    Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx}; Mean_OnlyPrior_Matrix(rowIdx_1, :)];
                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx}; Coh_OnlyPrior.cohspctrm(rowIdx_1, :)];
                end
            end
        end


           
            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx} = Mean_OnlyPrior_Matrix(rowIdx_1, :);

                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = Coh_OnlyPrior.cohspctrm(rowIdx_1, :);

                    else
                    % Append the new row of data to the existing data
                    Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx}; Mean_OnlyPrior_Matrix(rowIdx_1, :)];

                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx}; Coh_OnlyPrior.cohspctrm(rowIdx_1, :)];

                    end

                end
           end

           % Populate OnlyPretone_PSI_grid with data 

           for rowIdx_2 = 1:size(Mean_OnlyPretone_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = PSI_OnlyPretone.labelcmb{rowIdx_2, 1};  
            acLabel  = PSI_OnlyPretone.labelcmb{rowIdx_2, 2};  
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx} = Mean_OnlyPretone_Matrix(rowIdx_2, :);
                  
                    Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx} = Coh_OnlyPretone.cohspctrm(rowIdx_2, :);
                  
                    else
                    % Append the new row of data to the existing data
                    Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx} = [Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx}; Mean_OnlyPretone_Matrix(rowIdx_2, :)];
                   
                    Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx}; Coh_OnlyPretone.cohspctrm(rowIdx_2, :)];
                   
                    end                 
                end
           end
           
           % Populate ZDiff_Coherence_PSI_grid with data 

           % for rowIdx_1 = 1:size(zDiff_Matrix, 1)
           % 
           %  % Extract channel pair from label
           %  pfcLabel = Z_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
           %  acLabel  = Z_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
           % 
           %  % Remove the 'D' component and isolate the channel labels
           %  pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
           %  acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
           % 
           %  % Find indices in the grid
           %  pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
           %  acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));
           % 
           % 
           % 
           %      end
           end
            

         end
     end



%% Plot RAW heatmaps

% Plot RAW PSI heatmaps: Prior + Pretone side by side

% % Replace NaNs with 0 
% Mean_OnlyPrior_averagedGrid(isnan(Mean_OnlyPrior_averagedGrid)) = 0;
% Mean_OnlyPretone_averagedGrid(isnan(Mean_OnlyPretone_averagedGrid)) = 0;

% Convert each cell to a scalar mean, safely
Mean_OnlyPrior_Coh_averagedGrid = nan(20, 20);  % initialize grid with NaNs

for r = 1:20
    for c = 1:20
        data = Mean_OnlyPrior_Coh_grid{r, c};
        if ~isempty(data)
            Mean_OnlyPrior_Coh_averagedGrid(r, c) = mean(data(:), 'omitnan');
        end
    end
end

Mean_OnlyPretone_Coh_averagedGrid = nan(20, 20);  % initialize grid with NaNs

for r = 1:20
    for c = 1:20
        data = Mean_OnlyPretone_Coh_grid{r, c};  % ← use the *raw cell grid*
        if ~isempty(data)
            Mean_OnlyPretone_Coh_averagedGrid(r, c) = mean(data(:), 'omitnan');
        end
    end
end

Mean_OnlyPrior_PSI_averagedGrid = nan(20, 20);  % initialize grid with NaNs

for r = 1:20
    for c = 1:20
        data = Mean_OnlyPrior_PSI_grid{r, c};  % ← use the *raw cell grid*        
        if ~isempty(data)
            Mean_OnlyPrior_PSI_averagedGrid(r, c) = mean(data(:), 'omitnan');
        end
    end
end

Mean_OnlyPretone_PSI_averagedGrid = nan(20, 20);  % initialize grid with NaNs

for r = 1:20
    for c = 1:20
        data = Mean_OnlyPretone_PSI_grid{r, c};  % ← use the *raw cell grid*
        if ~isempty(data)
            Mean_OnlyPretone_PSI_averagedGrid(r, c) = mean(data(:), 'omitnan');
        end
    end
end


% % Average over rows in each cell
% Mean_OnlyPrior_averagedGrid   = cellfun(@(x) mean(x, 1, 'omitnan'), Mean_OnlyPrior_PSI_grid, 'UniformOutput', false);
% Mean_OnlyPretone_averagedGrid = cellfun(@(x) mean(x, 1, 'omitnan'), Mean_OnlyPretone_PSI_grid, 'UniformOutput', false);
% 
% % Convert to matrix (each cell is 1×nFreq → extract scalar mean for heatmap)
% Mean_OnlyPrior_averagedGrid   = cellfun(@(x) mean(x), Mean_OnlyPrior_averagedGrid, 'UniformOutput', false);
% Mean_OnlyPrior_averagedGrid   = cell2mat(reshape(Mean_OnlyPrior_averagedGrid, 20, 20));
% 
% Mean_OnlyPretone_averagedGrid = cellfun(@(x) mean(x), Mean_OnlyPretone_averagedGrid, 'UniformOutput', false);
% Mean_OnlyPretone_averagedGrid = cell2mat(reshape(Mean_OnlyPretone_averagedGrid, 20, 20));

% Count how many grid cells are empty
num_empty_cells = sum(cellfun(@isempty, Mean_OnlyPrior_PSI_grid), 'all');
fprintf('Empty grid cells in Mean_OnlyPrior_PSI_grid: %d out of 400\n', num_empty_cells);

 
globalMin_PSI = -0.005;
globalMax_PSI = 0.005;

figure;

% Subplot 1: Prior
subplot(1, 2, 1);
imagesc(Mean_OnlyPrior_averagedGrid);
colormap(jet(256));
clim([globalMin_PSI, globalMax_PSI]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPrior_PSI_averagedGrid);

xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW PSI - OnlyPrior - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;

% 
% %  plot boxes around individual NaN cells ---
% hold on;
% 
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
%             % Plot a black square around this cell
%             x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
%             y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
%             plot(x, y, 'k-', 'LineWidth', 1.5);
%         end
%     end
% end
% 
% hold off;


% Subplot 2: Pretone
subplot(1, 2, 2);
imagesc(Mean_OnlyPretone_averagedGrid);
colormap(jet(256));
clim([globalMin_PSI, globalMax_PSI]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPretone_PSI_averagedGrid);


xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW PSI - OnlyPretone - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


% %  plot boxes around individual NaN cells ---
% hold on;
% 
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
%             % Plot a black square around this cell
%             x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
%             y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
%             plot(x, y, 'k-', 'LineWidth', 1.5);
%         end
%     end
% end
% 
% hold off;

% Plot RAW Coherence heatmaps: Prior + Pretone side by side

% % Replace NaNs with 0
% Mean_OnlyPrior_Coh_averagedGrid(isnan(Mean_OnlyPrior_Coh_averagedGrid)) = 0;
% Mean_OnlyPretone_Coh_averagedGrid(isnan(Mean_OnlyPretone_Coh_averagedGrid)) = 0;

% Mean_OnlyPrior_Coh_averagedGrid   = cellfun(@(x) mean(x), Mean_OnlyPrior_Coh_grid, 'UniformOutput', false);
% Mean_OnlyPrior_Coh_averagedGrid   = cell2mat(reshape(Mean_OnlyPrior_Coh_averagedGrid, 20, 20));
% 
% Mean_OnlyPretone_Coh_averagedGrid = cellfun(@(x) mean(x), Mean_OnlyPretone_Coh_grid, 'UniformOutput', false);
% Mean_OnlyPretone_Coh_averagedGrid = cell2mat(reshape(Mean_OnlyPretone_Coh_averagedGrid, 20, 20));



globalMin_Coh = 0.005;
globalMax_Coh = max([Mean_OnlyPrior_Coh_averagedGrid(:); Mean_OnlyPretone_Coh_averagedGrid(:)]);

figure;

% Subplot 1: Prior
subplot(1, 2, 1);
imagesc(Mean_OnlyPrior_Coh_averagedGrid);
colormap(jet(256));
clim([globalMin_Coh, globalMax_Coh]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPrior_Coh_averagedGrid);


xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW Coherence - OnlyPrior - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;



% %  plot boxes around individual NaN cells ---
% hold on;
% 
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
%             % Plot a black square around this cell
%             x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
%             y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
%             plot(x, y, 'k-', 'LineWidth', 1.5);
%         end
%     end
% end
% 
% hold off;


% Subplot 2: Pretone
subplot(1, 2, 2);
imagesc(Mean_OnlyPretone_Coh_averagedGrid);
colormap(jet(256));
clim([globalMin_Coh, globalMax_Coh]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPretone_Coh_averagedGrid);


xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW Coherence - OnlyPretone - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


% %  plot boxes around individual NaN cells ---
% hold on;
% 
% for pfcIdx = 1:20
%     for acIdx = 1:20
%         if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
%             % Plot a black square around this cell
%             x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
%             y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
%             plot(x, y, 'k-', 'LineWidth', 1.5);
%         end
%     % end
% end

hold off;

%% Plot the average full PSI spectra across sessions with 95% confidence interval around them.

% Initialize an empty array to store all trials
PSI_OnlyPrior_all_trials = [];

% Iterate through each cell in preCueOnset_PSI_grid
for i = 1:numel(Z_OnlyPrior_PSI_grid)
    % Concatenate the trials from the current cell
    PSI_OnlyPrior_all_trials = [PSI_OnlyPrior_all_trials; Z_OnlyPrior_PSI_grid{i}];
end

% Initialize an empty array to store all trials
PSI_OnlyPretone_all_trials = [];

% Iterate through each cell in testToneOnset_PSI_grid
for b = 1:numel(Z_OnlyPretone_PSI_grid)
    % Concatenate the trials from the current cell
    PSI_OnlyPretone_all_trials = [PSI_OnlyPretone_all_trials; Z_OnlyPretone_PSI_grid{b}];
end

% 

% Compute mean and standard error across trials
mean_data_Prior = mean(PSI_OnlyPrior_all_trials, 1);
mean_data_Pretone = mean(PSI_OnlyPretone_all_trials, 1);

% Apply absolute values to the means
mean_data_Prior = abs(mean_data_Prior);
mean_data_Pretone = abs(mean_data_Pretone);

% Calculate standard error
SE_Prior = std(PSI_OnlyPrior_all_trials, 0, 1) / sqrt(size(PSI_OnlyPrior_all_trials, 1));
SE_Pretone = std(PSI_OnlyPretone_all_trials, 0, 1) / sqrt(size(PSI_OnlyPretone_all_trials, 1));

% Calculate 95% confidence intervals
CI_Prior = tinv(0.975, size(PSI_OnlyPrior_all_trials, 1) - 1) * SE_Prior;
CI_Pretone = tinv(0.975, size(PSI_OnlyPretone_all_trials, 1) - 1) * SE_Pretone;

% RAW PSI Spectra
freq = 1:101;
figure;
hold on;

errorbar(freq, mean(PSI_OnlyPrior_all_trials, 1), ...
         std(PSI_OnlyPrior_all_trials, 0, 1) / sqrt(size(PSI_OnlyPrior_all_trials,1)), ...
         'b', 'LineWidth', 1.5);

errorbar(freq, mean(PSI_OnlyPretone_all_trials, 1), ...
         std(PSI_OnlyPretone_all_trials, 0, 1) / sqrt(size(PSI_OnlyPretone_all_trials,1)), ...
         'r', 'LineWidth', 1.5);

xlabel('Frequency (Hz)');
ylabel('RAW PSI');
xlim([4 35]);
title(['RAW PSI Spectra: Prior vs Pretone']);
legend('Prior', 'Pretone');
grid on;
box on;
hold off;
% 
% RAW Coherence Spectra
figure;
hold on;

errorbar(freq, mean(Coh_OnlyPrior_all_trials, 1), ...
         std(Coh_OnlyPrior_all_trials, 0, 1) / sqrt(size(Coh_OnlyPrior_all_trials,1)), ...
         'b', 'LineWidth', 1.5);

errorbar(freq, mean(Coh_OnlyPretone_all_trials, 1), ...
         std(Coh_OnlyPretone_all_trials, 0, 1) / sqrt(size(Coh_OnlyPretone_all_trials,1)), ...
         'r', 'LineWidth', 1.5);

xlabel('Frequency (Hz)');
ylabel('RAW Coherence');
xlim([4 35]);
title(['RAW Coherence Spectra: Prior vs Pretone']);
legend('Prior', 'Pretone');
grid on;
box on;
hold off;


% Count true zeros in OnlyPrior_theta
num_zeros_Prior = sum(PSI_OnlyPrior_all_trials == 0);

% Count true zeros in OnlyPretone_theta
num_zeros_Pretone = sum(PSI_OnlyPretone_all_trials == 0);

% Band-wise collect all values for each channel pair, and collapse into column

%% RAW PSI Histogram
if strcmp(Frequency_Band,'theta') == 1
    OnlyPrior_theta = PSI_OnlyPrior_all_trials(:, 5:8);
    OnlyPretone_theta = PSI_OnlyPretone_all_trials(:, 5:8);

    OnlyPrior_theta = OnlyPrior_theta(:);
    OnlyPretone_theta = OnlyPretone_theta(:);

    binEdges = linspace(-10, 10, 50);

    figure;
    hold on;

    histogram(OnlyPrior_theta, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    histogram(OnlyPretone_theta, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    xlabel('RAW PSI Value');
    ylabel('Count');
    xlim([-10 10]);
    legend({'OnlyPrior', 'OnlyPretone'});
    title('RAW PSI: Prior vs Pretone');
    box on;
    set(gca, 'FontSize', 14);
    hold off;
end

% RAW Coherence Histogram
if strcmp(Frequency_Band,'theta') == 1
    OnlyPrior_Coh_theta = Coh_OnlyPrior_all_trials(:, 5:8);
    OnlyPretone_Coh_theta = Coh_OnlyPretone_all_trials(:, 5:8);

    OnlyPrior_Coh_theta = OnlyPrior_Coh_theta(:);
    OnlyPretone_Coh_theta = OnlyPretone_Coh_theta(:);

    binEdges = linspace(0, 1, 50); % Coherence is positive, range [0,1]

    figure;
    hold on;

    histogram(OnlyPrior_Coh_theta, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    histogram(OnlyPretone_Coh_theta, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    xlabel('RAW Coherence Value');
    ylabel('Count');
    xlim([0 1]);
    legend({'OnlyPrior', 'OnlyPretone'});
    title('RAW Coherence: Prior vs Pretone');
    box on;
    set(gca, 'FontSize', 14);
    hold off;
end


%% plotting histogram 

[~, bins] = hist([OnlyPrior_PFC_2_AC; abs(OnlyPrior_AC_2_PFC); OnlyPretone_PFC_2_AC; abs(OnlyPretone_AC_2_PFC)], 100);
edges = linspace(0, 10, 10);                % 100 bins from 0 to 0.1
bins  = edges(1:end-1) + diff(edges)/2;       % Bin centers


% Compute histograms for each dataset
OnlyPrior_PFC_2_AC_counts   = hist(OnlyPrior_PFC_2_AC, bins);
OnlyPrior_AC_2_PFC_counts   = hist(abs(OnlyPrior_AC_2_PFC), bins);
OnlyPretone_PFC_2_AC_counts = hist(OnlyPretone_PFC_2_AC, bins);
OnlyPretone_AC_2_PFC_counts = hist(abs(OnlyPretone_AC_2_PFC), bins);

% Normalize counts to percentages
OnlyPrior_PFC_2_AC_counts   = OnlyPrior_PFC_2_AC_counts / sum(OnlyPrior_PFC_2_AC_counts) * 100;
OnlyPrior_AC_2_PFC_counts   = OnlyPrior_AC_2_PFC_counts / sum(OnlyPrior_AC_2_PFC_counts) * 100;
OnlyPretone_PFC_2_AC_counts = OnlyPretone_PFC_2_AC_counts / sum(OnlyPretone_PFC_2_AC_counts) * 100;
OnlyPretone_AC_2_PFC_counts = OnlyPretone_AC_2_PFC_counts / sum(OnlyPretone_AC_2_PFC_counts) * 100;

% Plot histograms

figure();

hold on;

% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   OnlyPrior_PFC_2_AC_counts,   'FaceColor', [1, 0.6, 0.6],   'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins, -OnlyPrior_AC_2_PFC_counts,   'FaceColor', [0.6, 0.6, 1],   'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
bar(bins,   OnlyPretone_PFC_2_AC_counts, 'FaceColor', [0.8, 0.1, 0.1], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins, -OnlyPretone_AC_2_PFC_counts, 'FaceColor', [0.2, 0.4, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
xlim([-10 10])
%axis([-.025 .025 -65 65]) 
% Adjust legend to include all conditions
leg = legend({'LED PFC→AC', 'LED AC→PFC', 'Pretone PFC→AC', 'Pretone AC→PFC'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size
% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;
