%% Expectation Analysis
% Updated: CR 11/30/25
% Overall, this code compares the observed phase slope index values during a target tone between two expectation conditions (OnlyPrior and OnlyPretone). 

%% Step 1: Book-Keeping 

Frequency_Band    = 'theta';
Behavior          = 'correct';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];                                              % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius'};                                                                                                 % Options: 'MrCassius' or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_24_PSI_Expectation_Analysis_All_Congruency';             % Options: Where the output of the PSI_connectivity_test.m type file is saved                    
sessions          = dir(fullfile(rootdir, '19*'));
savedir           = 'C:\Users\Corey Roach\Desktop\Figure_PSI_Expectation_Analysis';                                          % Options: Where you want to save your figures
savefig           = 'no';                                                                                                   % Options: 'yes' or 'no'

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

%% Step 5: Plot the average PSI value of each channel pair [HEAT MAP]

    % Replace NaNs with zero
    OnlyPrior_averagedGrid(isnan(OnlyPrior_averagedGrid))   = 0;
    OnlyPretone_averagedGrid(isnan(OnlyPretone_averagedGrid)) = 0;
    
    % Determine the global color range (keep these consistent across figs)
    globalMin = -0.02;
    globalMax =  0.02;
    
    % Clamp values to the global color range
    averagedGrid_1_clamped = max(globalMin, min(globalMax, OnlyPrior_averagedGrid));
    averagedGrid_2_clamped = max(globalMin, min(globalMax, OnlyPretone_averagedGrid));



if strcmp(savefig,'yes') == 1
    
    % ==============================
    % Export settings and save path
    % ==============================
    cmap   = jet(256);
    dpi    = 600;
    inchW  = 4; inchH = 4;
    channels = 3:22;                      % adjust if your grid differs
    tickPos  = 1:numel(channels);
    
    % ==============================
    % Heat map 1: OnlyPrior
    % ==============================
    fig1 = figure('Units','inches','Position',[1 1 inchW inchH], 'Color','w');
    imagesc(averagedGrid_1_clamped, [globalMin globalMax]);
    axis square; box on; colormap(cmap); colorbar;
    xlabel('AC Channels'); ylabel('PFC Channels');
    title(sprintf('OnlyPrior-%s-%s', Frequency_Band, Animal), 'Interpreter','none');
    xticks(tickPos); yticks(tickPos);
    xticklabels(arrayfun(@(x) sprintf('ch%02d', x), channels, 'UniformOutput', false));
    yticklabels(arrayfun(@(x) sprintf('ch%02d', x), channels, 'UniformOutput', false));
    set(gca,'Layer','top','LineWidth',0.75,'FontName','Times New Roman','FontSize',10);
    
    base1 = sprintf('%s_%s_OnlyPrior_PSI_Heatmap', Animal, Frequency_Band);
    exportgraphics(fig1, fullfile(savedir, base1 + ".png"), 'Resolution', dpi, 'BackgroundColor','white');
    exportgraphics(fig1, fullfile(savedir, base1 + ".tif"), 'Resolution', dpi, 'BackgroundColor','white');
    exportgraphics(fig1, fullfile(savedir, base1 + ".pdf"), 'ContentType','vector', 'BackgroundColor','white');
    close(fig1);
    
    % ==============================
    % Heat map 2: OnlyPretone
    % ==============================
    fig2 = figure('Units','inches','Position',[1 1 inchW inchH], 'Color','w');
    imagesc(averagedGrid_2_clamped, [globalMin globalMax]);
    axis square; box on; colormap(cmap); colorbar;
    xlabel('AC Channels'); ylabel('PFC Channels');
    title(sprintf('OnlyPretone-%s-%s', Frequency_Band, Animal), 'Interpreter','none');
    xticks(tickPos); yticks(tickPos);
    xticklabels(arrayfun(@(x) sprintf('ch%02d', x), channels, 'UniformOutput', false));
    yticklabels(arrayfun(@(x) sprintf('ch%02d', x), channels, 'UniformOutput', false));
    set(gca,'Layer','top','LineWidth',0.75,'FontName','Times New Roman','FontSize',10);
    
    base2 = sprintf('%s_%s_OnlyPretone_PSI_Heatmap', Animal, Frequency_Band);
    exportgraphics(fig2, fullfile(savedir, base2 + ".png"), 'Resolution', dpi, 'BackgroundColor','white');
    exportgraphics(fig2, fullfile(savedir, base2 + ".tif"), 'Resolution', dpi, 'BackgroundColor','white');
    exportgraphics(fig2, fullfile(savedir, base2 + ".pdf"), 'ContentType','vector', 'BackgroundColor','white');
    close(fig2);
    
    fprintf('Saved PNG, TIF (600 dpi), and vector PDF to:\n%s\n', savedir);

end

%% Step 6: Plot with Mean of Means + 95% CI [SUMMARY OF HEAT MAP]

    abs_OnlyPrior   = abs(OnlyPrior_averagedGrid);
    abs_OnlyPretone = abs(OnlyPretone_averagedGrid);
    
    % Mask out zeros (placeholder values)
    mask_nonzero_OnlyPrior   = OnlyPrior_averagedGrid ~= 0;
    mask_nonzero_OnlyPretone = OnlyPretone_averagedGrid ~= 0;
    
    % Extract directional values
    PFC_2_AC_OnlyPrior   = abs_OnlyPrior(mask_nonzero_OnlyPrior & OnlyPrior_averagedGrid   > 0);
    AC_2_PFC_OnlyPrior   = abs_OnlyPrior(mask_nonzero_OnlyPrior & OnlyPrior_averagedGrid   < 0);
    PFC_2_AC_OnlyPretone = abs_OnlyPretone(mask_nonzero_OnlyPretone & OnlyPretone_averagedGrid > 0);
    AC_2_PFC_OnlyPretone = abs_OnlyPretone(mask_nonzero_OnlyPretone & OnlyPretone_averagedGrid < 0);
    
    group_data = {PFC_2_AC_OnlyPrior, AC_2_PFC_OnlyPrior, PFC_2_AC_OnlyPretone, AC_2_PFC_OnlyPretone};
    means = cellfun(@mean, group_data);
    SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
    CIs   = SEMs * 1.96;  % 95% CI
    
    bar_colors = [222 30 38; 71 133 197; 222 30 38; 71 133 197]/255;
    
    dpi   = 600;
    inchW = 4; inchH = 4;

    
if strcmp(savefig,'yes') == 1

    % ==============================
    % Make figure (4" x 4") and plot
    % ==============================
    fig = figure('Color','w','Units','inches','Position',[1, 1, inchW, inchH]);
    hold on;
    
    % Bars
    for i = 1:4
        bar(i, means(i), 'FaceColor', bar_colors(i,:), 'EdgeColor', 'none');
    end
    
    % 95% CI error bars
    errorbar(1:4, means, CIs, 'k.', 'LineWidth', 1.25);
    
    % Axes formatting
    set(gca, 'XTick', [1.5 3.5], ...
             'XTickLabel', {'OnlyPrior','OnlyPretone'}, ...
             'FontSize', 10, 'FontName', 'Times New Roman', ...
             'LineWidth', 0.75, 'Layer','top');
    ylabel('Mean |PSI Value|', 'FontSize', 10, 'FontName', 'Times New Roman');
    box off
    
    % Y-axis ticks (no scientific notation)
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.YAxis.TickLabelFormat = '%.4f';
    
    % Legend (top-left, no box)
    lg = legend({'PFC→AC', 'AC→PFC'}, 'Location','northwest', ...
                'FontSize', 9, 'FontName','Times New Roman', 'Box','off');
    uistack(lg,'top');
    
    % ==============================
    % Exports (PNG/TIF 600 dpi + vector PDF)
    % ==============================
    base = sprintf('%s_%s_MeanOfMeans_95CI_Summary', Animal, Frequency_Band);
    exportgraphics(fig, fullfile(savedir, base + ".png"), 'Resolution', dpi, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(savedir, base + ".tif"), 'Resolution', dpi, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(savedir, base + ".pdf"), 'ContentType','vector', 'BackgroundColor','white');
    close(fig);
    
    fprintf('Saved PNG, TIF (600 dpi), and vector PDF to:\n%s\n', savedir);

end


%% Step 8: Independent nonparametric 2×2 (SRH) + ALL 6 pairwise tests
% Descriptives reported as mean ± 95% CI (t-interval). Holm-Bonferroni fixed.

%--------------------- Assemble data (independent groups) ------------------
PSI_values = [PFC_2_AC_OnlyPrior(:); AC_2_PFC_OnlyPrior(:); ...
              PFC_2_AC_OnlyPretone(:); AC_2_PFC_OnlyPretone(:)];

Expectation = [repmat({'OnlyPrior'},   numel(PFC_2_AC_OnlyPrior)+numel(AC_2_PFC_OnlyPrior), 1); ...
               repmat({'OnlyPretone'}, numel(PFC_2_AC_OnlyPretone)+numel(AC_2_PFC_OnlyPretone), 1)];

Direction = [repmat({'PFC→AC'}, numel(PFC_2_AC_OnlyPrior), 1); ...
             repmat({'AC→PFC'}, numel(AC_2_PFC_OnlyPrior), 1); ...
             repmat({'PFC→AC'}, numel(PFC_2_AC_OnlyPretone), 1); ...
             repmat({'AC→PFC'}, numel(AC_2_PFC_OnlyPretone), 1)];

Y = PSI_values;  % |PSI| values per your upstream computation

%---------------- SRH omnibus (rank-based 2×2 ANOVA; independent) ----------
[Exp_g, ~] = grp2idx(Expectation);
[Dir_g, ~] = grp2idx(Direction);
N  = numel(Y);
R  = tiedrank(Y);

a = numel(unique(Exp_g));  % 2
b = numel(unique(Dir_g));  % 2

R_A  = accumarray(Exp_g, R, [a 1], @sum, 0);       n_A  = accumarray(Exp_g, 1, [a 1], @sum, 0);
R_B  = accumarray(Dir_g, R, [b 1], @sum, 0);       n_B  = accumarray(Dir_g, 1, [b 1], @sum, 0);
R_AB = accumarray([Exp_g Dir_g], R, [a b], @sum, 0); n_AB = accumarray([Exp_g Dir_g], 1, [a b], @sum, 0);

% Tie correction
[~,~,ties] = unique(Y);
tie_counts = accumarray(ties,1);
C = 1 - sum(tie_counts.^3 - tie_counts) / (N^3 - N);

H_A = (12/(N*(N+1))) * sum((R_A.^2)./n_A) - 3*(N+1);
H_B = (12/(N*(N+1))) * sum((R_B.^2)./n_B) - 3*(N+1);
H_AB= (12/(N*(N+1))) * sum(sum((R_AB.^2)./n_AB)) - H_A - H_B - 3*(N+1);

H_A  = H_A  / C;
H_B  = H_B  / C;
H_AB = H_AB / C;

p_A  = 1 - chi2cdf(H_A,  a-1);
p_B  = 1 - chi2cdf(H_B,  b-1);
p_AB = 1 - chi2cdf(H_AB, (a-1)*(b-1));

SRH_tbl = table(...
    ["Expectation (OnlyPrior vs OnlyPretone)"; "Direction (PFC→AC vs AC→PFC)"; "Expectation×Direction"], ...
    [H_A; H_B; H_AB], ...
    [a-1; b-1; (a-1)*(b-1)], ...
    [p_A; p_B; p_AB], ...
    'VariableNames', {'Effect','H','df','p'});

disp('Scheirer–Ray–Hare nonparametric 2×2 ANOVA (independent groups):');
disp(SRH_tbl);

%---------------- Pairwise independent tests with effect sizes -------------
% 4 cells (Expectation × Direction)
cell_names = { ...
  'OnlyPrior · PFC→AC', ...
  'OnlyPrior · AC→PFC', ...
  'OnlyPretone · PFC→AC', ...
  'OnlyPretone · AC→PFC' };

cell_data = { ...
  PFC_2_AC_OnlyPrior(:), ...
  AC_2_PFC_OnlyPrior(:), ...
  PFC_2_AC_OnlyPretone(:), ...
  AC_2_PFC_OnlyPretone(:) };

% All 6 pairwise comparisons among the 4 cells
pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
m = size(pairs,1);

pair_p  = nan(m,1);
pair_Z  = nan(m,1);
pair_n1 = nan(m,1);
pair_n2 = nan(m,1);
pair_rb = nan(m,1);
pair_HL = nan(m,1);

for ii = 1:m
    i = pairs(ii,1); j = pairs(ii,2);
    x = cell_data{i}; y = cell_data{j};

    [p, ~, stats] = ranksum(x, y);
    pair_p(ii)  = p;
    pair_Z(ii)  = stats.zval;
    n1 = numel(x); n2 = numel(y);
    pair_n1(ii) = n1; pair_n2(ii) = n2;

    % rank-biserial from Mann–Whitney U
    Rx = stats.ranksum;              % sum of ranks for x
    Ux = Rx - n1*(n1+1)/2;
    pair_rb(ii) = 2*Ux/(n1*n2) - 1;

    % Hodges–Lehmann (median of all pairwise differences)
    pair_HL(ii) = median(x - y.', 'all');
end

[pHolm6, ord6] = holm_bonferroni(pair_p);

Pairwise6_tbl = table( ...
    cell_names(pairs(:,1))', cell_names(pairs(:,2))', ...
    pair_n1, pair_n2, pair_Z, pair_p, pHolm6, pair_rb, pair_HL, ...
    'VariableNames', {'Group1','Group2','n1','n2','Z','p_raw','p_Holm','rank_biserial','HL_shift'});

disp('ALL 6 pairwise contrasts (independent), Holm-corrected across 6:');
disp(Pairwise6_tbl(ord6,:));

%----------------- Descriptives as mean ± 95% CI (t-interval) --------------
labels_ED = ["OnlyPrior","OnlyPrior","OnlyPretone","OnlyPretone"];
dirs_ED   = ["PFC→AC","AC→PFC","PFC→AC","AC→PFC"];

n_vec = zeros(4,1); mean_vec = zeros(4,1);
ciL_vec = zeros(4,1); ciU_vec = zeros(4,1); sd_vec = zeros(4,1);

for ii = 1:4
    x = cell_data{ii};
    n_vec(ii) = numel(x);
    [mhat, ci, sdhat] = mean_ci(x);  % t-interval only
    mean_vec(ii) = mhat; ciL_vec(ii) = ci(1); ciU_vec(ii) = ci(2); sd_vec(ii) = sdhat;
end

MeanCI_tbl = table(labels_ED', dirs_ED', n_vec, mean_vec, ciL_vec, ciU_vec, sd_vec, ...
    'VariableNames', {'Expectation','Direction','n','Mean','CI_L','CI_U','SD'});

disp('Cell-wise descriptives (mean ± 95% CI, t-interval):');
disp(MeanCI_tbl);

%---------------------- Return everything in a struct ----------------------
Res = struct('Animal', Animal, 'Band', Frequency_Band, 'Behavior', Behavior, ...
    'SRH_tbl', SRH_tbl, ...
    'Pairwise6_tbl', Pairwise6_tbl, ...
    'MeanCI_tbl', MeanCI_tbl);

assignin('base', ['Res_' Animal], Res);

