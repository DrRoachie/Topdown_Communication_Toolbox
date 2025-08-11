
% Culmination of Sophia_Psychometric/Chronometric series. Calculates and
% plots psychometric and chronometric curves for specified animal and
% condition. 
% 
% Psychometric chunk takes data from Epoc_Cut_2 directory and
% chronometric chunk takes data from DDM excel table. Variables
% 'total_proportion_H' and 'total_RTs' are used in
% Psycho_Chrono_Wilcoxon_Test.m.

%% Define data

Animal      = 'MrM';              % 'MrCassius' or 'MrM'
Epoch       = 'testToneOnset';    % 'testToneOnset' or 'preCueOnset' 
Condition   = 'pretone';            % 'prior' or 'pretone'
datadir     = fullfile('D:\04_Epoc_Cut', Animal, extractBefore(Epoch, 'Onset'));
ddmdir      = 'D:\03_DDM_Decision_Times';
ddm_fName   = '20210511_audiDeci_monkeyBeh_DT_13-Jun-2023';
sessions    = dir(fullfile(datadir, '19*'));


%% generate psychometric curve

total_proportion_H  = [];               % size: prior x SNR x session (number of priors x 9 x number of sessions)

for k = 1:length(sessions)

    % load session data
    RecDate = sessions(k).name;
    fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
    if strcmp(Condition, 'prior')
        load(fullfile(datadir, RecDate, fName), 'SNR', 'prior', 'choice');
    elseif strcmp(Condition, 'pretone')
        load(fullfile(datadir, RecDate, fName), 'SNR', 'pretone', 'choice');
        prior = pretone;
    else
        error('Invalid prior type. Prior type must be prior or pretone.');
    end

    SNR_list        = unique(SNR);      % x-axis of psychometric curve
    prior_list      = unique(prior);    % list of all priors: (H, L, N)
    if strcmp(Condition, 'pretone')
        prior_list(prior_list == 'N') = []; % remove 'N' pretone - since our experiment involves pretone only, neutral pretone is not applicable
    end

    session_prop_H  = []; % size: prior x SNR (number of priors x 9)
    
    for p = 1:length(prior_list) % for each prior H, L, N

        prior_val       = prior_list(p);
        if strcmp(Condition, 'prior')
            prior_indices   = strcmp(prior, prior_val); % get indices of trials with current prior value
        elseif strcmp(Condition, 'pretone')
            prior_indices   = prior == prior_val; % prior_val is a char and requires '==' to compare
        end

        for s = 1:length(SNR_list) % for each unique SNR value

            SNR_val = SNR_list(s);
            SNR_indices = prior_indices & (SNR == SNR_val); % get indices of trials with current prior and SNR value

            count_H     = sum(choice(SNR_indices) == 'H'); % count the number of 'H' choices for current prior and SNR value
            session_prop_H(p,s) = count_H/(length(choice(SNR_indices)));
    
        end
    end
    total_proportion_H = cat(3,total_proportion_H, session_prop_H);
end


mean_proportion_H  = mean(total_proportion_H, 3);
std_proportion_H   = std(total_proportion_H, 0, 3);
n                  = size(total_proportion_H, 3); % number of observations

% Calculate standard error
standard_error     = std_proportion_H / sqrt(n);

% Use Z-score for 95% confidence level (1.96)
CI_lower           = mean_proportion_H - 1.96 * standard_error;
CI_upper           = mean_proportion_H + 1.96 * standard_error;

% Plot psychometric curve with confidence intervals
figure; hold on;

% Define light colors for each prior condition
gray_blue       = [0.00, 0.45, 0.70];  % Light Blue
gray_green      = [0.00, 0.60, 0.50]; % Light Green
gray_yellow     = [0.55, 0.55, 0.00]; % Dandelion Yellow

% Calculate error bars using confidence intervals
CI_error_H = mean_proportion_H - CI_lower; % Upper and lower errors are symmetric

% Plot for High Prior (Light Blue)
h = errorbar(SNR_list, mean_proportion_H(1,:), CI_error_H(1,:), CI_error_H(1,:), ...
    '.-', 'Color', gray_blue, 'LineWidth', 3);
h.CapSize = 0;

% Plot for Low Prior (Light Green)
l = errorbar(SNR_list, mean_proportion_H(2,:), CI_error_H(2,:), CI_error_H(2,:), ...
    '.-', 'Color', gray_green, 'LineWidth', 3);
l.CapSize = 0;

% Check for Neutral Prior condition (Dandelion Yellow)
if strcmp(Condition, 'prior')
    n = errorbar(SNR_list, mean_proportion_H(3,:), CI_error_H(3,:), CI_error_H(3,:), ...
        '.-', 'Color', gray_yellow, 'LineWidth', 3);
    n.CapSize = 0;
end

axis square;                  % Make the plot area square
set(gcf, 'Position', [100, 100, 600, 600]);  % Set figure window size (in pixels)

% Customize axes and labels
ylim([0, 1]);
set(gca, 'XTick', [-2, -1, 0, 1, 2]);  % Set specific x-tick positions

% Set y-axis ticks to 0, 0.5, and 1 with font size 24
%set(gca, 'YTick', [0 0.5 1], 'YTickLabel', {'0', '0.5', '1'}, 'FontSize', 24);


%Show specific x-axis tick marks for 'pretone' condition
if strcmp(Condition, 'pretone')
    set(gca, 'XTick', [-2, -1, 0, 1, 2]);  % Set specific x-tick positions
else
    % % Remove x-axis tick marks and labels for 'prior' condition
    % set(gca, 'XTick', []);  % Remove tick marks
    % set(gca, 'XTickLabel', []);  % Remove labels
end

%Set font to Arial for the entire plot
set(gca, 'FontName', 'Arial');


%% generate chronometric curve

% Adjust defined data
if strcmp(Animal, 'MrCassius')
    Animal2 = 'MrC';
elseif strcmp(Animal, 'MrM')
    Animal2 = Animal;
end

if strcmp(Condition, 'prior')
    Condition2 = 'priorOnly';
elseif strcmp(Condition, 'pretone')
    Condition2 = 'pretone_pLH';
end


% load in DDM table
DDM_table   = readtable(fullfile(ddmdir, ddm_fName));
DDM_table   = DDM_table(strcmp(DDM_table.subject, Animal2) & ...    % 'MrM' or 'MrC'
                        strcmp(DDM_table.ttype, Condition2) & ...   % prior only trials
                        DDM_table.success == 1, ...                 % correct only trials
                        {'subject','session','SNR','prior','ptC','ttype','RT'});
if strcmp(Condition, 'pretone')
    DDM_table = DDM_table(DDM_table.prior == 0, :); % pretone only requires LED to be neutral
end

sessions    = unique(DDM_table.session);
SNR_list    = unique(DDM_table.SNR);
if strcmp(Condition, 'prior')
    prior_list = unique(DDM_table.prior);
else
    prior_list = unique(DDM_table.ptC);
end

total_RTs           = []; % size: prior x SNR x session (number of priors x 9 x number of sessions)

for k = 1:length(sessions)
    
    session_num = sessions(k);
    session_data = DDM_table(DDM_table.session == session_num,:);

    SNR     = session_data.SNR;
    RT      = session_data.RT;
    if strcmp(Condition, 'prior')
        prior = session_data.prior;
    elseif strcmp(Condition, 'pretone')
        prior = session_data.ptC;
    end

    RTs = []; % size: prior x SNR
    
    for p = 1:length(prior_list) % for each prior

        prior_val       = prior_list(p);
        if strcmp(Condition, 'prior')
            prior_indices   = prior == prior_val; % get indices of trials with current prior value (-2, 0, 2)
        elseif strcmp(Condition, 'pretone')
            prior_indices   = strcmp(prior, prior_val); % prior_val is a cell and requires 'strcmp' to compare (HHH, LLL)
        end

        for s = 1:length(SNR_list) % for each unique SNR value

            SNR_val = SNR_list(s);
            SNR_indices = prior_indices & (SNR == SNR_val); % get indices of trials with current prior and SNR value

            mean_RT = nanmean(RT(SNR_indices)); % calculate the mean reaction time for current prior and SNR value
            RTs(p,s) = mean_RT;
    
        end
    end
    total_RTs = cat(3,total_RTs,RTs);
end

% calculate mean reaction time and standard error
mean_RTs    = nanmean(total_RTs,3);
standard_error_RTs     = nanstd(total_RTs, 0, 3) / sqrt(size(total_RTs, 3));

% plot chronometric curve
figure; hold on;

% Define light colors for each prior condition
gray_blue       = [0.00, 0.45, 0.70];  % Light Blue
gray_green      = [0.00, 0.60, 0.50]; % Light Green
gray_yellow     = [0.55, 0.55, 0.00]; % Dandelion Yellow

% Points to connect - first half and second half
first_half       = [1, 2, 3, 4];
second_half      = [6, 7, 8, 9];
first_and_second = [first_half, second_half];

% Plot first and second half
if strcmp(Condition, 'prior')
    % Scatter plot
    l = errorbar(SNR_list(first_and_second), mean_RTs(1, first_and_second), standard_error_RTs(1, first_and_second), ...
        '.', 'Color', gray_green, 'LineWidth', 3);  % Low Prior (Light Green)
    n = errorbar(SNR_list(first_and_second), mean_RTs(2, first_and_second), standard_error_RTs(2, first_and_second), ...
        '.', 'Color', gray_yellow, 'LineWidth', 3);  % Neutral Prior (Dandelion Yellow)
    h = errorbar(SNR_list(first_and_second), mean_RTs(3, first_and_second), standard_error_RTs(3, first_and_second), ...
        '.', 'Color', gray_blue, 'LineWidth', 3);  % High Prior (Light Blue)

    % Remove error bar caps
    l.CapSize = 0;
    n.CapSize = 0;
    h.CapSize = 0;

    % Connecting lines
    plot(SNR_list(first_half), mean_RTs(1, first_half), '-', 'Color', gray_green, 'LineWidth', 3);  % Low Prior
    plot(SNR_list(second_half), mean_RTs(1, second_half), '-', 'Color', gray_green, 'LineWidth', 3);

    plot(SNR_list(first_half), mean_RTs(2, first_half), '-', 'Color', gray_yellow, 'LineWidth', 3);  % Neutral Prior
    plot(SNR_list(second_half), mean_RTs(2, second_half), '-', 'Color', gray_yellow, 'LineWidth', 3);

    plot(SNR_list(first_half), mean_RTs(3, first_half), '-', 'Color', gray_blue, 'LineWidth', 3);  % High Prior
    plot(SNR_list(second_half), mean_RTs(3, second_half), '-', 'Color', gray_blue, 'LineWidth', 3);
    set(gca, 'XTick', [-2, -1, 0, 1, 2]);  % Set specific x-tick positions

    axis square;                  % Make the plot area square
    set(gcf, 'Position', [100, 100, 600, 600]);  % Set figure window size (in pixels)

elseif strcmp(Condition, 'pretone')
    % Scatter plot for pretone condition
    h = errorbar(SNR_list(first_and_second), mean_RTs(1, first_and_second), standard_error_RTs(1, first_and_second), ...
        '.', 'Color', gray_blue, 'LineWidth', 3);  % High Pretone
    l = errorbar(SNR_list(first_and_second), mean_RTs(2, first_and_second), standard_error_RTs(2, first_and_second), ...
        '.', 'Color', gray_green, 'LineWidth', 3);  % Low Pretone

    % Remove error bar caps
    h.CapSize = 0;
    l.CapSize = 0;

    % Connecting lines
    plot(SNR_list(first_half), mean_RTs(1, first_half), '-', 'Color', gray_blue, 'LineWidth', 3);  % High Pretone
    plot(SNR_list(second_half), mean_RTs(1, second_half), '-', 'Color', gray_blue, 'LineWidth', 3);

    plot(SNR_list(first_half), mean_RTs(2, first_half), '-', 'Color', gray_green, 'LineWidth', 3);  % Low Pretone
    plot(SNR_list(second_half), mean_RTs(2, second_half), '-', 'Color', gray_green, 'LineWidth', 3);
    set(gca, 'XTick', [-2, -1, 0, 1, 2]);  % Set specific x-tick positions

end

axis square;                  % Make the plot area square
set(gcf, 'Position', [100, 100, 600, 600]);  % Set figure window size (in pixels)

ylim([400, 1200]);
axis square;                  % Make the plot area square
set(gcf, 'Position', [100, 100, 600, 600]);  % Set figure window size (in pixels)

%Show specific x-axis tick marks for 'pretone' condition
if strcmp(Condition, 'pretone')
    set(gca, 'XTick', [-2, -1, 0, 1, 2]);  % Set specific x-tick positions
else
    % % Remove x-axis tick marks and labels for 'prior' condition
    % set(gca, 'XTick', []);  % Remove tick marks
    % set(gca, 'XTickLabel', []);  % Remove labels
end

% Set font to Arial for the entire plot
set(gca, 'FontName', 'Arial');

%% Binomial Generalized Linear Model with Logit Function to influence of LED on decision behavior and reaction time

if strcmp(Condition, 'prior')

sessions    = dir(fullfile(datadir, '19*'));

% Initialize storage
all_choice = [];
all_LED    = {};
all_SNR    = [];
all_Session = [];

for k = 1:length(sessions)
    % Load session data
    RecDate = sessions(k).name;
    fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
    load(fullfile(datadir, RecDate, fName), 'SNR', 'prior', 'choice');
    
    % Skip bad sessions if needed
    if isempty(choice)
        continue
    end

    % Prepare LED label
    prior_labels = cell(size(prior));
    prior_labels(strcmp(prior, 'H')) = {'Blue'};
    prior_labels(strcmp(prior, 'L')) = {'Green'};
    prior_labels(strcmp(prior, 'N')) = {'Yellow'};

    % Store
    all_choice   = [all_choice; double(choice == 'H')];  % Convert to 1=High, 0=Low
    all_LED      = [all_LED; prior_labels(:)];
    all_SNR      = [all_SNR; SNR(:)];
    all_Session  = [all_Session; repmat(k, length(SNR), 1)];
end

% Create table
T = table(all_choice, categorical(all_LED), all_SNR, categorical(all_Session), ...
          'VariableNames', {'Choice','LED','SNR','Session'});
% Set Yellow as reference level
LED_cat = categorical(all_LED);
LED_cat = reordercats(LED_cat, {'Yellow', 'Blue', 'Green'});

% Create table with new categorical ordering
T = table(all_choice, LED_cat, all_SNR, categorical(all_Session), ...
          'VariableNames', {'Choice','LED','SNR','Session'});

% Split data
high_idx = T.SNR > 0;   % High target trials
low_idx  = T.SNR < 0;   % Low target trials

T_high = T(high_idx, :);
T_low  = T(low_idx, :);

% Fit separate models
glme_high = fitglme(T_high, 'Choice ~ LED * SNR + (1|Session)', ...
                    'Distribution','Binomial','Link','logit');

glme_low = fitglme(T_low, 'Choice ~ LED * SNR + (1|Session)', ...
                   'Distribution','Binomial','Link','logit');

disp(glme_high);
disp(glme_low);

% Reaction Time Analysis: Linear Mixed-Effects Model (LMM)

% Exclude trials where SNR == 0
DDM_table = DDM_table(DDM_table.SNR ~= 0, :);

% Add TargetType variable based on SNR sign
DDM_table.TargetType = repmat({'HighTarget'}, height(DDM_table), 1);
DDM_table.TargetType(DDM_table.SNR < 0) = {'LowTarget'};
DDM_table.TargetType = categorical(DDM_table.TargetType);

% Map LED labels based on prior values (0 = Yellow, 2 = Blue, -2 = Green)
LED_labels = cell(size(DDM_table.prior));
LED_labels(DDM_table.prior == 0) = {'Yellow'};
LED_labels(DDM_table.prior == 2) = {'Blue'};
LED_labels(DDM_table.prior == -2) = {'Green'};
DDM_table.LED = categorical(LED_labels);

% Keep only relevant variables
RT_table = DDM_table(:, {'RT','LED','SNR','TargetType','session'});
RT_table.Session = categorical(RT_table.session);

% Set reference levels
RT_table.LED = reordercats(RT_table.LED, {'Yellow','Blue','Green'});
RT_table.TargetType = reordercats(RT_table.TargetType, {'LowTarget','HighTarget'});

% Fit Linear Mixed Model: RT ~ LED * TargetType + (1|Session)
lme_RT = fitlme(RT_table, 'RT ~ LED * TargetType + (1|Session)');

% Display model summary
disp(lme_RT);

% Optional: ANOVA to test interaction significance
anova_RT = anova(lme_RT);
disp(anova_RT);

% Optional: Plot Marginal Means (Estimated RTs)
[predTable,CI] = predict(lme_RT);

end

%% Pretone Influence on Choice: Binomial GLMM with Logit Link

sessions = dir(fullfile(datadir, '19*'));

% Initialize storage
all_choice  = [];
all_Pretone = {}; 
all_SNR     = [];
all_Session = [];

for k = 1:length(sessions)
    % Load session data
    RecDate = sessions(k).name;
    fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
    
    if strcmp(Condition, 'pretone')
        load(fullfile(datadir, RecDate, fName), 'SNR', 'pretone', 'choice');
        prior = pretone; % rename for consistency
    else
        error('This code is for pretone condition only.');
    end

    % Store everything first (do not exclude yet)
    all_choice  = [all_choice; double(choice == 'H')];
    all_Pretone = [all_Pretone; cellstr(prior(:))];
    all_SNR     = [all_SNR; SNR(:)];
    all_Session = [all_Session; repmat(k, length(SNR), 1)];
end

% >>> Global exclusion of 'N' trials after pooling
valid_idx = ~strcmp(all_Pretone, 'N');
all_choice  = all_choice(valid_idx);
all_Pretone = all_Pretone(valid_idx);
all_SNR     = all_SNR(valid_idx);
all_Session = all_Session(valid_idx);

% Create table
T = table(all_choice, categorical(all_Pretone), all_SNR, categorical(all_Session), ...
          'VariableNames', {'Choice','Pretone','SNR','Session'});

% Set 'L' as reference
T.Pretone = reordercats(T.Pretone, {'L','H'});



%% Split into High and Low target trials
T_highPretone = T(T.SNR > 0, :);
T_lowPretone  = T(T.SNR < 0, :);

% Fit GLMM for High Target trials
glme_highPretone = fitglme(T_highPretone, ...
    'Choice ~ Pretone * SNR + (1|Session)', ...
    'Distribution','Binomial','Link','logit');

disp(glme_highPretone);

% Fit GLMM for Low Target trials
glme_lowPretone = fitglme(T_lowPretone, ...
    'Choice ~ Pretone * SNR + (1|Session)', ...
    'Distribution','Binomial','Link','logit');

disp(glme_lowPretone);


%% Pretone Influence on Reaction Time: Linear Mixed-Effects Model (LMM)

% Load DDM table
DDM_table = readtable(fullfile(ddmdir, ddm_fName));

% Filter for correct subject, pretone condition, and correct trials
DDM_table = DDM_table(strcmp(DDM_table.subject, Animal2) & ...
                      strcmp(DDM_table.ttype, 'pretone_pLH') & ...
                      DDM_table.success == 1, :);

% Exclude SNR == 0 trials
DDM_table = DDM_table(DDM_table.SNR ~= 0, :);

% Add TargetType based on SNR sign
DDM_table.TargetType = repmat({'HighTarget'}, height(DDM_table), 1);
DDM_table.TargetType(DDM_table.SNR < 0) = {'LowTarget'};
DDM_table.TargetType = categorical(DDM_table.TargetType);

% Map Pretone triplet codes ('HHH', 'LLL') to 'H' and 'L'
Pretone_labels = cell(size(DDM_table.ptC));
Pretone_labels(strcmp(DDM_table.ptC, 'HHH')) = {'H'};
Pretone_labels(strcmp(DDM_table.ptC, 'LLL')) = {'L'};
DDM_table.Pretone = categorical(Pretone_labels);

% Remove rows with missing or empty pretone labels (safety check)
DDM_table = DDM_table(~cellfun(@isempty, Pretone_labels), :);

% Build Reaction Time analysis table
RT_table = DDM_table(:, {'RT', 'Pretone', 'TargetType', 'session'});
RT_table.Session = categorical(RT_table.session);

% Set reference levels
RT_table.Pretone = reordercats(RT_table.Pretone, {'L', 'H'});
RT_table.TargetType = reordercats(RT_table.TargetType, {'LowTarget', 'HighTarget'});

%% Fit Linear Mixed-Effects Model
lme_pretoneRT = fitlme(RT_table, 'RT ~ Pretone * TargetType + (1|Session)');

% Display model summary
disp(lme_pretoneRT);

% ANOVA for significance of main effects and interaction
anova_pretoneRT = anova(lme_pretoneRT);
disp(anova_pretoneRT);

%% Optional: Predict Marginal Means (if you’d like plots later)
[predTable, CI] = predict(lme_pretoneRT);

%% SNR-wise Chi-Squared Test 

% Assumes you already have these arrays: all_choice, all_Pretone, all_SNR
% all_choice: vector of binary values (1 = High choice, 0 = Low)
% all_Pretone: cell array of 'H' or 'L' pretone values
% all_SNR: vector of SNR values

% Assumes: all_choice (0/1), all_Pretone ('H'/'L'), all_SNR (numeric)

SNR_list = unique(all_SNR);
p_values = nan(size(SNR_list));
chi2_stats = nan(size(SNR_list));

fprintf('SNR\tChi2\tp-value\n');
fprintf('---------------------------\n');

for i = 1:length(SNR_list)
    snr_val = SNR_list(i);
    
    idx = all_SNR == snr_val;
    choices = all_choice(idx);
    pretones = all_Pretone(idx);
    
    % Build contingency table: rows = Pretone, columns = Choice
    [tbl,~,~,labels] = crosstab(pretones, choices);
    
    % Run chi-squared test of independence
    [chi2stat, p] = chi2cont(tbl);
    
    % Store and display
    p_values(i) = p;
    chi2_stats(i) = chi2stat;
    fprintf('%d\t%.2f\t%.4f\n', snr_val, chi2stat, p);
end

%%
% Sort p-values in ascending order
[p_sorted, sort_idx] = sort(p_values);
n = length(p_sorted);

% Apply Benjamini-Hochberg correction
fdr_corrected = p_sorted .* n ./ (1:n);  % BH formula: p * n / rank
fdr_corrected = min(fdr_corrected, 1);   % Cap at 1

% Enforce monotonicity (non-decreasing correction)
for i = n-1:-1:1
    fdr_corrected(i) = min(fdr_corrected(i), fdr_corrected(i+1));
end

% Reorder to match original SNR order
fdr_corrected_pvals = nan(size(p_values));
fdr_corrected_pvals(sort_idx) = fdr_corrected;

% Display results
fprintf('\nSNR\tChi2\tp-value\tFDR-corrected\n');
fprintf('---------------------------------------------\n');
for i = 1:length(SNR_list)
    fprintf('%.6f\t%.2f\t%.4f\t%.4f\n', ...
        SNR_list(i), ...
        chi2_stats(i), ...
        p_values(i), ...
        fdr_corrected_pvals(i));
end


%%

% Unique SNR values
SNR_list = unique(all_SNR);
p_values_RT = nan(size(SNR_list));
chi2_stats_RT = nan(size(SNR_list));

fprintf('SNR\tChi2\tp-value\n');
fprintf('-----------------------------\n');

for i = 1:length(SNR_list)
    snr_val = SNR_list(i);

    % Get trials at current SNR
    idx = all_SNR == snr_val;
    rt_vals = all_RT(idx);
    pretone_vals = all_Pretone(idx);

    % Skip if too few values
    if sum(idx) < 10
        continue
    end

    % Bin RTs as 'Fast' or 'Slow' based on median split
    median_rt = median(rt_vals);
    rt_bin = repmat("Slow", size(rt_vals));
    rt_bin(rt_vals <= median_rt) = "Fast";

    % Build contingency table (Pretone x RT bin)
    [tbl,~,~,~] = crosstab(pretone_vals, rt_bin);

    % Skip if not enough categories
    if any(size(tbl) < 2)
        continue
    end

    % Chi-squared test
    [chi2stat, p] = chi2cont(tbl);

    % Store and print results
    chi2_stats_RT(i) = chi2stat;
    p_values_RT(i) = p;

    fprintf('%.2f\t%.2f\t%.4f\n', snr_val, chi2stat, p);
end


%%
function [chi2stat, p] = chi2cont(tbl)
% chi2cont - Chi-squared test of independence for contingency table
    expected = sum(tbl, 2) * sum(tbl, 1) / sum(tbl(:));
    chi2stat = sum((tbl - expected).^2 ./ expected, 'all');
    df = (size(tbl,1)-1) * (size(tbl,2)-1);
    p = 1 - chi2cdf(chi2stat, df);
end

%%

function out = ternary(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
