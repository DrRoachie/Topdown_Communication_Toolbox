
% Culmination of Sophia_Psychometric/Chronometric series. Calculates and
% plots psychometric and chronometric curves for specified animal and
% condition. 
% 
% Psychometric chunk takes data from Epoc_Cut_2 directory and
% chronometric chunk takes data from DDM excel table. Variables
% 'total_proportion_H' and 'total_RTs' are used in
% Psycho_Chrono_Wilcoxon_Test.m.

%% Define data

Animal      = 'MrCassius';              % 'MrCassius' or 'MrM'
Epoch       = 'testToneOnset';    % 'testToneOnset' or 'preCueOnset' 
Condition   = 'pretone';            % 'prior' or 'pretone'
datadir     = fullfile('D:\04_Epoc_Cut_bipolar', Animal, extractBefore(Epoch, 'Onset'));
ddmdir      = 'D:\03_DDM_Decision_Times';
ddm_fName   = '20210511_audiDeci_monkeyBeh_DT_13-Jun-2023';
sessions    = dir(fullfile(datadir, '19*'));


%% Psychometric Curve Calculation and Visualization

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
% CI_lower           = mean_proportion_H - 1.96 * standard_error;
% CI_upper           = mean_proportion_H + 1.96 * standard_error;

CI_lower           = mean_proportion_H - std_proportion_H;
CI_upper           = mean_proportion_H + std_proportion_H;

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


%% Chronometric Curve Calculation and Visualization 

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


%% LED Prior Analysis

% Binomial Generalized Linear Model with Logit Function to influence of LED on decision behavior and reaction time

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



end

% 
% %% Pretone Analysis Version 1
% 
% % Pretone → Choice (GLMM) with effect-coded Pretone (no neutral baseline)
% 
% if strcmp(Condition, 'pretone')
% 
% sessions = dir(fullfile(datadir, '19*'));
% 
% % Collect trial-wise data
% all_choice  = [];
% all_Pretone = {};
% all_SNR     = [];
% all_Session = [];
% 
% for k = 1:length(sessions)
%     RecDate = sessions(k).name;
%     fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
% 
%     % This block is for Condition == 'pretone'
%     load(fullfile(datadir, RecDate, fName), 'SNR', 'pretone', 'choice');
%     prior = pretone;  % rename for consistency
% 
%     % Append (don't exclude 'N' yet)
%     all_choice  = [all_choice; double(choice == 'H')];
%     all_Pretone = [all_Pretone; cellstr(prior(:))];
%     all_SNR     = [all_SNR; SNR(:)];
%     all_Session = [all_Session; repmat(k, length(SNR), 1)];
% end
% 
% % Exclude neutral/silent pretone 'N'
% valid_idx     = ~strcmp(all_Pretone, 'N');
% all_choice    = all_choice(valid_idx);
% all_Pretone   = all_Pretone(valid_idx);
% all_SNR       = all_SNR(valid_idx);
% all_Session   = all_Session(valid_idx);
% 
% % Build table
% T = table(all_choice, categorical(all_Pretone), all_SNR, categorical(all_Session), ...
%           'VariableNames', {'Choice','Pretone','SNR','Session'});
% 
% % Effect-code Pretone: H = +0.5, L = -0.5  (sum-to-zero; no privileged baseline)
% Pretone_c = -0.5 * ones(height(T),1);
% Pretone_c(T.Pretone == 'H') = +0.5;
% T.Pretone_c = Pretone_c;
% 
% % Split by target side (SNR sign)
% T_highPretone = T(T.SNR > 0, :);   % High target
% T_lowPretone  = T(T.SNR < 0, :);   % Low target
% 
% % Fit GLMMs (binomial, logit)
% glme_highPretone = fitglme(T_highPretone, ...
%     'Choice ~ Pretone_c * SNR + (1|Session)', ...
%     'Distribution','Binomial','Link','logit');
% 
% glme_lowPretone = fitglme(T_lowPretone, ...
%     'Choice ~ Pretone_c * SNR + (1|Session)', ...
%     'Distribution','Binomial','Link','logit');
% 
% disp(glme_highPretone);
% disp(glme_lowPretone);
% 
% % Interpretation tips:
% %  - Intercept: grand-mean log-odds of High choice at SNR=0 (because Pretone_c sums to zero).
% %  - Pretone_c: symmetric H vs L contrast (positive ⇒ H pretone increases odds of High choice vs L).
% %  - Pretone_c:SNR: how that H–L contrast changes with SNR (nonlinearities can be added later via poly(SNR,2)).
% 
% 
% % Pretone → Reaction Time (LMM) using Congruency and AbsSNR
% 
% % Load DDM table
% DDM_table = readtable(fullfile(ddmdir, ddm_fName));
% DDM_table = DDM_table(strcmp(DDM_table.subject, Animal2) & ...
%                       strcmp(DDM_table.ttype, 'pretone_pLH') & ...
%                       DDM_table.success == 1, :);
% 
% % Exclude SNR == 0 trials
% DDM_table = DDM_table(DDM_table.SNR ~= 0, :);
% 
% % Target side (from SNR sign)
% TargetType = repmat({'HighTarget'}, height(DDM_table), 1);
% TargetType(DDM_table.SNR < 0) = {'LowTarget'};
% DDM_table.TargetType = categorical(TargetType);
% 
% % Map ptC triplets to single-letter pretone labels
% Pretone_labels = cell(size(DDM_table.ptC));
% Pretone_labels(strcmp(DDM_table.ptC, 'HHH')) = {'H'};
% Pretone_labels(strcmp(DDM_table.ptC, 'LLL')) = {'L'};
% DDM_table = DDM_table(~cellfun(@isempty, Pretone_labels), :);  % drop any missing
% DDM_table.Pretone = categorical(Pretone_labels(~cellfun(@isempty, Pretone_labels)));
% 
% % Congruency: (H & HighTarget) or (L & LowTarget)
% isHighTarget   = (DDM_table.TargetType == 'HighTarget');
% isHighPretone  = (DDM_table.Pretone     == 'H');
% Congruent      = (isHighTarget & isHighPretone) | (~isHighTarget & ~isHighPretone);
% DDM_table.Congruent = categorical(Congruent, [false true], {'Incongruent','Congruent'});
% 
% % Evidence magnitude
% DDM_table.AbsSNR = abs(DDM_table.SNR);
% 
% % Build RT table
% RT_table = DDM_table(:, {'RT','Congruent','AbsSNR','session','TargetType'});
% RT_table.Session = categorical(RT_table.session);
% 
% % Base model: RT ~ Congruent * AbsSNR + (1|Session)
% lme_pretoneRT = fitlme(RT_table, 'RT ~ Congruent * AbsSNR + (1|Session)');
% 
% % If you want to adjust for any baseline High vs Low target RT difference:
% % lme_pretoneRT = fitlme(RT_table, 'RT ~ Congruent * AbsSNR + TargetType + (1|Session)');
% 
% disp(lme_pretoneRT);
% anova_pretoneRT = anova(lme_pretoneRT);
% disp(anova_pretoneRT);
% 
% % Interpretation tips:
% %  - Congruent (main effect): average RT facilitation (Congruent faster than Incongruent) at AbsSNR=0.
% %  - Congruent:AbsSNR: how the congruency benefit scales with evidence strength.
% %  - If you included TargetType: that term soaks up any baseline High vs Low RT asymmetry.
% 
% end 


% %% Pretone Analysis Version 2
% 
% % Pretone → Choice (GLMM) with effect-coded Pretone (no neutral baseline)
% 
% if strcmp(Condition, 'pretone')
% 
% sessions = dir(fullfile(datadir, '19*'));
% 
% % Collect trial-wise data
% all_choice  = [];
% all_Pretone = {};
% all_SNR     = [];
% all_Session = [];
% 
% for k = 1:length(sessions)
%     RecDate = sessions(k).name;
%     fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
% 
%     % This block is for Condition == 'pretone'
%     load(fullfile(datadir, RecDate, fName), 'SNR', 'pretone', 'choice');
%     prior = pretone;  % rename for consistency
% 
%     % Append (don't exclude 'N' yet)
%     all_choice  = [all_choice; double(choice == 'H')];
%     all_Pretone = [all_Pretone; cellstr(prior(:))];
%     all_SNR     = [all_SNR; SNR(:)];
%     all_Session = [all_Session; repmat(k, length(SNR), 1)];
% end
% 
% % Exclude neutral/silent pretone 'N'
% valid_idx     = ~strcmp(all_Pretone, 'N');
% all_choice    = all_choice(valid_idx);
% all_Pretone   = all_Pretone(valid_idx);
% all_SNR       = all_SNR(valid_idx);
% all_Session   = all_Session(valid_idx);
% 
% % Build table
% T = table(all_choice, categorical(all_Pretone), all_SNR, categorical(all_Session), ...
%           'VariableNames', {'Choice','Pretone','SNR','Session'});
% 
% % Effect-code Pretone: H = +0.5, L = -0.5  (sum-to-zero; no privileged baseline)
% Pretone_c = -0.5 * ones(height(T),1);
% Pretone_c(T.Pretone == 'H') = +0.5;
% T.Pretone_c = Pretone_c;
% 
% % Split by target side (SNR sign)
% T_highPretone = T(T.SNR > 0, :);   % High target
% T_lowPretone  = T(T.SNR < 0, :);   % Low target
% 
% % ===== PRETONE GLMM with 4 SNR categories (rank-safe, no audit) =====
% 
% % Expected in T: Choice (0/1), Pretone_c (H=+0.5, L=-0.5), SNR (double), Session (categorical)
% 
% % ---- Settings -----------------------------------------------------------
% tol = 1e-6;           % tolerance for matching exact fractional SNRs
% zero_policy = 'drop'; % 'drop' | 'high' | 'low'  (where to route SNR==0)
% 
% % ---- Define SNR categories ---------------------------------------------
% High_LowTone_vals  = [-11/6, -5/3];     % strong toward Low target
% Low_LowTone_vals   = [ -3/2, -5/4];     % weaker toward Low target
% Low_HighTone_vals  = [  5/4,  3/2];     % weaker toward High target
% High_HighTone_vals = [  5/3, 11/6];     % strong toward High target
% 
% is_HLT = ismembertol(T.SNR, High_LowTone_vals,  tol);
% is_LLT = ismembertol(T.SNR, Low_LowTone_vals,   tol);
% is_LHT = ismembertol(T.SNR, Low_HighTone_vals,  tol);
% is_HHT = ismembertol(T.SNR, High_HighTone_vals, tol);
% 
% % Build SNR_cat strings
% SNR_cat = strings(height(T),1);
% SNR_cat(is_HLT) = "High_LowTone";
% SNR_cat(is_LLT) = "Low_LowTone";
% SNR_cat(is_LHT) = "Low_HighTone";
% SNR_cat(is_HHT) = "High_HighTone";
% 
% % Keep only rows that match one of the four bins
% keep = SNR_cat ~= "";
% T = T(keep, :);
% SNR_cat = SNR_cat(keep);
% T.SNR_cat = categorical(SNR_cat);  % only present categories are created
% 
% % ---- Split by target side and route zeros per policy -------------------
% isZero       = abs(T.SNR) <= tol;
% isHighTarget = T.SNR >  tol;
% isLowTarget  = T.SNR < -tol;
% 
% switch lower(zero_policy)
%     case 'drop'
%         Th = T(isHighTarget, :);
%         Tl = T(isLowTarget,  :);
%     case 'high'
%         Th = T(isHighTarget | isZero, :);
%         Tl = T(isLowTarget,           :);
%     case 'low'
%         Th = T(isHighTarget,          :);
%         Tl = T(isLowTarget  | isZero, :);
%     otherwise
%         error('zero_policy must be ''drop'', ''high'', or ''low''.');
% end
% 
% % ---- Remove unused levels and set reference levels per split -----------
% % High-target split should contain: Low_HighTone (ref), High_HighTone
% if ~isempty(Th)
%     Th.SNR_cat = removecats(Th.SNR_cat);
%     ch = categories(Th.SNR_cat);
%     % Put Low_HighTone as reference if present
%     if ismember('Low_HighTone', ch)
%         Th.SNR_cat = reordercats(Th.SNR_cat, ['Low_HighTone'; setdiff(ch,'Low_HighTone','stable')]);
%     end
% end
% 
% % Low-target split should contain: Low_LowTone (ref), High_LowTone
% if ~isempty(Tl)
%     Tl.SNR_cat = removecats(Tl.SNR_cat);
%     cl = categories(Tl.SNR_cat);
%     % Put Low_LowTone as reference if present
%     if ismember('Low_LowTone', cl)
%         Tl.SNR_cat = reordercats(Tl.SNR_cat, ['Low_LowTone'; setdiff(cl,'Low_LowTone','stable')]);
%     end
% end
% 
% % ---- Choose rank-safe formulas (fallback if a factor collapses) --------
% form_high = 'Choice ~ Pretone_c * SNR_cat + (1|Session)';
% form_low  = 'Choice ~ Pretone_c * SNR_cat + (1|Session)';
% 
% if isempty(Th) || numel(categories(Th.SNR_cat)) < 2, form_high = 'Choice ~ Pretone_c + (1|Session)'; end
% if ~isempty(Th) && numel(unique(Th.Pretone_c)) < 2, form_high = 'Choice ~ 1 + (1|Session)'; end
% 
% if isempty(Tl) || numel(categories(Tl.SNR_cat)) < 2, form_low  = 'Choice ~ Pretone_c + (1|Session)'; end
% if ~isempty(Tl) && numel(unique(Tl.Pretone_c)) < 2, form_low  = 'Choice ~ 1 + (1|Session)'; end
% 
% % ---- Fit GLMMs ---------------------------------------------------------
% if ~isempty(Th)
%     glme_high4 = fitglme(Th, form_high, 'Distribution','Binomial','Link','logit', ...
%                          'DummyVarCoding','reference');
%     disp(glme_high4);
% end
% if ~isempty(Tl)
%     glme_low4  = fitglme(Tl, form_low,  'Distribution','Binomial','Link','logit', ...
%                          'DummyVarCoding','reference');
%     disp(glme_low4);
% end
% 
% % ---- Print simple Pretone effects within each present bin --------------
% % High-target
% if exist('glme_high4','var')
%     cn = glme_high4.CoefficientNames;
%     b  = glme_high4.Coefficients.Estimate;
%     V  = glme_high4.CoefficientCovariance;
% 
%     % Reference bin name (first category of Th.SNR_cat)
%     refcats_high = categories(Th.SNR_cat);
%     if ~isempty(refcats_high)
%         ref_high = char(refcats_high{1});   % convert to char for %s
%     else
%         ref_high = 'REF';
%     end
% 
%     % Pretone effect within reference bin
%     if any(strcmp(cn,'Pretone_c'))
%         iPret = strcmp(cn,'Pretone_c');
%         L = zeros(size(b)); L(iPret) = 1;
%         est = L'*b; se = sqrt(L'*V*L); z = est/se; p = 2*normcdf(-abs(z));
%         fprintf('High-target: Pretone effect in %s: est=%.3f, SE=%.3f, z=%.2f, p=%.3g, OR=%.2f\n', ...
%                 ref_high, est, se, z, p, exp(est));
%     end
% 
%     % Pretone effect within High_HighTone (if present)
%     if any(strcmp(cn,'Pretone_c:SNR_cat_High_HighTone'))
%         iPret = strcmp(cn,'Pretone_c');
%         iInt  = strcmp(cn,'Pretone_c:SNR_cat_High_HighTone');
%         L = zeros(size(b)); L(iPret)=1; L(iInt)=1;
%         est = L'*b; se = sqrt(L'*V*L); z = est/se; p = 2*normcdf(-abs(z));
%         fprintf('High-target: Pretone effect in High_HighTone: est=%.3f, SE=%.3f, z=%.2f, p=%.3g, OR=%.2f\n', ...
%                 est, se, z, p, exp(est));
%     end
% end
% 
% % Low-target
% if exist('glme_low4','var')
%     cn = glme_low4.CoefficientNames;
%     b  = glme_low4.Coefficients.Estimate;
%     V  = glme_low4.CoefficientCovariance;
% 
%     % Reference bin name (first category of Tl.SNR_cat)
%     refcats_low = categories(Tl.SNR_cat);
%     if ~isempty(refcats_low)
%         ref_low = char(refcats_low{1});
%     else
%         ref_low = 'REF';
%     end
% 
%     % Pretone effect within reference bin
%     if any(strcmp(cn,'Pretone_c'))
%         iPret = strcmp(cn,'Pretone_c');
%         L = zeros(size(b)); L(iPret) = 1;
%         est = L'*b; se = sqrt(L'*V*L); z = est/se; p = 2*normcdf(-abs(z));
%         fprintf('Low-target:  Pretone effect in %s: est=%.3f, SE=%.3f, z=%.2f, p=%.3g, OR=%.2f\n', ...
%                 ref_low, est, se, z, p, exp(est));
%     end
% 
%     % Pretone effect within High_LowTone (if present)
%     if any(strcmp(cn,'Pretone_c:SNR_cat_High_LowTone'))
%         iPret = strcmp(cn,'Pretone_c');
%         iInt  = strcmp(cn,'Pretone_c:SNR_cat_High_LowTone');
%         L = zeros(size(b)); L(iPret)=1; L(iInt)=1;
%         est = L'*b; se = sqrt(L'*V*L); z = est/se; p = 2*normcdf(-abs(z));
%         fprintf('Low-target:  Pretone effect in High_LowTone: est=%.3f, SE=%.3f, z=%.2f, p=%.3g, OR=%.2f\n', ...
%                 est, se, z, p, exp(est));
%     end
% end
% 
% end 

%% Pretone Analysis Version 3

% Pretone → Choice (GLMM) with effect-coded Pretone (no neutral baseline)

if strcmp(Condition, 'pretone')

sessions = dir(fullfile(datadir, '19*'));

% Collect trial-wise data
all_choice  = [];
all_Pretone = {};
all_SNR     = [];
all_Session = [];

for k = 1:length(sessions)
    RecDate = sessions(k).name;
    fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');

    % This block is for Condition == 'pretone'
    load(fullfile(datadir, RecDate, fName), 'SNR', 'pretone', 'choice');
    prior = pretone;  % rename for consistency

    % Append (don't exclude 'N' yet)
    all_choice  = [all_choice; double(choice == 'H')];
    all_Pretone = [all_Pretone; cellstr(prior(:))];
    all_SNR     = [all_SNR; SNR(:)];
    all_Session = [all_Session; repmat(k, length(SNR), 1)];
end

% Exclude neutral/silent pretone 'N'
valid_idx     = ~strcmp(all_Pretone, 'N');
all_choice    = all_choice(valid_idx);
all_Pretone   = all_Pretone(valid_idx);
all_SNR       = all_SNR(valid_idx);
all_Session   = all_Session(valid_idx);

% Build table
T = table(all_choice, categorical(all_Pretone), all_SNR, categorical(all_Session), ...
          'VariableNames', {'Choice','Pretone','SNR','Session'});

% Effect-code Pretone: H = +0.5, L = -0.5  (sum-to-zero; no privileged baseline)
Pretone_c = -0.5 * ones(height(T),1);
Pretone_c(T.Pretone == 'H') = +0.5;
T.Pretone_c = Pretone_c;

% Split by target side (SNR sign)
T_highPretone = T(T.SNR > 0, :);   % High target
T_lowPretone  = T(T.SNR < 0, :);   % Low target

% ===== Pretone GLMM with symmetric piecewise SNR (zero shared) =====
% Expected in T: Choice (0/1), Pretone_c (H=+0.5, L=-0.5), SNR (double), Session (categorical)

tol = 1e-6;

% Construct piecewise SNR (nonnegative) and a zero indicator
T.SNR_pos = max(T.SNR, 0);
T.SNR_neg = max(-T.SNR, 0);      % magnitude for negative SNRs
T.ZeroInd = double(abs(T.SNR) <= tol);  % 1 only for truly-zero SNR trials

% Rank-sane guard: drop rows not on your SNR grid if you want (optional)
% (Skip this if you want to keep everything)
% keep = ismembertol(T.SNR, [-11/6,-5/3,-3/2,-5/4,0,5/4,3/2,5/3,11/6], tol);
% T = T(keep,:);

% Fit: separate slopes for +/- SNR; zero is a shared knot (intercept).
% Include ZeroInd so 0 can have its own offset, and interact Pretone with everything.
form_piece = 'Choice ~ Pretone_c * (SNR_pos + SNR_neg + ZeroInd) + (1|Session)';

glme_pw = fitglme(T, form_piece, ...
    'Distribution','Binomial','Link','logit', ...
    'DummyVarCoding','reference');

disp(glme_pw);

% ---- Simple Pretone (H vs L) contrasts at selected SNRs ----
% For a given SNR s:
%  if s > 0:  Δβ = β_P + s*β_P:SNR_pos
%  if s < 0:  Δβ = β_P + |s|*β_P:SNR_neg
%  if s = 0:  Δβ = β_P + β_P:ZeroInd
cn = glme_pw.CoefficientNames;
b  = glme_pw.Coefficients.Estimate;
V  = glme_pw.CoefficientCovariance;

iP     = strcmp(cn,'Pretone_c');
iP_pos = strcmp(cn,'Pretone_c:SNR_pos');
iP_neg = strcmp(cn,'Pretone_c:SNR_neg');
iP_zer = strcmp(cn,'Pretone_c:ZeroInd');

% robust to missing terms
if ~any(iP_pos), iP_pos = false(size(b)); end
if ~any(iP_neg), iP_neg = false(size(b)); end
if ~any(iP_zer), iP_zer = false(size(b)); end

zcrit = 1.96;
s_grid_neg  = [-1.8333, -1.6667, -1.5, -1.25];
s_grid_pos  = [ 1.25, 1.5, 1.6667, 1.8333];

fprintf('\n--- Simple Pretone (H vs L) effects with piecewise SNR ---\n');

% Negative SNR side (Low-target)
fprintf('Negative SNR (Low-target side):\n');
for s = s_grid_neg
    sabs = abs(s);
    L = zeros(size(b));
    L(iP)     = 1;
    L(iP_neg) = sabs;
    est = L'*b; se = sqrt(L'*V*L);
    ci = est + [-1 1]*zcrit*se;
    z = est/se; p = 2*normcdf(-abs(z));
    fprintf('  SNR=%6.4f:  Δβ=%.3f  SE=%.3f  z=%.2f  p=%.3g  OR=%.2f  (95%% CI %.2f–%.2f)\n', ...
            s, est, se, z, p, exp(est), exp(ci(1)), exp(ci(2)));
end

% Zero SNR (shared)
L = zeros(size(b)); L(iP)=1; L(iP_zer)=1;
est = L'*b; se = sqrt(L'*V*L); ci = est + [-1 1]*zcrit*se; z = est/se; p = 2*normcdf(-abs(z));
fprintf('  SNR= 0.0000:  Δβ=%.3f  SE=%.3f  z=%.2f  p=%.3g  OR=%.2f  (95%% CI %.2f–%.2f)\n', ...
        est, se, z, p, exp(est), exp(ci(1)), exp(ci(2)));

% Positive SNR side (High-target)
fprintf('Positive SNR (High-target side):\n');
for s = s_grid_pos
    L = zeros(size(b));
    L(iP)     = 1;
    L(iP_pos) = s;
    est = L'*b; se = sqrt(L'*V*L);
    ci = est + [-1 1]*zcrit*se;
    z = est/se; p = 2*normcdf(-abs(z));
    fprintf('  SNR=%6.4f:  Δβ=%.3f  SE=%.3f  z=%.2f  p=%.3g  OR=%.2f  (95%% CI %.2f–%.2f)\n', ...
            s, est, se, z, p, exp(est), exp(ci(1)), exp(ci(2)));
end

end 

