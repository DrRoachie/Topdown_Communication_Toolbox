
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
datadir     = fullfile('D:\04_Epoc_Cut_bipolar', Animal, extractBefore(Epoch, 'Onset'));
ddmdir      = 'D:\03_DDM_Decision_Times';
ddm_fName   = '20210511_audiDeci_monkeyBeh_DT_13-Jun-2023';
sessions    = dir(fullfile(datadir, '19*'));

% ------- Figure export setup -------
savedir = 'C:\Users\Corey Roach\Desktop\Figure_Behavioral_Summary';  

% helper to sanitize strings for filenames
sanitize = @(s) regexprep(s, '[^\w\-]+','_');
baseTag  = sprintf('%s_%s_%s', sanitize(Animal), sanitize(Epoch), sanitize(Condition));

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

% % Plot psychometric curve with confidence intervals
% figure; hold on;

% ------------------- PSYCHOMETRIC: plot & export -------------------
figP = figure('Color','w','Visible','on'); 
axP  = axes(figP); hold(axP,'on');

% Colors
gray_blue   = [0.00, 0.45, 0.70];  % High
gray_green  = [0.00, 0.60, 0.50];  % Low
gray_yellow = [0.55, 0.55, 0.00];  % Neutral

% CI errors
CI_error_H = mean_proportion_H - CI_lower;  % symmetric

% High
h = errorbar(axP, SNR_list, mean_proportion_H(1,:), CI_error_H(1,:), CI_error_H(1,:), ...
    '.-', 'Color', gray_blue, 'LineWidth', 3, 'MarkerSize', 18);  h.CapSize = 0;
% Low
l = errorbar(axP, SNR_list, mean_proportion_H(2,:), CI_error_H(2,:), CI_error_H(2,:), ...
    '.-', 'Color', gray_green, 'LineWidth', 3, 'MarkerSize', 18); l.CapSize = 0;

% Neutral (only for 'prior')
if strcmpi(Condition,'prior') && size(mean_proportion_H,1) >= 3
    nplt = errorbar(axP, SNR_list, mean_proportion_H(3,:), CI_error_H(3,:), CI_error_H(3,:), ...
        '.-', 'Color', gray_yellow, 'LineWidth', 3, 'MarkerSize', 18); nplt.CapSize = 0;
end

% Axes/limits/ticks
axis(axP,'square');
ylim(axP,[0 1]);
xlim(axP,[min(SNR_list) max(SNR_list)]);
set(axP,'XTick',[-2 -1 0 1 2]);

% Typography & no scientific notation
set(axP,'FontName','Arial','FontSize',18,'LineWidth',1.5);
xtickformat(axP,'%.2f');   % ensure no sci-notation
ytickformat(axP,'%.2f');
axP.XAxis.Exponent = 0; axP.YAxis.Exponent = 0;

% Labels (optional – remove/edit to taste)
xlabel(axP,'SNR','FontName','Arial','FontSize',18);
ylabel(axP,'Proportion High','FontName','Arial','FontSize',18);

% Make exactly 4" × 4"
set(figP,'Units','inches'); pos=get(figP,'Position'); pos(3:4)=[4 4]; set(figP,'Position',pos);
set(axP,'LooseInset', max(get(axP,'TightInset'), 0.03*[1 1 1 1])); % compact margins

% Export @ 1200 DPI (PNG + TIF)
basenameP = [baseTag '_psychometric'];
exportgraphics(figP, fullfile(savedir,[basenameP '.png']), 'Resolution',1200, 'BackgroundColor','white');
exportgraphics(figP, fullfile(savedir,[basenameP '.tif']), 'Resolution',1200, 'BackgroundColor','white');



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


% % plot chronometric curve
% figure; hold on;
% 
% Define light colors for each prior condition
% gray_blue       = [0.00, 0.45, 0.70];  % Light Blue
% gray_green      = [0.00, 0.60, 0.50]; % Light Green
% gray_yellow     = [0.55, 0.55, 0.00]; % Dandelion Yellow
% 
% Points to connect - first half and second half
% first_half       = [1, 2, 3, 4];
% second_half      = [6, 7, 8, 9];
% first_and_second = [first_half, second_half];

% ------------------- CHRONOMETRIC: plot & export -------------------
figC = figure('Color','w','Visible','on'); 
axC  = axes(figC); hold(axC,'on');

% Colors
gray_blue   = [0.00, 0.45, 0.70];  % High/HHH
gray_green  = [0.00, 0.60, 0.50];  % Low/LLL
gray_yellow = [0.55, 0.55, 0.00];  % Neutral

% Index blocks you use
first_half       = [1, 2, 3, 4];
second_half      = [6, 7, 8, 9];
first_and_second = [first_half, second_half];

if strcmpi(Condition,'prior')
    l = errorbar(axC, SNR_list(first_and_second), mean_RTs(1, first_and_second), ...
        standard_error_RTs(1, first_and_second), '.', 'Color', gray_green, 'LineWidth', 3, 'MarkerSize', 18); l.CapSize=0;
    n = errorbar(axC, SNR_list(first_and_second), mean_RTs(2, first_and_second), ...
        standard_error_RTs(2, first_and_second), '.', 'Color', gray_yellow, 'LineWidth', 3, 'MarkerSize', 18); n.CapSize=0;
    h = errorbar(axC, SNR_list(first_and_second), mean_RTs(3, first_and_second), ...
        standard_error_RTs(3, first_and_second), '.', 'Color', gray_blue, 'LineWidth', 3, 'MarkerSize', 18); h.CapSize=0;

    plot(axC, SNR_list(first_half),   mean_RTs(1, first_half),   '-', 'Color', gray_green,  'LineWidth',3);
    plot(axC, SNR_list(second_half),  mean_RTs(1, second_half),  '-', 'Color', gray_green,  'LineWidth',3);
    plot(axC, SNR_list(first_half),   mean_RTs(2, first_half),   '-', 'Color', gray_yellow, 'LineWidth',3);
    plot(axC, SNR_list(second_half),  mean_RTs(2, second_half),  '-', 'Color', gray_yellow, 'LineWidth',3);
    plot(axC, SNR_list(first_half),   mean_RTs(3, first_half),   '-', 'Color', gray_blue,   'LineWidth',3);
    plot(axC, SNR_list(second_half),  mean_RTs(3, second_half),  '-', 'Color', gray_blue,   'LineWidth',3);

elseif strcmpi(Condition,'pretone')
    h = errorbar(axC, SNR_list(first_and_second), mean_RTs(1, first_and_second), ...
        standard_error_RTs(1, first_and_second), '.', 'Color', gray_blue, 'LineWidth', 3, 'MarkerSize',18); h.CapSize=0;
    l = errorbar(axC, SNR_list(first_and_second), mean_RTs(2, first_and_second), ...
        standard_error_RTs(2, first_and_second), '.', 'Color', gray_green,'LineWidth', 3, 'MarkerSize',18); l.CapSize=0;

    plot(axC, SNR_list(first_half),   mean_RTs(1, first_half),   '-', 'Color', gray_blue,  'LineWidth',3);
    plot(axC, SNR_list(second_half),  mean_RTs(1, second_half),  '-', 'Color', gray_blue,  'LineWidth',3);
    plot(axC, SNR_list(first_half),   mean_RTs(2, first_half),   '-', 'Color', gray_green, 'LineWidth',3);
    plot(axC, SNR_list(second_half),  mean_RTs(2, second_half),  '-', 'Color', gray_green, 'LineWidth',3);
end

% Axes/limits/ticks
axis(axC,'square');
ylim(axC,[400 1200]);
xlim(axC,[min(SNR_list) max(SNR_list)]);
set(axC,'XTick',[-2 -1 0 1 2]);

% Typography & no scientific notation
set(axC,'FontName','Arial','FontSize',18,'LineWidth',1.5);
xtickformat(axC,'%.2f');   % no sci-notation on SNR
ytickformat(axC,'%.0f');   % integer ms on RT
axC.XAxis.Exponent = 0; axC.YAxis.Exponent = 0;

% Labels (optional – remove/edit to taste)
xlabel(axC,'SNR','FontName','Arial','FontSize',18);
ylabel(axC,'Reaction Time (ms)','FontName','Arial','FontSize',18);

% Make exactly 4" × 4"
set(figC,'Units','inches'); pos=get(figC,'Position'); pos(3:4)=[4 4]; set(figC,'Position',pos);
set(axC,'LooseInset', max(get(axC,'TightInset'), 0.03*[1 1 1 1])); % compact margins

% Export @ 1200 DPI (PNG + TIF)
basenameC = [baseTag '_chronometric'];
exportgraphics(figC, fullfile(savedir,[basenameC '.png']), 'Resolution',1200, 'BackgroundColor','white');
exportgraphics(figC, fullfile(savedir,[basenameC '.tif']), 'Resolution',1200, 'BackgroundColor','white');


% %% LED Prior Statistical Analysis
% 
% % Binomial Generalized Linear Model with Logit Function to influence of LED on decision behavior and reaction time
% 
% if strcmp(Condition, 'prior')
% 
% sessions    = dir(fullfile(datadir, '19*'));
% 
% % Initialize storage
% all_choice = [];
% all_LED    = {};
% all_SNR    = [];
% all_Session = [];
% 
% for k = 1:length(sessions)
%     % Load session data
%     RecDate = sessions(k).name;
%     fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
%     load(fullfile(datadir, RecDate, fName), 'SNR', 'prior', 'choice');
% 
%     % Skip bad sessions if needed
%     if isempty(choice)
%         continue
%     end
% 
%     % Prepare LED label
%     prior_labels = cell(size(prior));
%     prior_labels(strcmp(prior, 'H')) = {'Blue'};
%     prior_labels(strcmp(prior, 'L')) = {'Green'};
%     prior_labels(strcmp(prior, 'N')) = {'Yellow'};
% 
%     % Store
%     all_choice   = [all_choice; double(choice == 'H')];  % Convert to 1=High, 0=Low
%     all_LED      = [all_LED; prior_labels(:)];
%     all_SNR      = [all_SNR; SNR(:)];
%     all_Session  = [all_Session; repmat(k, length(SNR), 1)];
% end
% 
% % Create table
% T = table(all_choice, categorical(all_LED), all_SNR, categorical(all_Session), ...
%           'VariableNames', {'Choice','LED','SNR','Session'});
% % Set Yellow as reference level
% LED_cat = categorical(all_LED);
% LED_cat = reordercats(LED_cat, {'Yellow', 'Blue', 'Green'});
% 
% % Create table with new categorical ordering
% T = table(all_choice, LED_cat, all_SNR, categorical(all_Session), ...
%           'VariableNames', {'Choice','LED','SNR','Session'});
% 
% % Split data
% high_idx = T.SNR > 0;   % High target trials
% low_idx  = T.SNR < 0;   % Low target trials
% 
% T_high = T(high_idx, :);
% T_low  = T(low_idx, :);
% 
% % Fit separate models
% glme_high = fitglme(T_high, 'Choice ~ LED * SNR + (1|Session)', ...
%                     'Distribution','Binomial','Link','logit');
% 
% glme_low = fitglme(T_low, 'Choice ~ LED * SNR + (1|Session)', ...
%                    'Distribution','Binomial','Link','logit');
% 
% disp(glme_high);
% disp(glme_low);
% 
% % Reaction Time Analysis: Linear Mixed-Effects Model (LMM)
% 
% % Exclude trials where SNR == 0
% DDM_table = DDM_table(DDM_table.SNR ~= 0, :);
% 
% % Add TargetType variable based on SNR sign
% DDM_table.TargetType = repmat({'HighTarget'}, height(DDM_table), 1);
% DDM_table.TargetType(DDM_table.SNR < 0) = {'LowTarget'};
% DDM_table.TargetType = categorical(DDM_table.TargetType);
% 
% % Map LED labels based on prior values (0 = Yellow, 2 = Blue, -2 = Green)
% LED_labels = cell(size(DDM_table.prior));
% LED_labels(DDM_table.prior == 0) = {'Yellow'};
% LED_labels(DDM_table.prior == 2) = {'Blue'};
% LED_labels(DDM_table.prior == -2) = {'Green'};
% DDM_table.LED = categorical(LED_labels);
% 
% % Keep only relevant variables
% RT_table = DDM_table(:, {'RT','LED','SNR','TargetType','session'});
% RT_table.Session = categorical(RT_table.session);
% 
% % Set reference levels
% RT_table.LED = reordercats(RT_table.LED, {'Yellow','Blue','Green'});
% RT_table.TargetType = reordercats(RT_table.TargetType, {'LowTarget','HighTarget'});
% 
% % Fit Linear Mixed Model: RT ~ LED * TargetType + (1|Session)
% lme_RT = fitlme(RT_table, 'RT ~ LED * TargetType + (1|Session)');
% 
% % Display model summary
% disp(lme_RT);
% 
% % Optional: ANOVA to test interaction significance
% anova_RT = anova(lme_RT);
% disp(anova_RT);
% 
% 
% 
% end
% 
% %% Pretone Statistical Analysis 
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
% % ===== Pretone GLMM with symmetric piecewise SNR (zero shared) =====
% % Expected in T: Choice (0/1), Pretone_c (H=+0.5, L=-0.5), SNR (double), Session (categorical)
% 
% tol = 1e-6;
% 
% % Construct piecewise SNR (nonnegative) and a zero indicator
% T.SNR_pos = max(T.SNR, 0);
% T.SNR_neg = max(-T.SNR, 0);      % magnitude for negative SNRs
% T.ZeroInd = double(abs(T.SNR) <= tol);  % 1 only for truly-zero SNR trials
% 
% 
% % Fit: separate slopes for +/- SNR; zero is a shared knot (intercept).
% % Include ZeroInd so 0 can have its own offset, and interact Pretone with everything.
% form_piece = 'Choice ~ Pretone_c * (SNR_pos + SNR_neg + ZeroInd) + (1|Session)';
% 
% glme_pw = fitglme(T, form_piece, ...
%     'Distribution','Binomial','Link','logit', ...
%     'DummyVarCoding','reference');
% 
% disp(glme_pw);
% 
% % ---- Simple Pretone (H vs L) contrasts at selected SNRs ----
% % For a given SNR s:
% %  if s > 0:  Δβ = β_P + s*β_P:SNR_pos
% %  if s < 0:  Δβ = β_P + |s|*β_P:SNR_neg
% %  if s = 0:  Δβ = β_P + β_P:ZeroInd
% cn = glme_pw.CoefficientNames;
% b  = glme_pw.Coefficients.Estimate;
% V  = glme_pw.CoefficientCovariance;
% 
% iP     = strcmp(cn,'Pretone_c');
% iP_pos = strcmp(cn,'Pretone_c:SNR_pos');
% iP_neg = strcmp(cn,'Pretone_c:SNR_neg');
% iP_zer = strcmp(cn,'Pretone_c:ZeroInd');
% 
% % robust to missing terms
% if ~any(iP_pos), iP_pos = false(size(b)); end
% if ~any(iP_neg), iP_neg = false(size(b)); end
% if ~any(iP_zer), iP_zer = false(size(b)); end
% 
% zcrit = 1.96;
% s_grid_neg  = [-1.8333, -1.6667, -1.5, -1.25];
% s_grid_pos  = [ 1.25, 1.5, 1.6667, 1.8333];
% 
% fprintf('\n--- Simple Pretone (H vs L) effects with piecewise SNR ---\n');
% 
% % Negative SNR side (Low-target)
% fprintf('Negative SNR (Low-target side):\n');
% for s = s_grid_neg
%     sabs = abs(s);
%     L = zeros(size(b));
%     L(iP)     = 1;
%     L(iP_neg) = sabs;
%     est = L'*b; se = sqrt(L'*V*L);
%     ci = est + [-1 1]*zcrit*se;
%     z = est/se; p = 2*normcdf(-abs(z));
%     fprintf('  SNR=%6.4f:  Δβ=%.3f  SE=%.3f  z=%.2f  p=%.3g  OR=%.2f  (95%% CI %.2f–%.2f)\n', ...
%             s, est, se, z, p, exp(est), exp(ci(1)), exp(ci(2)));
% end
% 
% % Zero SNR (shared)
% L = zeros(size(b)); L(iP)=1; L(iP_zer)=1;
% est = L'*b; se = sqrt(L'*V*L); ci = est + [-1 1]*zcrit*se; z = est/se; p = 2*normcdf(-abs(z));
% fprintf('  SNR= 0.0000:  Δβ=%.3f  SE=%.3f  z=%.2f  p=%.3g  OR=%.2f  (95%% CI %.2f–%.2f)\n', ...
%         est, se, z, p, exp(est), exp(ci(1)), exp(ci(2)));
% 
% % Positive SNR side (High-target)
% fprintf('Positive SNR (High-target side):\n');
% for s = s_grid_pos
%     L = zeros(size(b));
%     L(iP)     = 1;
%     L(iP_pos) = s;
%     est = L'*b; se = sqrt(L'*V*L);
%     ci = est + [-1 1]*zcrit*se;
%     z = est/se; p = 2*normcdf(-abs(z));
%     fprintf('  SNR=%6.4f:  Δβ=%.3f  SE=%.3f  z=%.2f  p=%.3g  OR=%.2f  (95%% CI %.2f–%.2f)\n', ...
%             s, est, se, z, p, exp(est), exp(ci(1)), exp(ci(2)));
% end
% 
% 
% % Reaction Time Analysis for Pretone (LMM)
% % Goal: mirror the LED RT analysis with Pretone (HHH vs LLL) as the predictor.
% % Model: RT ~ Pretone * TargetType + (1|Session)
% 
% % --- Reload a clean DDM subset for safety (subject, ttype=pretone_pLH, correct only) ---
% DDM_pre = readtable(fullfile(ddmdir, ddm_fName));
% DDM_pre = DDM_pre(strcmp(DDM_pre.subject, Animal2) & ...
%                   strcmp(DDM_pre.ttype, 'pretone_pLH') & ...
%                   DDM_pre.success == 1, ...
%                   {'subject','session','SNR','prior','ptC','ttype','RT'});
% 
% % For pretone analyses, LED/prior must be neutral (0)
% DDM_pre = DDM_pre(DDM_pre.prior == 0, :);
% 
% % Exclude SNR == 0 (to separate target sides cleanly, like in the LED analysis)
% DDM_pre = DDM_pre(DDM_pre.SNR ~= 0, :);
% 
% % --- Define TargetType by SNR sign (LowTarget vs HighTarget) ---
% DDM_pre.TargetType = repmat({'HighTarget'}, height(DDM_pre), 1);
% DDM_pre.TargetType(DDM_pre.SNR < 0) = {'LowTarget'};
% DDM_pre.TargetType = categorical(DDM_pre.TargetType);
% 
% % --- Pretone labels from ptC (your data uses 'HHH'/'LLL') ---
% % Keep only HHH/LLL trials; drop anything unexpected
% keep_idx = ismember(DDM_pre.ptC, {'HHH','LLL'});
% DDM_pre = DDM_pre(keep_idx, :);
% 
% % Build Pretone categorical with LLL as reference (so HHH effects are relative to LLL)
% Pretone = categorical(DDM_pre.ptC);
% Pretone = reordercats(Pretone, {'LLL','HHH'});
% 
% % --- Assemble RT table (parallel to the LED version) ---
% RT_table_pre = table( ...
%     DDM_pre.RT, ...
%     Pretone, ...
%     DDM_pre.SNR, ...
%     DDM_pre.TargetType, ...
%     categorical(DDM_pre.session), ...
%     'VariableNames', {'RT','Pretone','SNR','TargetType','Session'});
% 
% % Ensure reference levels match the LED analysis approach
% RT_table_pre.Pretone   = reordercats(RT_table_pre.Pretone, {'LLL','HHH'});
% RT_table_pre.TargetType = reordercats(RT_table_pre.TargetType, {'LowTarget','HighTarget'});
% 
% % --- Fit Linear Mixed Model: RT ~ Pretone * TargetType + (1|Session) ---
% lme_RT_pre = fitlme(RT_table_pre, 'RT ~ Pretone * TargetType + (1|Session)');
% 
% % Display model summary and ANOVA table (interaction tests)
% disp(lme_RT_pre);
% anova_RT_pre = anova(lme_RT_pre);
% disp(anova_RT_pre);
% 
% % Optional: quick interpretation hints (commented)
% % - The Pretone_HHH coefficient is the mean RT difference (HHH vs LLL) on the reference TargetType (LowTarget).
% % - Pretone_HHH:TargetType_HighTarget tests whether the HHH vs LLL RT effect changes for HighTarget trials.
% 
% 
% end 

