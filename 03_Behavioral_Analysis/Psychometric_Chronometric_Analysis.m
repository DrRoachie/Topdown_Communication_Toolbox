
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
Condition   = 'prior';            % 'prior' or 'pretone'

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
% Customize axes and labels
% min_y = min(mean_RTs(:));  % Minimum y-value from the data
% max_y = max(mean_RTs(:));  % Maximum y-value from the data
% mid_y = (min_y + max_y) / 2;  % Midpoint of the y-range
% 
% % Round min_y down to the nearest increment of 100, and mid_y/max_y up
% min_y_rounded = floor(min_y / 100) * 100;
% mid_y_rounded = ceil(mid_y / 100) * 100;
% max_y_rounded = ceil(max_y / 100) * 100 + 100;  % Add 100 to the rounded max_y
% 
% ylim([min_y_rounded, max_y_rounded]);  % Set y-axis limits based on rounded values
% 
% axis square;                  % Make the plot area square
% set(gcf, 'Position', [100, 100, 600, 600]);  % Set figure window size (in pixels)
% 
% % Set y-axis ticks to rounded minimum, midpoint, and maximum, with a font size of 24
% set(gca, 'YTick', [min_y_rounded, mid_y_rounded, max_y_rounded], ...
%     'YTickLabel', {num2str(min_y_rounded), num2str(mid_y_rounded), num2str(max_y_rounded)}, 'FontSize', 24);

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

