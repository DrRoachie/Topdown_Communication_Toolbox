
% For a specified frequency band, loops through all sessions and counts the 
% number of times the p value for coherence for all channel pairs. 
% Heatmap is generated to visualize the spatial
% distribution of the channel pairs that significantly modulated by
% context. A shuffled analysis is carried out that compares the number of
% significant channel pairs across epochs. A second analysis is conducted
% that compares the number of significant interactions between superficial,
% middle, and, deep layers.

%% Define data

Frequency_Band   = 'theta';         % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
animals          = {'MrM'};         % 'MrCassius' and/or 'MrM'
preCue_rootdir   = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\zz_METADATA\2024_10_15_PreCue_Correct\01_Belt';
testTone_rootdir = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\zz_METADATA\2024_09_27_TestTone_Correct\01_Belt';
%savedir          = 'put path to save folder here';

preCue_sessions   = dir(fullfile(preCue_rootdir, '19*'));
testTone_sessions = dir(fullfile(testTone_rootdir, '19*'));

%% Define Data arrays 

% Total Counts
preCue_ChanPair_Array_Coh_sig       = zeros(20,20);   
preCue_ChanPair_Array_Coh_total     = zeros(20, 20); 
testTone_ChanPair_Array_Coh_sig     = zeros(20,20);  
testTone_ChanPair_Array_Coh_total   = zeros(20, 20);  

% Session-wise count in the top-middle-deep layers

preCue_superficial_Channel_Counts   = [];
preCue_middle_Channel_Counts        = [];
preCue_deep_Channel_Counts          = [];

testTone_superficial_Channel_Counts = [];
testTone_middle_Channel_Counts      = [];
testTone_deep_Channel_Counts        = [];

%% Get Session-wise PreCue Counts 

for i = 1:length(preCue_sessions)

    RecDate = preCue_sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        
        if exist(fullfile(preCue_rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
     
        preCue_files = dir(fullfile(preCue_rootdir, RecDate, Animal, 'preCue', ['*' Frequency_Band '*Coh_data.mat']));
           
        else
            continue;
        end

            if ~isempty(preCue_files)
                fprintf('Processing session %s: ', num2str(RecDate));
            end
            
            preCue_p_count      = 0;
            preCue_upper_count  = 0;
            preCue_middle_count = 0;
            preCue_deep_count   = 0;
    
            for k = 1:length(preCue_files) % runs if there are any processed channel pairs in given frequency band
            
                file_name = preCue_files(k).name;
    
                % pull channel pair info from file name
                tokens = regexp(file_name, '(.*?)_', 'tokens');
                if isempty(tokens)
                    continue;
                end
    
                Send_Cort   = tokens{4}{1}; 
                Send_Num    = tokens{5}{1};
                Rec_Cort    = tokens{6}{1};
                Rec_Num     = tokens{7}{1};
    
                Send_Chan   = [Send_Cort '_' Send_Num];
                Rec_Chan    = [Rec_Cort '_' Rec_Num];
    
                % pull all present data files for channel pair
                
                preCue_channel_files = dir(fullfile(preCue_rootdir, RecDate, Animal, 'preCue', [extractBefore(file_name, 'Coh') '*']));
              
                % check if Maris statistical test has been run on the channel pair yet (pair can have Coh_data/Granger_data file without
                % having statistical test results stored in that file since RunMaris_v3.m appends data from statistical test)
               
                if length(preCue_channel_files) ~= 15      % number of files for a channel pair is 15 after Maris runs
                    continue;
                end
                
                % load p value of frequency band
                preCue_datadir = fullfile(preCue_rootdir, RecDate, Animal, 'preCue', file_name);
    
                    % load in p-value
                    load(preCue_datadir,'pvalue_c');
                    pvalue_c = pvalue_c.(Frequency_Band);
        
                    % increase passing p-value count at appropriate index if p-value < 0.05
                    preCue_pfc_idx    = str2double(extractAfter(Send_Num, 'ch')) - 2;    % channels start at ch03
                    preCue_ac_idx     = str2double(extractAfter(Rec_Num, 'ch'))  - 2;    % adjust index to fit all pairs in 20x20 array
    
                    % Increment total appearances for this pair
                    preCue_ChanPair_Array_Coh_total(preCue_pfc_idx, preCue_ac_idx) = preCue_ChanPair_Array_Coh_total(preCue_pfc_idx, preCue_ac_idx) + 1;
    
                    if pvalue_c <= 0.05
                        preCue_ChanPair_Array_Coh_sig(preCue_pfc_idx, preCue_ac_idx) = preCue_ChanPair_Array_Coh_sig(preCue_pfc_idx, preCue_ac_idx) + 1;
                        preCue_p_count = preCue_p_count + 1;
                      
    
                         % Layer-specific counting (based on AC channel involvement)
                            if strcmp(Send_Cort, 'AC') || strcmp(Rec_Cort, 'AC')
                        
                                if strcmp(Rec_Cort, 'AC')
                                    ac_channel = str2double(extractAfter(Rec_Num, 'ch'));
                                elseif strcmp(Send_Cort, 'AC')
                                    ac_channel = str2double(extractAfter(Send_Num, 'ch'));
                                else
                                    ac_channel = NaN;
                                end
                        
                                if ~isnan(ac_channel)
                                    if ac_channel <= 7
                                        preCue_upper_count = preCue_upper_count + 1;
                                    elseif ac_channel >= 8 && ac_channel <= 13
                                        preCue_middle_count = preCue_middle_count + 1;
                                    elseif ac_channel >= 14
                                        preCue_deep_count = preCue_deep_count + 1;
                                    end
                                end
                            end
    
                    end
            end

                    if ~isempty(preCue_files)
                    disp([num2str(preCue_p_count) ' passing p-values found.']);
                    end

                    preCue_superficial_Channel_Counts(end+1)  = preCue_upper_count;
                    preCue_middle_Channel_Counts(end+1)       = preCue_middle_count;
                    preCue_deep_Channel_Counts(end+1)         = preCue_deep_count;

    
    end
end


%% Get Session-wise testTone Counts 

for i = 1:length(testTone_sessions)

    RecDate = testTone_sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        
        if exist(fullfile(testTone_rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
     
        testTone_files = dir(fullfile(testTone_rootdir, RecDate, Animal, 'testTone', ['*' Frequency_Band '*Coh_data.mat']));
           
        else
            continue;
        end

            if ~isempty(testTone_files)
                fprintf('Processing session %s: ', num2str(RecDate));
            end
            
            testTone_p_count      = 0;
            testTone_upper_count  = 0;
            testTone_middle_count = 0;
            testTone_deep_count   = 0;
    
            for k = 1:length(testTone_files) % runs if there are any processed channel pairs in given frequency band
            
                file_name = testTone_files(k).name;
    
                % pull channel pair info from file name
                tokens = regexp(file_name, '(.*?)_', 'tokens');
                if isempty(tokens)
                    continue;
                end
    
                Send_Cort   = tokens{4}{1}; 
                Send_Num    = tokens{5}{1};
                Rec_Cort    = tokens{6}{1};
                Rec_Num     = tokens{7}{1};
    
                Send_Chan   = [Send_Cort '_' Send_Num];
                Rec_Chan    = [Rec_Cort '_' Rec_Num];
    
                % pull all present data files for channel pair
                
                testTone_channel_files = dir(fullfile(testTone_rootdir, RecDate, Animal, 'testTone', [extractBefore(file_name, 'Coh') '*']));
              
                % check if Maris statistical test has been run on the channel pair yet (pair can have Coh_data/Granger_data file without
                % having statistical test results stored in that file since RunMaris_v3.m appends data from statistical test)
               
                if length(testTone_channel_files) ~= 15      % number of files for a channel pair is 15 after Maris runs
                    continue;
                end
                
                % load p value of frequency band
                testTone_datadir = fullfile(testTone_rootdir, RecDate, Animal, 'testTone', file_name);
    
                    % load in p-value
                    load(testTone_datadir,'pvalue_c');
                    pvalue_c = pvalue_c.(Frequency_Band);
        
                    % increase passing p-value count at appropriate index if p-value < 0.05
                    testTone_pfc_idx    = str2double(extractAfter(Send_Num, 'ch')) - 2;    % channels start at ch03
                    testTone_ac_idx     = str2double(extractAfter(Rec_Num, 'ch'))  - 2;    % adjust index to fit all pairs in 20x20 array
    
                    % Increment total appearances for this pair
                    testTone_ChanPair_Array_Coh_total(testTone_pfc_idx, testTone_ac_idx) = testTone_ChanPair_Array_Coh_total(testTone_pfc_idx, testTone_ac_idx) + 1;
    
                    if pvalue_c <= 0.05
                        testTone_ChanPair_Array_Coh_sig(testTone_pfc_idx, testTone_ac_idx) = testTone_ChanPair_Array_Coh_sig(testTone_pfc_idx, testTone_ac_idx) + 1;
                        testTone_p_count = testTone_p_count + 1;  
    
                         % Layer-specific counting (based on AC channel involvement)
                            if strcmp(Send_Cort, 'AC') || strcmp(Rec_Cort, 'AC')
                        
                                if strcmp(Rec_Cort, 'AC')
                                    ac_channel = str2double(extractAfter(Rec_Num, 'ch'));
                                elseif strcmp(Send_Cort, 'AC')
                                    ac_channel = str2double(extractAfter(Send_Num, 'ch'));
                                else
                                    ac_channel = NaN;
                                end
                        
                                if ~isnan(ac_channel)
                                    if ac_channel <= 7
                                        testTone_upper_count = testTone_upper_count + 1;
                                    elseif ac_channel >= 8 && ac_channel <= 13
                                        testTone_middle_count = testTone_middle_count + 1;
                                    elseif ac_channel >= 14
                                        testTone_deep_count = testTone_deep_count + 1;
                                    end
                                end
                            end
    
                    end
            end

                    if ~isempty(testTone_files)
                    disp([num2str(testTone_p_count) ' passing p-values found.']);
                    end

                    testTone_superficial_Channel_Counts(end+1)  = testTone_upper_count;
                    testTone_middle_Channel_Counts(end+1)       = testTone_middle_count;
                    testTone_deep_Channel_Counts(end+1)         = testTone_deep_count;

    
    end
end

%% Wilcoxon signed-rank tests for preCue 

data_preCue = [preCue_superficial_Channel_Counts; 
               preCue_middle_Channel_Counts; 
               preCue_deep_Channel_Counts];  % 3 x N

% Run pairwise Wilcoxon tests
[p1_preCue, ~] = signrank(data_preCue(1,:), data_preCue(2,:));  % Superficial vs Middle
[p2_preCue, ~] = signrank(data_preCue(1,:), data_preCue(3,:));  % Superficial vs Deep
[p3_preCue, ~] = signrank(data_preCue(2,:), data_preCue(3,:));  % Middle vs Deep

% Bonferroni correction
pvals_preCue_raw = [p1_preCue, p2_preCue, p3_preCue];
pvals_preCue_corrected = min(pvals_preCue_raw * 3, 1);  % Cap at 1


%% Wilcoxon signed-rank tests for testTone (row vectors) ---

data_testTone = [testTone_superficial_Channel_Counts; 
                 testTone_middle_Channel_Counts; 
                 testTone_deep_Channel_Counts];  % 3 x N

% Run pairwise Wilcoxon tests
[p1_testTone, ~] = signrank(data_testTone(1,:), data_testTone(2,:));
[p2_testTone, ~] = signrank(data_testTone(1,:), data_testTone(3,:));
[p3_testTone, ~] = signrank(data_testTone(2,:), data_testTone(3,:));

% Bonferroni correction
pvals_testTone_raw = [p1_testTone, p2_testTone, p3_testTone];
pvals_testTone_corrected = min(pvals_testTone_raw * 3, 1);



%% Bar Plot for preCue Epoch 

figure;
hold on;

mean_vals = mean(data_preCue, 2);  % mean across sessions
ci95 = 1.96 * std(data_preCue, 0, 2) ./ sqrt(size(data_preCue, 2));  % 95% CI

bar(mean_vals, 'FaceColor', [0.5 0.7 0.9], 'EdgeColor', 'k');
% errorbar(1:3, mean_vals, ci95, 'k.', 'LineWidth', 1.5);

% Overlay session data
for i = 1:3
    scatter(repmat(i, 1, size(data_preCue,2)), data_preCue(i,:), ...
        30, 'k', 'filled', 'jitter','on', 'jitterAmount', 0.1);
end

xticks(1:3);
xticklabels({'Superficial', 'Middle', 'Deep'});
ylabel('# of Significant Pairs');
title('preCue Epoch');
ylim([0, max(mean_vals + ci95)*1.3]);

% Significance annotations
sigPairs = [1 2; 1 3; 2 3];
ypos = max(mean_vals + ci95)*1.1;
offset = 0.05 * max(mean_vals + ci95);

for i = 1:3
    if pvals_preCue_corrected(i) <= 0.05
        x = sigPairs(i,:);
        y = [ypos ypos];
        plot(x, y, 'k', 'LineWidth', 1.2);

        if pvals_preCue_corrected(i) <= 0.01
            astr = '***';
        else
            astr = '*';
        end
        text(mean(x), ypos + offset, astr, 'FontSize', 14, ...
            'HorizontalAlignment', 'center');
        ypos = ypos + offset * 2;
    end
end
box off

%% Bar Plot for testTone Epoch 

figure;
hold on;
mean_vals = mean(data_testTone, 2);
ci95 = 1.96 * std(data_testTone, 0, 2) ./ sqrt(size(data_testTone, 2));

bar(mean_vals, 'FaceColor', [0.5 0.9 0.7], 'EdgeColor', 'k');
%errorbar(1:3, mean_vals, ci95, 'k.', 'LineWidth', 1.5);

% Overlay session data
for i = 1:3
    scatter(repmat(i, 1, size(data_testTone,2)), data_testTone(i,:), ...
        30, 'k', 'filled', 'jitter','on', 'jitterAmount', 0.1);
end

xticks(1:3);
xticklabels({'Superficial', 'Middle', 'Deep'});
ylabel('# of Significant Pairs');
title('testTone Epoch');
ylim([0, max(mean_vals + ci95)*1.3]);

% Significance annotations
sigPairs = [1 2; 1 3; 2 3];
ypos = max(mean_vals + ci95)*1.1;
offset = 0.05 * max(mean_vals + ci95);

for i = 1:3
    if pvals_testTone_corrected(i) <= 0.05
        x = sigPairs(i,:);
        y = [ypos ypos];
        plot(x, y, 'k', 'LineWidth', 1.2);

        if pvals_testTone_corrected(i) <= 0.01
            astr = '***';
        else
            astr = '*';
        end
        text(mean(x), ypos + offset, astr, 'FontSize', 14, ...
            'HorizontalAlignment', 'center');
        ypos = ypos + offset * 2;
    end
end
box off


%% Plots

figure

% configure grid
s = [20 20]; % [y x]
xrange = [3 22]; % imagesc only needs the endpoints
yrange = [3 22];
dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);

% plot heatmap
imagesc(xrange,yrange, preCue_ChanPair_Array_Coh_sig); hold on
clim([0 8]);   % Set minimum color value to 0 and maximum to 10
colorbar
cb = colorbar;
ylabel(cb, 'Number of Significant Interactions', 'FontSize', 12)
cb.Label.Rotation = 270;

% adjust coloring
hm = mesh(xg,yg,zeros(s+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

title([Frequency_Band, '-', Animal, '-', 'preCue'],...
    'Fontsize', 12);
xlabel('AC Channels', 'FontSize', 12);
ylabel('PFC Channels', 'FontSize', 12);

ax = gca; % Get current axis
ax.FontSize = 12; % Set font size for tick marks

%%

figure 

% configure grid
s = [20 20]; % [y x]
xrange = [3 22]; % imagesc only needs the endpoints
yrange = [3 22];
dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);

% plot heatmap
imagesc(xrange,yrange, testTone_ChanPair_Array_Coh_sig); hold on
clim([0 8]);   % Set minimum color value to 0 and maximum to 10
colorbar
cb = colorbar;
ylabel(cb, 'Number of Significant Interactions', 'FontSize', 12)
cb.Label.Rotation = 270;

% adjust coloring
hm = mesh(xg,yg,zeros(s+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';

title([ Frequency_Band, '-', Animal, '-','testTone'],...
    'Fontsize', 12);
xlabel('AC Channels', 'FontSize', 12);
ylabel('PFC Channels', 'FontSize', 12);

ax = gca; % Get current axis
ax.FontSize = 12; % Set font size for tick marks

%% 

n_iter = 1000;
within_epoch_bootstrap_diffs = zeros(n_iter, 400);  % 100 iterations x 400 channel pairs

for i = 1:n_iter
    % Shuffle within each array
    testTone_shuffled = reshape(testTone_ChanPair_Array_Coh_sig(randperm(400)), 20, 20);
    preCue_shuffled   = reshape(preCue_ChanPair_Array_Coh_sig(randperm(400)), 20, 20);

    % Difference of shuffled values
    diff_array = testTone_shuffled - preCue_shuffled;

    % Store flattened result
    within_epoch_bootstrap_diffs(i, :) = diff_array(:);
end

%%

across_epoch_null_diffs = zeros(n_iter, 400);

% Flatten both arrays
test_flat = testTone_ChanPair_Array_Coh_sig(:);
pre_flat  = preCue_ChanPair_Array_Coh_sig(:);
combined  = [test_flat; pre_flat];  % 800 total values

for i = 1:n_iter
    shuffled = combined(randperm(800));
    
    group1 = reshape(shuffled(1:400), 20, 20);  % simulated testTone
    group2 = reshape(shuffled(401:end), 20, 20);  % simulated preCue
    
    diff_array = group1 - group2;
    across_epoch_null_diffs(i, :) = diff_array(:);
end

%%

figure;
hold on;

% Flatten the matrices to vectors
within_diffs = within_epoch_bootstrap_diffs(:);
null_diffs   = across_epoch_null_diffs(:);

histogram(null_diffs, 20, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability', 'DisplayName', 'Null Distribution');
histogram(within_diffs, 20, 'FaceColor', [0.2 0.2 0.8], 'Normalization', 'probability', 'DisplayName', 'Within-Epoch Bootstrap');

xlabel('Difference (Shuffled TestTone - Shuffled PreCue)');
ylabel('Proportion');
legend;
title('Bootstrap Difference vs Null Distribution');


%%

% Collapse each iteration to a mean difference
mean_within_diffs = mean(within_epoch_bootstrap_diffs, 2);  % 1000 x 1
mean_null_diffs   = mean(across_epoch_null_diffs, 2);       % 1000 x 1

% Mann-Whitney test on bootstrapped means
[p, h, stats] = ranksum(mean_within_diffs, mean_null_diffs);
fprintf('Mann-Whitney test on mean iteration diffs: p = %.4f\n', p);


%% save all figures and statistical outputs 


