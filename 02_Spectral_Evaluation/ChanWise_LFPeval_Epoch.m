
%% Epoch_Connectivity_Test

% This code tests how connectivity between PFC and AC changes as a function
% of epoch (LED illumination versus testTone presentation). For each
% session, PSI during preCueOnset and testTone onset is calculated. Then
% the session-wise PSI values are visualized and statistically compared.
% Though code for OnlyPretone Trials (neural LED followed by pretones) is
% present, we chose to do this analysis on OnlyPrior trials (informatative
% LED followed by silence). 
 

%% Set parameters for the preCue Condition

Condition         = 'OnlyPrior';                                                  % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];   % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 

%% Calculate PSI Spectra during the preCueOnset epoch; PriorOnly Trials, Collapsed Across Behavior (Correct and Wrong) and SNR (High and Low)
% ahead of the development of this code, we found that there was not a
% strong influence of Behavior/SNR on PSI, so we collapse across them here.
% 

    % Step 1: Get list of modulated channels 

    % Because theta/alpha & beta have the chance to use different channels,
    % you have to make the code Frequency_Band dependent. 

    if strcmp(Frequency_Band, 'theta') || strcmp(Frequency_Band, 'alpha') % theta and alpha have the same shared pairs according to LFP_spec_eval, so they can use the same set/

        session_info = {'MrCassius', '190330';
                        'MrCassius', '190404';
                        'MrCassius', '190413';
                        'MrCassius', '190416';
                        'MrM', '190417';
                        'MrCassius', '190418';
                        'MrCassius', '190419';
                        'MrCassius', '190421';
                        'MrM', '190422';
                        'MrCassius', '190423';
                        'MrM', '190425';
                        'MrM', '190427';
                        'MrM', '190502';
                        'MrM', '190514';
                        'MrCassius', '190515';
                        'MrCassius', '190517';
                        'MrM', '190525';
                        'MrM', '190527';
                        'MrM', '190530';
                        'MrCassius', '190531';
                        'MrM', '190601';
                        'MrCassius', '190603';
                        'MrM', '190604';
                        'MrCassius', '190605';
                        'MrCassius', '190703';
                        'MrM', '190704';
                        'MrM', '190709';
                        'MrCassius', '190711';
                        'MrCassius', '190713';
                        'MrM', '190717';
                        'MrM', '190719';
                        'MrCassius', '190720';
                        'MrM', '190722';
                        'MrCassius', '190723';
                        'MrCassius', '190725';};
     end


    if strcmp(Frequency_Band, 'beta')

        session_info = {'MrCassius', '190404';
                        'MrCassius', '190413';
                        'MrCassius', '190416';
                        'MrCassius', '190418';
                        'MrCassius', '190419';
                        'MrM', '190422';
                        'MrCassius', '190517';
                        'MrM', '190525';
                        'MrM', '190527';
                        'MrM', '190601';
                        'MrCassius', '190703';
                        'MrCassius', '190713';
                        'MrCassius', '190718';
                        'MrM', '190719';
                        'MrCassius', '190723';};

    end

    
%% Initialize counters for each of the 20 possible channels that pass the
% phase shuffles

Channelwise_preCue_sigchans_PFC     = zeros(20,1);
Channelwise_preCue_sigchans_AC      = zeros(20,1);
Channelwise_testTone_sigchans_PFC   = zeros(20,1);
Channelwise_testTone_sigchans_AC    = zeros(20,1);

Channelwise_preCue_LFPpower_PFC     = cell(20, 1);
Channelwise_preCue_LFPpower_AC      = cell(20, 1);
Channelwise_testTone_LFPpower_PFC   = cell(20, 1);
Channelwise_testTone_LFPpower_AC    = cell(20, 1);

%%

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 

        % Get PreCue Significant Channels from Phase Shuffle

        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir = fullfile('D:\04_Epoc_Cut', Animal, 'preCue', RecDate);                                 % epoch data, with the time chunk that you want isolated
        chandir = fullfile('D:\05_Significant_Channels_epoch', Animal, 'preCue', RecDate);     % the sig channels that are output LFP_Spectral_Analysis
        savedir        = 'C:\Users\Corey Roach\Documents\00_DATA\LFP_Epoch_Analysis';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')  % make folders if they don't exist already
           mkdir(savedir);
        end   

        preCue_correct_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'preCueOnset', 'OnlyPrior', 'Correct', Frequency_Band);
        preCue_correct_files = dir(fullfile(chandir, preCue_correct_SigChans_fn));

        if  isempty(preCue_correct_files)

            % generate null output files if there are no sig. channels for Condition_1
            % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_preCue_data.txt', Animal, RecDate, 'OnlyPrior', RecDate, Frequency_Band, 'Correct'));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');

        else % get the labels for the PFC channels and AC channels for correct trials 

        preCue_correct_SigChans = load(fullfile(chandir, preCue_correct_SigChans_fn), 'significant_channels');
        preCue_correct_modifiedCellArray   = regexprep(preCue_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        preCue_correct_PFC_chans  =  preCue_correct_modifiedCellArray(startsWith(preCue_correct_modifiedCellArray, '*PFC'));
        preCue_correct_AC_chans   =  preCue_correct_modifiedCellArray(startsWith(preCue_correct_modifiedCellArray, '*AC'));

        end
       

        % STEP 2: Load preCue data 

        % extract session information
        baseFileName = sessions.name;
        fullFileName = fullfile(sessions.folder, baseFileName);
        fprintf(1, 'Now reading %s/n', fullFileName);
        V = load(fullfile(fullFileName));

        % set parameters required for fieldtrip functions and data selection
        params.choice = V.choice;
        params.err = V.err;
        params.pretone = V.pretone;
        params.pretoneLength = V.pretoneLength;
        params.prior = cell2char(V.prior); 
        params.SNR = V.SNR;

        % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

        V = ConvertStim(V, RecDate);

        % Adds new field in data structure for OnlyPrior Trials that
        % indicates congruence between target and LED

            if strcmp(Condition,'OnlyPrior') == 1

             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.prior));  % Preallocate as a cell array

             % Loop through each element in V.prior and V.target

                 for i = 1:length(V.prior)
                    if strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'H')
                        V.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'L')
                        V.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(V.prior{i}, 'N')
                        V.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'H'))
                        V.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

                 params.congruency = V.congruency;

            end

          % Adds new field in data structure for OnlyPretone Trials that
          % indicates congruence between target and pretone

            if strcmp(Condition,'OnlyPretone') == 1

             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
             V.pretone = cellstr(V.pretone);

             % Loop through each element in V.pretone and V.target

                 for i = 1:length(V.pretone)
                    if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                        V.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                        V.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(V.pretone{i}, 'N')
                        V.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                        V.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

                 params.congruency = V.congruency;

            end

            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = [];                                    % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = [];                                   % Select all trials regarless of congruency
            preCueOnset_data     = selectData(V.data,params,iSelect);     % 

            % Step 2b: Load testToneOnset data 

            % Set directory and get a list of all files in the folder with the desired file name pattern.
            datadir = fullfile('D:\04_Epoc_Cut', Animal, 'testTone', RecDate);
            chandir = fullfile('D:\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);     % the sig channels that are output LFP_Spectral_Analysis

            sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
            addpath(genpath(datadir));
            addpath(genpath(chandir));

            testTone_correct_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'testToneOnset', 'OnlyPrior', 'Correct', Frequency_Band);
            testTone_correct_files = dir(fullfile(chandir, testTone_correct_SigChans_fn));

            if  isempty(testTone_correct_files)

            % generate null output files if there are no sig. channels for Condition_1
            % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_testTone_data.txt', Animal, RecDate, 'OnlyPrior', RecDate, Frequency_Band, 'Correct'));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');

            else % get the labels for the PFC channels and AC channels for correct trials 

            testTone_correct_SigChans          = load(fullfile(chandir, testTone_correct_SigChans_fn), 'significant_channels');
            testTone_correct_modifiedCellArray = regexprep(testTone_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
    
            % Separate elements that start with *PFC and *AC
            testTone_correct_PFC_chans  =  testTone_correct_modifiedCellArray(startsWith(testTone_correct_modifiedCellArray, '*PFC'));
            testTone_correct_AC_chans   =  testTone_correct_modifiedCellArray(startsWith(testTone_correct_modifiedCellArray, '*AC'));

             % extract session information
             baseFileName = sessions.name;
             fullFileName = fullfile(sessions.folder, baseFileName);
             fprintf(1, 'Now reading %s/n', fullFileName);
             V = load(fullfile(fullFileName));

             % set parameters required for fieldtrip functions and data selection
             params.choice = V.choice;
             params.err = V.err;
             params.pretone = V.pretone;
             params.pretoneLength = V.pretoneLength;
             params.prior = cell2char(V.prior); 
             params.SNR = V.SNR;

             % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

             V = ConvertStim(V, RecDate);

             % Adds new field in data structure for OnlyPrior Trials that
             % indicates congruence between target and LED

                if strcmp(Condition,'OnlyPrior') == 1

                     % Initialize V.congruency as a cell array as a new selection factor for select data function  

                     V.congruency = cell(size(V.prior));  % Preallocate as a cell array

                     % Loop through each element in V.prior and V.target

                     for i = 1:length(V.prior)
                        if strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'H')
                            V.congruency{i} = 'congruent';  % If both are 'H'
                        elseif strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'L')
                            V.congruency{i} = 'congruent';  % If both are 'L'
                        elseif strcmp(V.prior{i}, 'N')
                            V.congruency{i} = 'neutral';  % If prior is 'N'
                        elseif (strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'H'))
                            V.congruency{i} = 'incongruent';  % If prior and target are different
                        end
                     end

                     params.congruency = V.congruency;

                end

            % Adds new field in data structure for OnlyPretone Trials that
            % indicates congruence between target and pretone

                if strcmp(Condition,'OnlyPretone') == 1

                 % Initialize V.congruency as a cell array as a new selection factor for select data function  

                 V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
                 V.pretone = cellstr(V.pretone);

                 % Loop through each element in V.pretone and V.target

                     for i = 1:length(V.pretone)
                        if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                            V.congruency{i} = 'congruent';  % If both are 'H'
                        elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                            V.congruency{i} = 'congruent';  % If both are 'L'
                        elseif strcmp(V.pretone{i}, 'N')
                            V.congruency{i} = 'neutral';  % If prior is 'N'
                        elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                            V.congruency{i} = 'incongruent';  % If prior and target are different
                        end
                     end

                     params.congruency = V.congruency;

                 end

            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = [];                                    % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = [];                                   % Select all trials regarless of congruency
            testToneOnset_data   = selectData(V.data,params,iSelect);   

            % Truncate the larger data frame to match smaller values. This
            %  ensures that the testTone epoch and the preCue
            % epoch have the same number of trials.

            % numTrials_preCue     = size(preCueOnset_data.trial, 2);
            % numTrials_testTone   = size(testToneOnset_data.trial, 2);
            % minTrials            = min(numTrials_preCue, numTrials_testTone);
            % 
            % preCueOnset_data.trial      = preCueOnset_data.trial(:, 1:minTrials);
            % preCueOnset_data.time       = preCueOnset_data.time(:, 1:minTrials);
            % preCueOnset_data.sampleinfo = preCueOnset_data.sampleinfo(1:minTrials, :);
            % 
            % testToneOnset_data.trial      = testToneOnset_data.trial(:, 1:minTrials);
            % testToneOnset_data.time       = testToneOnset_data.time(:, 1:minTrials);
            % testToneOnset_data.sampleinfo = testToneOnset_data.sampleinfo(1:minTrials, :);

            % calculate cross for the preCue epoch

            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.taper      = 'dpss';
            cfg.output     = 'fourier';
            cfg.tapsmofrq  = 4; 
            cfg.pad        = 1;
            cfg.foilim     = [0 100];

            freq_preCue         = ft_freqanalysis(cfg, preCueOnset_data);
            fd_preCue           = ft_freqdescriptives(cfg, freq_preCue);
            freq_testTone       = ft_freqanalysis(cfg, testToneOnset_data);
            fd_testTone         = ft_freqdescriptives(cfg, freq_testTone);

           
            % save_file_name = sprintf('%s_%s_%s_%s_PSI_data.mat', Animal, RecDate, 'preCueOnset', Frequency_Band);
            % save(fullfile(savedir, save_file_name), 'PSI_preCueOnset');

                    % This section takes the phase shuffled sig channels, and
                    % extracts the powerspectra from those channels. 

                    % === OnlyPrior PFC ===

                    for i = 1:length(preCue_correct_PFC_chans)
                        chNum = sscanf(preCue_correct_PFC_chans{i}, '*PFC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_preCue_sigchans_PFC(idx) = Channelwise_preCue_sigchans_PFC(idx) + 1; % keeps count
                            Channelwise_preCue_LFPpower_PFC{idx} = [Channelwise_preCue_LFPpower_PFC{idx}; fd_preCue.powspctrm(idx, :)];
                        end
                    end

                    % === OnlyPrior AC ===

                    for i = 1:length(preCue_correct_AC_chans)
                        chNum = sscanf(preCue_correct_AC_chans{i}, '*AC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_preCue_sigchans_AC(idx) = Channelwise_preCue_sigchans_AC(idx) + 1;
                            Channelwise_preCue_LFPpower_AC{idx} = [Channelwise_preCue_LFPpower_AC{idx}; fd_preCue.powspctrm(idx + 20, :)];
                        end
                    end

                    % === OnlyPretone PFC ===

                    for i = 1:length(testTone_correct_PFC_chans)
                        chNum = sscanf(testTone_correct_PFC_chans{i}, '*PFC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_testTone_sigchans_PFC(idx) = Channelwise_testTone_sigchans_PFC(idx) + 1;
                            Channelwise_testTone_LFPpower_PFC{idx} = [Channelwise_testTone_LFPpower_PFC{idx}; fd_testTone.powspctrm(idx, :)];
                        end
                    end

                    % === OnlyPretone AC ===

                    for i = 1:length(testTone_correct_AC_chans)
                        chNum = sscanf(testTone_correct_AC_chans{i}, '*AC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_testTone_sigchans_AC(idx) = Channelwise_testTone_sigchans_AC(idx) + 1;
                            Channelwise_testTone_LFPpower_AC{idx} = [Channelwise_testTone_LFPpower_AC{idx}; fd_testTone.powspctrm(idx + 20, :)];
                        end
                    end


            end
        
    end



%% Plot what percentage of sessions pass phase-shuffle evaluation per channel

propSessions_preCue_PFC     = (Channelwise_preCue_sigchans_PFC/rd)   * 100;
propSessions_preCue_AC      = (Channelwise_preCue_sigchans_AC/rd)    * 100;
propSessions_testTone_PFC   = (Channelwise_testTone_sigchans_PFC/rd) * 100;
propSessions_testTone_AC    = (Channelwise_testTone_sigchans_AC/rd)  * 100;

% Define custom RGB colors
blue   = [0.2 0.4 1];      % OnlyPrior
orange = [1 0.5 0];        % OnlyPretone

if strcmp(Frequency_Band,'theta') == 1
prop_range = [50 100];
end

if strcmp(Frequency_Band,'alpha') == 1
prop_range = [50 100];
end

if strcmp(Frequency_Band,'beta') == 1
prop_range = [1 100];
end

% ----- AC -----
figure();
hold on;
plot(propSessions_preCue_AC,  'Color', blue,   'LineWidth', 2);
plot(propSessions_testTone_AC,'Color', orange, 'LineWidth', 2);
ylim(prop_range);
xlabel('Channel');
ylabel('Proportion of Sessions (%)');
title('Auditory Cortex (AC)');
%legend({'OnlyPrior', 'OnlyPretone'});
hold off;

% ----- PFC -----
figure();
hold on;
plot(propSessions_preCue_PFC,  'Color', blue,   'LineWidth', 2);
plot(propSessions_testTone_PFC,'Color', orange, 'LineWidth', 2);
ylim(prop_range);
xlabel('Channel');
ylabel('Proportion of Sessions (%)');
title('Prefrontal Cortex (PFC)');
%legend({'OnlyPrior', 'OnlyPretone'});
hold off;

% %% Plot power spectra of LED epoch collapsed across trials, pooling across sessions and channels 
% 
% % --- Aggregation ---
% 
% preCue_LFPpower_AC_array = [];
% 
% for ch = 1:20
%     data = Channelwise_preCue_LFPpower_AC{ch};  % n x 101
%     band = data(:, 5:9);                      
%     preCue_LFPpower_AC_array = [preCue_LFPpower_AC_array; band];              % accumulate across channels
% end
% 
% % --- Statistics ---
% avg_power  = mean(preCue_LFPpower_AC_array, 1);                             % 1 x 47
% sem_power  = std(preCue_LFPpower_AC_array, 0, 1) / sqrt(size(preCue_LFPpower_AC_array, 1)); % SEM
% ci95       = 1.96 * sem_power;                               % 95% CI
% upper_ci   = avg_power + ci95;
% lower_ci   = avg_power - ci95;
% 
% freq_range = 5:9;
% 
% % --- Plot ---
% figure;
% sgtitle('CORRECT');
% 
% subplot(2,1,1);
% hold on;
% 
% ax1 = gca;
% set(ax1, 'Color', [1 1 1]);  % White background
% 
% % Red shaded 95% confidence interval
% fill([freq_range, fliplr(freq_range)], ...
%      [upper_ci, fliplr(lower_ci)], ...
%      [1 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');  % Light red fill
% 
% % Red mean line
% plot(freq_range, avg_power, 'r', 'LineWidth', 2);
% 
% xlabel('Frequency (Hz)');
% ylabel('LFP Mean Power');
% title('preCue Power Spectrum (Theta–Gamma)');
% %legend({'95% CI', 'Mean Power'}, 'Location', 'best');
% xlim([5 9]);
% grid on;
% hold off;
% 
% %% Plot power spectra of testTone epoch collapsed across trials, pooling across sessions and channels 
% 
% % --- Aggregation ---
% testTone_LFPpower_AC_array = [];
% 
% for ch = 1:20
%     data = Channelwise_testTone_LFPpower_AC{ch};  % n x 101
%     band = data(:, 5:9);                       % n x 47
%     testTone_LFPpower_AC_array = [testTone_LFPpower_AC_array; band];              % accumulate across channels
% end
% 
% % --- Statistics ---
% avg_power  = mean(testTone_LFPpower_AC_array, 1);                             % 1 x 47
% sem_power  = std(testTone_LFPpower_AC_array, 0, 1) / sqrt(size(testTone_LFPpower_AC_array, 1)); % SEM
% ci95       = 1.96 * sem_power;                               % 95% CI
% upper_ci   = avg_power + ci95;
% lower_ci   = avg_power - ci95;
% 
% freq_range = 5:9;
% 
% % --- Plot ---
% figure;
% sgtitle('CORRECT');
% 
% subplot(2,1,1);
% hold on;
% 
% ax1 = gca;
% set(ax1, 'Color', [1 1 1]);  % White background
% 
% % Red shaded 95% confidence interval
% fill([freq_range, fliplr(freq_range)], ...
%      [upper_ci, fliplr(lower_ci)], ...
%      [1 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');  % Light red fill
% 
% % Red mean line
% plot(freq_range, avg_power, 'r', 'LineWidth', 2);
% 
% xlabel('Frequency (Hz)');
% ylabel('LFP Mean Power');
% title('testTone Power Spectrum (Theta–Gamma)');
% %legend({'95% CI', 'Mean Power'}, 'Location', 'best');
% xlim([5 9]);
% grid on;
% hold off;

%% The mean power spectra across significant sessions per channel after collapsing across trial and frequency bins

            if strcmp(Frequency_Band,'theta') == 1
                freq_range = 5:9;
            end

            if strcmp(Frequency_Band,'alpha') == 1
                freq_range = 10:15;
            end

            if strcmp(Frequency_Band,'beta') == 1
                freq_range = 16:36;
            end

    % average across frequencies and sessions in power spectra from fft 
    
    % Preallocate 20x1 arrays
    Channelwise_preCue_LFPpowerAVG_AC     = zeros(20, 1);
    Channelwise_preCue_LFPpowerAVG_CI_AC  = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = freq_range;
    
    % pick up here, you need to average selectred_data row wise, then
    % calculates CI from n session number...

    % Loop through each channel
    for ch = 1:20
        data = Channelwise_preCue_LFPpower_AC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        %session = mean(selected_data, 2);                   % all values as a vector
        session = selected_data(:);

        % Compute mean and 95% CI
        mu = mean(session);
        se = std(session) / sqrt(length(session));
        ci = 1.96 * se;
    
        % Store results
        Channelwise_preCue_LFPpowerAVG_AC(ch)    = mu;
        Channelwise_preCue_LFPpowerAVG_CI_AC(ch) = ci;
        
    end

    % average across frequencies and sessions in power spectra from fft 
    
    % Preallocate 20x1 arrays
    Channelwise_testTone_LFPpowerAVG_AC     = zeros(20, 1);
    Channelwise_testTone_LFPpowerAVG_CI_AC  = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = freq_range;
    
    % pick up here, you need to average selectred_data row wise, then
    % calculates CI from n session number...

    % Loop through each channel
    for ch = 1:20
        data = Channelwise_testTone_LFPpower_AC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        session = selected_data(:);                   % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(session);
        se = std(session) / sqrt(length(session));
        ci = 1.96 * se;
    
        % Store results
        Channelwise_testTone_LFPpowerAVG_AC(ch)    = mu;
        Channelwise_testTone_LFPpowerAVG_CI_AC(ch) = ci;
        
    end

    % average across frequencies and sessions in power spectra from fft 
    
    % Preallocate 20x1 arrays
    Channelwise_preCue_LFPpowerAVG_PFC     = zeros(20, 1);
    Channelwise_preCue_LFPpowerAVG_CI_PFC  = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = freq_range;
    
    % pick up here, you need to average selectred_data row wise, then
    % calculates CI from n session number...

    % Loop through each channel
    for ch = 1:20
        data = Channelwise_preCue_LFPpower_PFC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        session = selected_data(:);                   % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(session);
        se = std(session) / sqrt(length(session));
        ci = 1.96 * se;
    
        % Store results
        Channelwise_preCue_LFPpowerAVG_PFC(ch)    = mu;
        Channelwise_preCue_LFPpowerAVG_CI_PFC(ch) = ci;
        
    end

    % average across frequencies and sessions in power spectra from fft 
    
    % Preallocate 20x1 arrays
    Channelwise_testTone_LFPpowerAVG_PFC     = zeros(20, 1);
    Channelwise_testTone_LFPpowerAVG_CI_PFC  = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = freq_range;
    
    % pick up here, you need to average selectred_data row wise, then
    % calculates CI from n session number...

    % Loop through each channel
    for ch = 1:20
        data = Channelwise_testTone_LFPpower_PFC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        session = selected_data(:);                  % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(session);
        se = std(session) / sqrt(length(session));
        ci = 1.96 * se;
    
        % Store results
        Channelwise_testTone_LFPpowerAVG_PFC(ch)    = mu;
        Channelwise_testTone_LFPpowerAVG_CI_PFC(ch) = ci;
        
    end


%%
% Define custom RGB colors
blue   = [0.2 0.4 1];      % OnlyPrior
orange = [1 0.5 0];        % OnlyPretone

if strcmp(Frequency_Band,'theta') == 1
freq_range = 5:9;
y_range = [175 4000];
x_range = [50 100];
end

if strcmp(Frequency_Band,'alpha') == 1
freq_range = 10:15;
y_range = [175 4000];
x_range = [50 100];
end

if strcmp(Frequency_Band,'beta') == 1
freq_range = 16:36;
y_range = [50 1000];
x_range = [1 100];
end


% --- Compute AC bounds ---
Channelwise_preCue_LFPpowerAVG_ACupper    = Channelwise_preCue_LFPpowerAVG_AC  + Channelwise_preCue_LFPpowerAVG_CI_AC;
Channelwise_preCue_LFPpowerAVG_AClower    = Channelwise_preCue_LFPpowerAVG_AC  - Channelwise_preCue_LFPpowerAVG_CI_AC;
Channelwise_testTone_LFPpowerAVG_ACupper  = Channelwise_testTone_LFPpowerAVG_AC + Channelwise_testTone_LFPpowerAVG_CI_AC;
Channelwise_testTone_LFPpowerAVG_AClower  = Channelwise_testTone_LFPpowerAVG_AC - Channelwise_testTone_LFPpowerAVG_CI_AC;

% --- Compute PFC bounds ---
Channelwise_preCue_LFPpowerAVG_PFCupper   = Channelwise_preCue_LFPpowerAVG_PFC  + Channelwise_preCue_LFPpowerAVG_CI_PFC;
Channelwise_preCue_LFPpowerAVG_PFClower   = Channelwise_preCue_LFPpowerAVG_PFC  - Channelwise_preCue_LFPpowerAVG_CI_PFC;
Channelwise_testTone_LFPpowerAVG_PFCupper = Channelwise_testTone_LFPpowerAVG_PFC + Channelwise_testTone_LFPpowerAVG_CI_PFC;
Channelwise_testTone_LFPpowerAVG_PFClower = Channelwise_testTone_LFPpowerAVG_PFC - Channelwise_testTone_LFPpowerAVG_CI_PFC;

% --- Plot: AC ---
figure();
hold on;
fill([1:20, fliplr(1:20)], [Channelwise_preCue_LFPpowerAVG_ACupper',   fliplr(Channelwise_preCue_LFPpowerAVG_AClower')],   blue,   'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([1:20, fliplr(1:20)], [Channelwise_testTone_LFPpowerAVG_ACupper', fliplr(Channelwise_testTone_LFPpowerAVG_AClower')], orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Channelwise_preCue_LFPpowerAVG_AC,   'Color', blue,   'LineWidth', 2);
plot(Channelwise_testTone_LFPpowerAVG_AC, 'Color', orange, 'LineWidth', 2);
xlabel('Channel');
ylabel('Mean LFP Power (4–8 Hz)');
title('Auditory Cortex (AC)');

ylim(y_range);
%legend({'', '', 'preCue', 'testTone'});
hold off;

% --- Plot: PFC ---
figure();
hold on;
fill([1:20, fliplr(1:20)], [Channelwise_preCue_LFPpowerAVG_PFCupper',   fliplr(Channelwise_preCue_LFPpowerAVG_PFClower')],   blue,   'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([1:20, fliplr(1:20)], [Channelwise_testTone_LFPpowerAVG_PFCupper', fliplr(Channelwise_testTone_LFPpowerAVG_PFClower')], orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Channelwise_preCue_LFPpowerAVG_PFC,   'Color', blue,   'LineWidth', 2);
plot(Channelwise_testTone_LFPpowerAVG_PFC, 'Color', orange, 'LineWidth', 2);
xlabel('Channel');
ylabel('Mean LFP Power (4–8 Hz)');
title('Prefrontal Cortex (PFC)');
%xlim(freq_range);
ylim(y_range);
%legend({'', '', 'preCue', 'testTone'});
hold off;


% Define custom RGB colors (use same as previous figures)
blue   = [0.2 0.4 1];      % OnlyPrior
orange = [1 0.5 0];        % OnlyPretone

% 
figure();
hold on;
scatter(propSessions_preCue_AC,   Channelwise_preCue_LFPpowerAVG_AC,   60, 'filled', 'MarkerFaceColor', blue);
scatter(propSessions_testTone_AC, Channelwise_testTone_LFPpowerAVG_AC, 60, 'filled', 'MarkerFaceColor', orange);
xlabel('Proportion of Sessions (%)');
ylabel('Mean LFP Power (4–8 Hz)');
title('Auditory Cortex (AC)');
%legend({'preCue', 'testTone'}, 'Location', 'best');
grid on;
xlim(x_range);
ylim(y_range);
hold off;

% ----- PFC: OnlyPrior vs. OnlyPretone -----
figure();
hold on;
scatter(propSessions_preCue_PFC,   Channelwise_preCue_LFPpowerAVG_PFC,   60, 'filled', 'MarkerFaceColor', blue);
scatter(propSessions_testTone_PFC, Channelwise_testTone_LFPpowerAVG_PFC, 60, 'filled', 'MarkerFaceColor', orange);
xlabel('Proportion of Sessions (%)');
ylabel('Mean LFP Power (4–8 Hz)');
title('Prefrontal Cortex (PFC)');
%legend({'preCue', 'testTone'}, 'Location', 'best');
grid on;
xlim(x_range);
%ylim([1 1000]);
hold off;

%% Statistics

% AC - OnlyPrior
[rho_AC_Prior, p_AC_Prior] = corr(propSessions_preCue_AC, Channelwise_preCue_LFPpowerAVG_AC, 'Type', 'Spearman');

% AC - OnlyPretone
[rho_AC_Pretone, p_AC_Pretone] = corr(propSessions_testTone_AC, Channelwise_testTone_LFPpowerAVG_AC, 'Type', 'Spearman');

% PFC - OnlyPrior
[rho_PFC_Prior, p_PFC_Prior] = corr(propSessions_preCue_PFC, Channelwise_preCue_LFPpowerAVG_PFC, 'Type', 'Spearman');

% PFC - OnlyPretone
[rho_PFC_Pretone, p_PFC_Pretone] = corr(propSessions_testTone_PFC, Channelwise_testTone_LFPpowerAVG_PFC, 'Type', 'Spearman');

% %%  Plot Spectrograms 
% 
% % Set one session to visualize (first session)
% Animal    = 'MrCassius';
% RecDate   = '190330';
% Epoch     = 'testToneOnset';
% Condition = 'OnlyPrior';
% Behavior  = 'Correct';
% 
% datadir   = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\02_Preprocessed', Animal, Epoch, RecDate);
% 
% % Get full spectrogram
% fprintf('Getting full spectrogram for %s - %s...\n', Animal, RecDate);
% tfreq = GetSpectrogram(datadir, Animal, RecDate, Epoch, Condition, Behavior);
% 
% % Get baseline and apply z-scoring
% tfreq_baseline = GetSpectrogram(datadir, Animal, RecDate, Epoch, Condition, Behavior);  % Use shorter toi inside function for baseline
% tfreq.zpowspctrm = zscore_spectrogram(tfreq, tfreq_baseline);
% 
% % Plot all PFC channels
% plot_spectrogram(tfreq_OnlyPrior, 'PFC');
% 
% % Plot all AC channels
% plot_spectrogram(tfreq_OnlyPretone, 'AC');