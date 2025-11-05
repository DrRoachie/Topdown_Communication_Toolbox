%% Congruency_Connectivity_Test_Epochs_OnlyPrior

% This code summarizes the LFP data used to carry out the connectivity study.
% Now modified to:
% - Process Epochs ('LED' and 'testTone')
% - Use OnlyPrior condition only
% - Save results to a new folder "SpecBalanced_Eval_Epoch"

%% Set parameters for the preCue Condition

Frequency_Band = 'theta';
SNR = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];

%% Session info

session_info = { 'MrCassius', '190330';
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
                 'MrCassius', '190725' };

Epochs_to_process = {'LED', 'testTone'};

%% Process Epochs

    % Session-wise processing
    for rd = 1:length(session_info(:,1))

        % Initialize storage for this session
        theta_power_struct = struct();
        labels_struct = struct();

        % Session info
        Animal  = session_info{rd,1};
        RecDate = session_info{rd,2};

        for iEpoch = 1:length(Epochs_to_process)
            EpochLabel = Epochs_to_process{iEpoch};
            if strcmp(EpochLabel, 'LED')
                EpochFolderName = 'preCue';
                FileEpochName   = 'preCueOnset';
            else
                EpochFolderName = 'testTone';
                FileEpochName = 'testToneOnset';
            end
       
        % Set directory and get a list of all files in the folder with the desired file name pattern
        datadir = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\02_Preprocessed', Animal, FileEpochName, RecDate);
        chandir = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_epoch', Animal, EpochFolderName, RecDate);
        savedir = 'C:\Users\auditory research la\OneDrive\Desktop\SpecBalanced_Eval_Epoch';
        % savedir = fullfile(savedir, RecDate, Animal, EpochLabel);
        savedir = fullfile(savedir, RecDate, Animal, 'OnlyPrior');


        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end

        % File name
        OnlyPrior_correct_SigChans_fn = sprintf('SigChannels_%s_%s_OnlyPrior_Correct_%s.mat', ...
            RecDate, FileEpochName, Frequency_Band);
        
        % Check file
        OnlyPrior_correct_files = dir(fullfile(chandir, OnlyPrior_correct_SigChans_fn));

        if isempty(OnlyPrior_correct_files)
            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', RecDate, Animal, EpochFolderName, 'OnlyPrior', 'Correct', Frequency_Band));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');
            fprintf('Skipping %s %s %s — no sig channels\n', RecDate, Animal, EpochLabel);
        continue

        else
            OnlyPrior_correct_SigChans = load(fullfile(chandir, OnlyPrior_correct_SigChans_fn), 'significant_channels');
            OnlyPrior_correct_modifiedCellArray = regexprep(OnlyPrior_correct_SigChans.significant_channels, '^(D1_|D2_|D3_|D4_)', '*');
            OnlyPrior_correct_PFC_chans = OnlyPrior_correct_modifiedCellArray(startsWith(OnlyPrior_correct_modifiedCellArray, '*PFC'));
            OnlyPrior_correct_AC_chans  = OnlyPrior_correct_modifiedCellArray(startsWith(OnlyPrior_correct_modifiedCellArray, '*AC'));
        end

        % Load data

        sessions = dir(fullfile(datadir,'*.mat'));
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % Assume one file per folder, or pick the first .mat file
        baseFileName = sessions(1).name;
        fullFileName = fullfile(sessions(1).folder, baseFileName);
        fullFileName = fullfile(sessions.folder, baseFileName);
        fprintf(1, 'Now reading %s\n', fullFileName);
        OnlyPrior_data = load(fullfile(fullFileName));

        % Set parameters required for fieldtrip functions and data selection
        OnlyPrior_params.choice = OnlyPrior_data.choice;
        OnlyPrior_params.err = OnlyPrior_data.err;
        OnlyPrior_params.pretone = OnlyPrior_data.pretone;
        OnlyPrior_params.pretoneLength = OnlyPrior_data.pretoneLength;
        OnlyPrior_params.prior = cell2char(OnlyPrior_data.prior);
        OnlyPrior_params.SNR = OnlyPrior_data.SNR;

        % Convert stim labels
        OnlyPrior_data = ConvertStim(OnlyPrior_data, RecDate);

        % Add congruency field
        OnlyPrior_data.congruency = cell(size(OnlyPrior_data.prior));
        OnlyPrior_data.pretone = cellstr(OnlyPrior_data.pretone);
        for i = 1:length(OnlyPrior_data.prior)
            if strcmp(OnlyPrior_data.prior{i}, 'H') && strcmp(OnlyPrior_data.target{i}, 'H')
                OnlyPrior_data.congruency{i} = 'congruent';
            elseif strcmp(OnlyPrior_data.prior{i}, 'L') && strcmp(OnlyPrior_data.target{i}, 'L')
                OnlyPrior_data.congruency{i} = 'congruent';
            elseif strcmp(OnlyPrior_data.prior{i}, 'N')
                OnlyPrior_data.congruency{i} = 'neutral';
            elseif (strcmp(OnlyPrior_data.prior{i}, 'H') && strcmp(OnlyPrior_data.target{i}, 'L')) || ...
                   (strcmp(OnlyPrior_data.prior{i}, 'L') && strcmp(OnlyPrior_data.target{i}, 'H'))
                OnlyPrior_data.congruency{i} = 'incongruent';
            end
        end
        OnlyPrior_params.congruency = OnlyPrior_data.congruency;

        % Select OnlyPrior trials
        iSelect = setStimulusCondition('OnlyPrior');
        iSelect.err = 'c';
        iSelect.SNR = SNR;
        iSelect.congruency = 'congruent';
        OnlyPrior_trials = selectData(OnlyPrior_data.data, OnlyPrior_params, iSelect);

        % Time-frequency analysis
        t_length = 2;
        t_win = 0.15;
        if strcmp(EpochLabel, 'LED')
            t_slidingwin = 0.1:0.01:0.75;  % 100 ms to 750 ms
        else
            t_slidingwin = -0.25:0.01:0.6;  % -250 ms to 600 ms
        end


        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.output = 'powandcsd';
        cfg.foi = 1:1/(t_length/2):100;
        cfg.taper = 'hanning';
        cfg.t_ftimwin = ones(1,length(cfg.foi)) .* t_win;
        cfg.toi = t_slidingwin;
        cfg.keeptrials = 'no';
        cfg.pad = t_length;

        % Time-freq
        tfreq = ft_freqanalysis(cfg, OnlyPrior_trials);
        
        % Theta power
        fidx = tfreq.freq >= 4 & tfreq.freq <= 8;
        tidx = tfreq.time >= 0 & tfreq.time <= 0.2;
        theta_power = squeeze(mean(mean(tfreq.powspctrm(:,fidx,tidx),3),2));
        
        % Store
        theta_power_struct.(EpochLabel) = theta_power;
        labels_struct.(EpochLabel) = tfreq.label;
        
        end 
    
         % === After both epochs ===
        if isfield(theta_power_struct, 'testTone') && isfield(theta_power_struct, 'LED')
            [labels_testTone_out, labels_LED_out] = SpecBalancedRangeBased( ...
                theta_power_struct.testTone, theta_power_struct.LED, ...
                labels_struct.testTone, labels_struct.LED)
        
            % === Split into PFC / AC ===
            testTone_pfc = labels_testTone_out(contains(labels_testTone_out, 'PFC'));
            testTone_ac  = labels_testTone_out(contains(labels_testTone_out, 'AC'));
            
            led_pfc = labels_LED_out(contains(labels_LED_out, 'PFC'));
            led_ac  = labels_LED_out(contains(labels_LED_out, 'AC'));

            % Ensure all are cells 
            if isempty(testTone_pfc), testTone_pfc = cell(0,1); end
            if isempty(testTone_ac),  testTone_ac  = cell(0,1); end
            if isempty(led_pfc),      led_pfc      = cell(0,1); end
            if isempty(led_ac),       led_ac       = cell(0,1); end

            OverlappingLabelsEpoch.OnlyPrior = struct();
            OverlappingLabelsEpoch.OnlyPrior.testTone = struct();
            OverlappingLabelsEpoch.OnlyPrior.LED = struct();

            % === Build output struct ===
            OverlappingLabelsEpoch.OnlyPrior.testTone.PFC = testTone_pfc(:);
            OverlappingLabelsEpoch.OnlyPrior.testTone.AC  = testTone_ac(:);
            OverlappingLabelsEpoch.OnlyPrior.LED.PFC = led_pfc(:);
            OverlappingLabelsEpoch.OnlyPrior.LED.AC  = led_ac(:);

            
            % === Plot theta power distributions ===
            all_vals = [theta_power_struct.testTone; theta_power_struct.LED];
            min_val = min(all_vals);
            max_val = max(all_vals);
            num_bins = 30;
            bin_edges = linspace(min_val, max_val, num_bins + 1);
            
            figure;
            hold on;
            histogram(theta_power_struct.testTone, bin_edges, ...
                'FaceAlpha', 0.5, 'DisplayName', 'OnlyPrior testTone');
            histogram(theta_power_struct.LED, bin_edges, ...
                'FaceAlpha', 0.5, 'DisplayName', 'OnlyPrior LED');
            xlabel('Mean Theta Power');
            ylabel('Channel Count');
            legend;
            title(sprintf('Theta Power: OnlyPrior testTone vs LED (%s %s)', RecDate, Animal));
            
            % === Save plot ===
            saveas(gcf, fullfile(savedir, sprintf('ThetaPower_OnlyPrior_testTone_vs_LED_%s_%s.png', RecDate, Animal)));
            close(gcf);
            
            % === Save struct ===
            save(fullfile(savedir, sprintf('OverlappingLabelsEpoch_OnlyPrior_%s_%s.mat', RecDate, Animal)), 'OverlappingLabelsEpoch');

        else
             fprintf('Skipping %s %s — missing one or both epochs\n', RecDate, Animal);
        
        end
    end














