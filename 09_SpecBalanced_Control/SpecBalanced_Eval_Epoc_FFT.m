%% SpecBalanced_FFT
% 
% This code computes spectrally balanced PFC and AC channels for connectivity analysis.
% 
% Summary:
% - Processes LED (preCue) and testTone epochs using OnlyPrior trials (correct, congruent).
% - Loads Epoc_Cut data: 200 ms time windows that match connectivity epochs.
% - Computes the Fourier spectrum using mtmfft (no time resolution, just frequency power).
% - Extracts mean theta band power (4–8 Hz) for each channel.
% - Uses SpecBalancedRangeBased to identify overlapping and non-overlapping channels.
% - Generates and saves histograms of theta power distributions.
% - Saves overlap labels for use in downstream connectivity analysis (PSI, coherence).
% 
% Output:
% - OverlappingLabelsEpoch_OnlyPrior_<RecDate>_<Animal>.mat
% - Theta power distribution plots for each session

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
                datadir = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\04_Epoc_Cut', Animal, EpochFolderName, RecDate);
                chandir = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_epoch', Animal, EpochFolderName, RecDate);
                savedir = 'C:\Users\auditory research la\OneDrive\Desktop\SpecBalanced_Eval_Epoch_FFT2';
                % savedir = fullfile(savedir, RecDate, Animal, EpochLabel);
                savedir = fullfile(savedir, RecDate, Animal, 'OnlyPrior');
        
        
                if ~exist(savedir, 'dir')
                    mkdir(savedir);
                end
        
        
                % Load data
        
                sessions = dir(fullfile(datadir,'*.mat'));
                addpath(genpath(datadir));
                addpath(genpath(chandir));
        
                % Assume one file per folder, or pick the first .mat file
                baseFileName = sessions(1).name;
                fullFileName = fullfile(sessions(1).folder, baseFileName);
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
        
                % Fourier Analysis 
                % No time-frequency windowing needed; mtmfft computes spectrum over full epoch window
                cfg = [];
                cfg.method = 'mtmfft';
                cfg.taper = 'dpss';
                cfg.output = 'pow';      % use pow since power is needed for balancing, not fourier coefficients
                cfg.tapsmofrq = 4; 
                cfg.pad = 1;
                cfg.foilim = [0 100];
                
                freq = ft_freqanalysis(cfg, OnlyPrior_trials);
                
                % Theta power
                theta_idx = freq.freq >= 4 & freq.freq <= 8;
                theta_power = mean(freq.powspctrm(:, theta_idx), 2);  % channels × 1
                
                % Store
                theta_power_struct.(EpochLabel) = theta_power;
                labels_struct.(EpochLabel) = freq.label;
                
                end 
 
         % === After both epochs ===
        if isfield(theta_power_struct, 'testTone') && isfield(theta_power_struct, 'LED')
            [labels_testTone_out, labels_LED_out] = SpecBalancedRangeBased( ...
                theta_power_struct.testTone, theta_power_struct.LED, ...
                labels_struct.testTone, labels_struct.LED);
        
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
            saveas(gcf, fullfile(savedir, sprintf('ThetaPower_OnlyPrior_testTone_vs_LED_FFT_%s_%s.png', RecDate, Animal)));
            close(gcf);
            
            % === Save struct ===
            save(fullfile(savedir, sprintf('OverlappingLabelsEpoch_OnlyPrior_FFT_%s_%s.mat', RecDate, Animal)), 'OverlappingLabelsEpoch');

        else
             fprintf('Skipping %s %s — missing one or both epochs\n', RecDate, Animal);
        
        end
    end
