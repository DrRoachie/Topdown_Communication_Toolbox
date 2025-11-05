
%% Congruency_Connectivity_Test

% This code summarizes the LFP data used to carry out the connectivity study.

% Outputs include:

% Percentage of sessions that have sig. over phase shuffled counterpart 
% Percentage of sessions that have significant evoked potential.
% Average power of those significantly evoked channels 

% This code was expanded in late May and early June to compare the
% spectragrams of the OnlyPrior/OnlyPretone conditions during the testTone
% epoch. The code extracts the chanenls with either overlapping LFP ranges
% or, if non-overlapping, the top and bottom quartile. The channel labels
% of this proceedure will be saved similiarly to the phase shuffle control,
% and those labels will be used to calculate PSI values. Those PSI values
% will then be evaluated both spatially and non-spatially. 

%% Set parameters for the preCue Condition

Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];   % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 

%% Calculate PSI Spectra during the testTone epoch; PriorOnly Trials, Collapsed Across Behavior (Correct and Wrong) and SNR (High and Low)
% ahead of the development of this code, we found that there was not a
% strong influence of SNR on PSI so we collapse across all levels of SNR

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


    % if strcmp(Frequency_Band, 'beta')
    % 
    %     % session_info = {'MrCassius', '190404';
    %     %                 'MrCassius', '190413';
    %     %                 'MrCassius', '190416';
    %     %                 'MrCassius', '190418';
    %     %                 'MrCassius', '190419';
    %     %                 'MrM', '190422';
    %     %                 'MrCassius', '190517';
    %     %                 'MrM', '190525';
    %     %                 'MrM', '190527';
    %     %                 'MrM', '190601';
    %     %                 'MrCassius', '190703';
    %     %                 'MrCassius', '190713';
    %     %                 'MrCassius', '190718';
    %     %                 'MrM', '190719';
    %     %                 'MrCassius', '190723';};



    end



%% Initiate arrays for meta data collection

% Channel-wise count for phase shuffle pass

% Initialize counters for each of the 20 possible channels
Channelwise_OnlyPrior_sigchans_PFC     = zeros(20,1);
Channelwise_OnlyPrior_sigchans_AC      = zeros(20,1);
Channelwise_OnlyPretone_sigchans_PFC   = zeros(20,1);
Channelwise_OnlyPretone_sigchans_AC    = zeros(20,1);

Channelwise_OnlyPrior_Sigevoked_PFC   = zeros(20, 1);
Channelwise_OnlyPrior_Sigevoked_AC    = zeros(20, 1);
Channelwise_OnlyPretone_Sigevoked_PFC = zeros(20, 1);
Channelwise_OnlyPretone_Sigevoked_AC  = zeros(20, 1);

Channelwise_OnlyPrior_LFPpower_PFC   = cell(20, 1);
Channelwise_OnlyPrior_LFPpower_AC    = cell(20, 1);
Channelwise_OnlyPretone_LFPpower_PFC = cell(20, 1);
Channelwise_OnlyPretone_LFPpower_AC  = cell(20, 1);

Channelwise_OnlyPrior_LFPevoked_PFC   = cell(20, 1);
Channelwise_OnlyPrior_LFPevoked_AC    = cell(20, 1);
Channelwise_OnlyPretone_LFPevoked_PFC = cell(20, 1);
Channelwise_OnlyPretone_LFPevoked_AC  = cell(20, 1);

% Store theta-band average evoked responses for each session
AllSessions_ThetaEvoked_OnlyPrior_PFC   = {};
AllSessions_ThetaEvoked_OnlyPretone_PFC = {};
AllSessions_SessionID_PFC               = {};  % keep track of session name
AllSessions_ThetaEvoked_Overlap_PFC = {};


%% saving top-cut LFP labels

% Initialize tfreq collection for OnlyPrior and OnlyPretone
SpecBalanced_Labels_OnlyPrior_PFC   = cell(20, 1);     % each cell = matrix for one PFC channel
SpecBalanced_Labels_OnlyPretone_PFC = cell(20, 1);

SpecBalanced_Labels_OnlyPrior_AC    = cell(20, 1);      % each cell = matrix for one AC channel
SpecBalanced_Labels_OnlyPretone_AC  = cell(20, 1);


%% Initialize outputs for SpecBalancedRangeBased

SpecBalanced_theta.OnlyPrior.PFC = {};
SpecBalanced_theta.OnlyPrior.AC  = {};
SpecBalanced_theta.OnlyPretone.PFC = {};
SpecBalanced_theta.OnlyPretone.AC  = {};





%% Session wise generate the meta data 
       

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 
        
     
        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir           = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\02_Preprocessed', Animal, 'testToneOnset', RecDate);                      % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);    % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\auditory research la\OneDrive\Desktop\Sonia_JunkPile';                 % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
        
        
        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention; only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')  % make folders if they don't exist already
           mkdir(savedir);
        end   

        % Get OnlyPrior Sig. Channels

        OnlyPrior_correct_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'testToneOnset', 'OnlyPrior', 'Correct', Frequency_Band);
        OnlyPrior_correct_files = dir(fullfile(chandir, OnlyPrior_correct_SigChans_fn));

        if  isempty(OnlyPrior_correct_files)

            % generate null output files if there are no sig. channels for Condition_1
            % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', RecDate, Animal, 'testToneOnset', 'OnlyPrior', 'Correct', Frequency_Band));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');

        else % get the labels for the PFC channels and AC channels for correct trials 

        OnlyPrior_correct_SigChans = load(fullfile(chandir, OnlyPrior_correct_SigChans_fn), 'significant_channels');
        OnlyPrior_correct_modifiedCellArray   = regexprep(OnlyPrior_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        OnlyPrior_correct_PFC_chans  =  OnlyPrior_correct_modifiedCellArray(startsWith(OnlyPrior_correct_modifiedCellArray, '*PFC'));
        OnlyPrior_correct_AC_chans   =  OnlyPrior_correct_modifiedCellArray(startsWith(OnlyPrior_correct_modifiedCellArray, '*AC'));

        end

        % Get OnlyPretone Sig. Channels

        OnlyPretone_correct_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'testToneOnset', 'OnlyPretone', 'Correct', Frequency_Band);
        OnlyPretone_correct_files = dir(fullfile(chandir, OnlyPretone_correct_SigChans_fn));

        if  isempty(OnlyPretone_correct_files)

          % generate null output files if there are no sig. channels for Condition_2
          % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

          save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_%s_data.txt', RecDate, Animal,'testToneOnset', 'OnlyPretone', Frequency_Band, 'Correct'));
          fid4 = fopen(save_file_name, 'w');
          fclose('all');

        else

        OnlyPretone_correct_SigChans = load(fullfile(chandir, OnlyPretone_correct_SigChans_fn), 'significant_channels');
        OnlyPretone_correct_modifiedCellArray   = regexprep(OnlyPretone_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        OnlyPretone_correct_PFC_chans  =  OnlyPretone_correct_modifiedCellArray(startsWith(OnlyPretone_correct_modifiedCellArray, '*PFC'));
        OnlyPretone_correct_AC_chans   =  OnlyPretone_correct_modifiedCellArray(startsWith(OnlyPretone_correct_modifiedCellArray, '*AC'));
      
        end        

         % STEP 2: Load in the data, and record numerical values to 'H' and
         % 'L' designation, and code each trials as congruent and
         % incongruent. 

        % extract session information
         baseFileName = sessions.name;
         fullFileName = fullfile(sessions.folder, baseFileName);
         fprintf(1, 'Now reading %s/n', fullFileName);
         OnlyPrior_data = load(fullfile(fullFileName));
         OnlyPretone_data = OnlyPrior_data; % Fork the data to avoid loading the data again...

          % set parameters required for fieldtrip functions and data selection
         OnlyPrior_params.choice        = OnlyPrior_data.choice;
         OnlyPrior_params.err           = OnlyPrior_data.err;
         OnlyPrior_params.pretone       = OnlyPrior_data.pretone;
         OnlyPrior_params.pretoneLength = OnlyPrior_data.pretoneLength;
         OnlyPrior_params.prior         = cell2char(OnlyPrior_data.prior); 
         OnlyPrior_params.SNR           = OnlyPrior_data.SNR;

         % set parameters required for fieldtrip functions and data selection
         OnlyPretone_params.choice        =  OnlyPretone_data.choice;
         OnlyPretone_params.err           =  OnlyPretone_data.err;
         OnlyPretone_params.pretone       =  OnlyPretone_data.pretone;
         OnlyPretone_params.pretoneLength =  OnlyPretone_data.pretoneLength;
         OnlyPretone_params.prior         =  cell2char(OnlyPretone_data.prior); 
         OnlyPretone_params.SNR           =  OnlyPretone_data.SNR;

         % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

         OnlyPrior_data   = ConvertStim(OnlyPrior_data, RecDate);
         OnlyPretone_data = ConvertStim(OnlyPretone_data, RecDate);

         % Adds new field in data structure for OnlyPrior Trials that
         % indicates congruence between target and LED

             % Initialize OnlyPrior_data.congruency as a cell array as a new selection factor for select data function  

             OnlyPrior_data.congruency = cell(size(OnlyPrior_data.prior));  % Preallocate as a cell array
             OnlyPrior_data.pretone    = cellstr(OnlyPrior_data.pretone);

             % Loop through each element in V.prior and V.target

                 for i = 1:length(OnlyPrior_data.prior)
                    if strcmp(OnlyPrior_data.prior{i}, 'H') && strcmp(OnlyPrior_data.target{i}, 'H')
                        OnlyPrior_data.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(OnlyPrior_data.prior{i}, 'L') && strcmp(OnlyPrior_data.target{i}, 'L')
                        OnlyPrior_data.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(OnlyPrior_data.prior{i}, 'N')
                       OnlyPrior_data.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(OnlyPrior_data.prior{i}, 'H') && strcmp(OnlyPrior_data.target{i}, 'L')) || (strcmp(OnlyPrior_data.prior{i}, 'L') && strcmp(OnlyPrior_data.target{i}, 'H'))
                       OnlyPrior_data.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

             OnlyPrior_params.congruency = OnlyPrior_data.congruency;

          % Adds new field in data structure for OnlyPretone Trials that indicates congruence between target and pretone
          % Initialize OnlyPretone_data.congruency as a cell array as a new selection factor for select data function  

            OnlyPretone_data.congruency = cell(size(OnlyPretone_data.pretone));  % Preallocate as a cell array
            OnlyPretone_data.pretone = cellstr(OnlyPretone_data.pretone);

             % Loop through each element in V.pretone and V.target

                 for i = 1:length(OnlyPretone_data.pretone)
                    if strcmp(OnlyPretone_data.pretone{i}, 'H') && strcmp(OnlyPretone_data.target{i}, 'H')
                        OnlyPretone_data.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(OnlyPretone_data.pretone{i}, 'L') && strcmp(OnlyPretone_data.target{i}, 'L')
                        OnlyPretone_data.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(OnlyPretone_data.pretone{i}, 'N')
                        OnlyPretone_data.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(OnlyPretone_data.pretone{i}, 'H') && strcmp(OnlyPretone_data.target{i}, 'L')) || (strcmp(OnlyPretone_data.pretone{i}, 'L') && strcmp(OnlyPretone_data.target{i}, 'H'))
                        OnlyPretone_data.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

            OnlyPretone_params.congruency = OnlyPretone_data.congruency;

            % Select the OnlyPrior trials 
            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = 'c';                                   % Select correct trials
            iSelect.SNR          = SNR;                                   % Select all SNR values
            iSelect.congruency   = 'congruent';                           % Select congruent trials
            OnlyPrior_trials     = selectData(OnlyPrior_data.data, OnlyPrior_params, iSelect);     % 

            % Select the OnlyPretone trials 
            iSelect              = setStimulusCondition('OnlyPretone');
            iSelect.err          = 'c';                                    % Select  correct trials
            iSelect.SNR          = SNR;                                    % Select all SNR values
            iSelect.congruency   = 'congruent';                            % Select congruent trials
            OnlyPretone_trials   = selectData(OnlyPretone_data.data, OnlyPretone_params,iSelect); 
            
            % % Spectral Analysis
            % 
            % cfg                 = [];
            % cfg.method          = 'mtmfft';
            % cfg.output          = 'fourier';
            % cfg.taper           = 'dpss';                                                           
            % cfg.tapsmofrq       = 4;      
            % cfg.pad             = 2.0;
            % cfg.foilim          = [1 100];
            % freq_OnlyPrior      = ft_freqanalysis(cfg, OnlyPrior_trials);
            % freq_OnlyPretone    = ft_freqanalysis(cfg, OnlyPretone_trials);
            % fd_OnlyPrior        = ft_freqdescriptives(cfg, freq_OnlyPrior);
            % fd_OnlyPretone      = ft_freqdescriptives(cfg, freq_OnlyPretone);
            % nTrial_OnlyPrior    = numel(OnlyPrior_trials.trial);
            % nTrial_OnlyPretone  = numel(OnlyPretone_trials.trial);
            % nTaper_OnlyPrior    = size(freq_OnlyPrior.fourierspctrm,1) / nTrial_OnlyPrior;     % number of tapers
            % nTaper_OnlyPretone  = size(freq_OnlyPretone.fourierspctrm,1) / nTrial_OnlyPretone; % number of tapers

            % Time frequency analysis

            t_length     = 2;                  % make the trial length 2 sec with zero-padding
            t_win        = 0.15;                  % 150 ms 
            t_slidingwin = -0.25:0.01:0.6; % 10-ms sliding window
        
            cfg               = [];
            cfg.method        = 'mtmconvol';
            cfg.output        = 'powandcsd';                          % Options:'powandcsd'; 'fourier'
            cfg.foi           = 1:1/(t_length/2):100;
            cfg.taper         = 'hanning';                            % Options: 'hanning'; 'dpss'
            cfg.t_ftimwin     = ones(1,length(cfg.foi)) .* t_win;     % 250-ms sliding window
            cfg.toi           = t_slidingwin; 
            cfg.keeptrials    = 'no';                                 % Options: 'no'; 'yes'
            cfg.pad           = t_length;
        
            tfreq_OnlyPrior     = ft_freqanalysis(cfg, OnlyPrior_trials);
            tfreq_OnlyPretone   = ft_freqanalysis(cfg, OnlyPretone_trials);

            % === Theta Power Calculation ===

            if isfield(tfreq_OnlyPrior, 'freq') && isfield(tfreq_OnlyPrior, 'time') && ...
               isfield(tfreq_OnlyPretone, 'freq') && isfield(tfreq_OnlyPretone, 'time')
            
                % Indices for theta band and early time
                theta_idx_prior   = tfreq_OnlyPrior.freq >= 4 & tfreq_OnlyPrior.freq <= 8;
                time_idx_prior    = tfreq_OnlyPrior.time >= 0 & tfreq_OnlyPrior.time <= 0.2;
            
                theta_idx_pretone = tfreq_OnlyPretone.freq >= 4 & tfreq_OnlyPretone.freq <= 8;
                time_idx_pretone  = tfreq_OnlyPretone.time >= 0 & tfreq_OnlyPretone.time <= 0.2;
        
                    
                % Extract theta power
                theta_power_OnlyPrior   = squeeze(mean(mean(tfreq_OnlyPrior.powspctrm(:, theta_idx_prior, time_idx_prior), 3), 2)); % averages across time, averages across frequency, removes single dimensions
                theta_power_OnlyPretone = squeeze(mean(mean(tfreq_OnlyPretone.powspctrm(:, theta_idx_pretone, time_idx_pretone), 3), 2));
            
                disp('Theta power extraction complete.')
        
            end

            pfc_idx = 1:20;
            ac_idx = 21:40; 
        
            theta_power_OnlyPrior_PFC = theta_power_OnlyPrior(pfc_idx);
            theta_power_OnlyPrior_AC = theta_power_OnlyPrior(ac_idx);
            theta_power_OnlyPretone_PFC = theta_power_OnlyPretone(pfc_idx);
            theta_power_OnlyPretone_AC = theta_power_OnlyPretone(ac_idx);

                                % === Calling SpecBalancedRangeBased ===

                        conditions = {'OnlyPrior', 'OnlyPretone'};
                        
                        % Channel labels (same across sessions)
                        pfc_labels = arrayfun(@(x) sprintf('*PFC_ch%d', x), 3:22, 'UniformOutput', false);
                        ac_labels  = arrayfun(@(x) sprintf('*AC_ch%d', x),  3:22, 'UniformOutput', false);
                        
                        % Struct for theta values
                        theta_power_struct.OnlyPrior.PFC    = theta_power_OnlyPrior_PFC;
                        theta_power_struct.OnlyPrior.AC     = theta_power_OnlyPrior_AC;
                        theta_power_struct.OnlyPretone.PFC  = theta_power_OnlyPretone_PFC;
                        theta_power_struct.OnlyPretone.AC   = theta_power_OnlyPretone_AC;
                        
                        % Loop through conditions and assign SpecBalanced labels
                        for i = 1:length(conditions)
                            cond = conditions{i};
                            theta_a = theta_power_struct.(cond).PFC;
                            theta_b = theta_power_struct.(cond).AC;
                            
                            labels_a = pfc_labels;
                            labels_b = ac_labels;
                        
                            [labels_a_out, labels_b_out] = SpecBalancedRangeBased(theta_a, theta_b, labels_a, labels_b);
                        
                            SpecBalanced_theta.(cond).PFC = [SpecBalanced_theta.(cond).PFC; labels_a_out(:)];
                            SpecBalanced_theta.(cond).AC  = [SpecBalanced_theta.(cond).AC;  labels_b_out(:)];
                        end

                        % === Save SpecBalanced_theta for this session ===
                        save_file_name = fullfile(savedir, sprintf('SpecBalancedTheta_%s_%s.mat', RecDate, Animal));
                        save(save_file_name, 'SpecBalanced_theta'); 
                       
                % === PLOTS ===
        
                % Only Prior Plot
                % bin size modification
                all_vals = [theta_power_OnlyPrior_PFC; theta_power_OnlyPrior_AC];
                min_val = min(all_vals);
                max_val = max(all_vals);
                num_bins = 15;
                bin_edges = linspace(min_val, max_val, num_bins +1);
        
                figure;
                hold on;
                histogram(theta_power_OnlyPrior_PFC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPrior-PFC');
                histogram(theta_power_OnlyPrior_AC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPrior-AC');
                xlabel('Mean Theta Power');
                ylabel('Channel Count');
                legend;
                title('Theta Power Distribution: OnlyPrior Trials PFC versus AC');

                % Save and close
                saveas(gcf, fullfile(savedir, sprintf('ThetaHistogram_OnlyPrior_%s_%s.png', RecDate, Animal)));
                close(gcf);
        
        
                % Only Pretone Plot
                all_vals = [theta_power_OnlyPretone_PFC; theta_power_OnlyPretone_AC];
                min_val = min(all_vals);
                max_val = max(all_vals);
                num_bins = 15;
                bin_edges = linspace(min_val, max_val, num_bins +1);
        
                figure;
                hold on;
                histogram(theta_power_OnlyPretone_PFC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPretone-PFC');
                histogram(theta_power_OnlyPretone_AC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPretone-AC');
                xlabel('Mean Theta Power');
                ylabel('Channel Count');
                legend;
                title('Theta Power Distribution: OnlyPretone Trials PFC versus AC');

                % Save and close
                saveas(gcf, fullfile(savedir, sprintf('ThetaHistogram_OnlyPretone_%s_%s.png', RecDate, Animal)));
                close(gcf);

                % Clear session-specific variables to avoid cross-session contamination
                clear tfreq_OnlyPrior tfreq_OnlyPretone
                clear theta_power_OnlyPrior* theta_power_OnlyPretone*
                clear SpecBalanced_theta_OnlyPrior SpecBalanced_theta_OnlyPretone


                
            % % Get index of baseline period and the evoked period. 
            % 
            % t            = -250:10:600;  % in milliseconds
            % t_sec        = t / 1000;  % from milliseconds to seconds
            % baseline_idx = find(t >= -250 & t < 0);     % pre-stimulus
            % poststim_idx = find(t >= 10 & t <= 250);    % post-stimulus window (example)
            % freq_idx     = find(tfreq_OnlyPrior.freq   >= 1 & tfreq_OnlyPrior.freq   <= 35); % Frequencies 1–35 Hz
            % 
            % nChan_OnlyPrior    = size(tfreq_OnlyPrior.powspctrm, 1);
            % p_vals_OnlyPrior   = nan(nChan_OnlyPrior, 1);  % One test per channel, collapsing over 1–35 Hz
            % nChan_OnlyPretone  = size(tfreq_OnlyPretone.powspctrm, 1);
            % p_vals_OnlyPretone = nan(nChan_OnlyPretone, 1);  % One test per channel, collapsing over 1–35 Hz
        
            % for ch = 1:nChan_OnlyPrior
            % 
            % % Average power over 1–35 Hz for each time bin
            % baseline_power_OnlyPrior = squeeze(mean(tfreq_OnlyPrior.powspctrm(ch, freq_idx, baseline_idx), 1));  % f x t → 1 x t
            % poststim_power_OnlyPrior = squeeze(mean(tfreq_OnlyPrior.powspctrm(ch, freq_idx, poststim_idx), 1));
            % 
            % % Collapse across time
            % baseline_vals_OnlyPrior = mean(baseline_power_OnlyPrior, 2);  % mean power at each frequency
            % poststim_vals_OnlyPrior = mean(poststim_power_OnlyPrior, 2);
            % 
            % % Paired t-test across frequencies (or use nonparametric test if preferred)
            % [~, p] = ttest2(poststim_vals_OnlyPrior, baseline_vals_OnlyPrior);
            % p_vals_OnlyPrior(ch) = p;
            % 
            % end
            % 
            % for ch = 1:nChan_OnlyPretone
            % % Average power over 1–35 Hz for each time bin
            % baseline_power_OnlyPretone = squeeze(mean(tfreq_OnlyPretone.powspctrm(ch, freq_idx, baseline_idx), 1));  % f x t → 1 x t
            % poststim_power_OnlyPretone = squeeze(mean(tfreq_OnlyPretone.powspctrm(ch, freq_idx, poststim_idx), 1));
            % 
            % % Collapse across time
            % baseline_vals_OnlyPretone = mean(baseline_power_OnlyPretone, 2);  % mean power at each frequency
            % poststim_vals_OnlyPretone = mean(poststim_power_OnlyPretone, 2);
            % 
            % % Paired t-test across frequencies (or use nonparametric test if preferred)
            % [~, p] = ttest2(poststim_vals_OnlyPretone, baseline_vals_OnlyPretone);
            % p_vals_OnlyPretone(ch) = p;
            % end
            % 
            % % --- Split into PFC (first 20 channels) and AC (last 20 channels) ---
            % 
            % p_vals_OnlyPrior_PFC    = p_vals_OnlyPrior(1:20);
            % p_vals_OnlyPrior_AC     = p_vals_OnlyPrior(21:40);
            % 
            % p_vals_OnlyPretone_PFC  = p_vals_OnlyPretone(1:20);
            % p_vals_OnlyPretone_AC   = p_vals_OnlyPretone(21:40);
            % 
            % % --- Logical indices of significant channels (p < 0.05) ---
            % 
            % sig_OnlyPrior_PFC    = p_vals_OnlyPrior_PFC   < 0.05;
            % sig_OnlyPrior_AC     = p_vals_OnlyPrior_AC    < 0.05;
            % 
            % sig_OnlyPretone_PFC  = p_vals_OnlyPretone_PFC < 0.05;
            % sig_OnlyPretone_AC   = p_vals_OnlyPretone_AC  < 0.05;
            % 
            %         % === OnlyPrior PFC ===
            % 
            %         for i = 1:length(OnlyPrior_correct_PFC_chans)
            %             chNum = sscanf(OnlyPrior_correct_PFC_chans{i}, '*PFC_ch%d');
            %             idx = chNum - 2;
            %             if idx >= 1 && idx <= 20
            %                 Channelwise_OnlyPrior_sigchans_PFC(idx) = Channelwise_OnlyPrior_sigchans_PFC(idx) + 1;
            %                 Channelwise_OnlyPrior_LFPpower_PFC{idx} = [Channelwise_OnlyPrior_LFPpower_PFC{idx}; fd_OnlyPrior.powspctrm(idx, :)];
            %             end
            %         end
            % 
            %         % === OnlyPrior AC ===
            % 
            %         for i = 1:length(OnlyPrior_correct_AC_chans)
            %             chNum = sscanf(OnlyPrior_correct_AC_chans{i}, '*AC_ch%d');
            %             idx = chNum - 2;
            %             if idx >= 1 && idx <= 20
            %                 Channelwise_OnlyPrior_sigchans_AC(idx) = Channelwise_OnlyPrior_sigchans_AC(idx) + 1;
            %                 Channelwise_OnlyPrior_LFPpower_AC{idx} = [Channelwise_OnlyPrior_LFPpower_AC{idx}; fd_OnlyPrior.powspctrm(idx + 20, :)];
            %             end
            %         end
            % 
            %         % === OnlyPretone PFC ===
            % 
            %         for i = 1:length(OnlyPretone_correct_PFC_chans)
            %             chNum = sscanf(OnlyPretone_correct_PFC_chans{i}, '*PFC_ch%d');
            %             idx = chNum - 2;
            %             if idx >= 1 && idx <= 20
            %                 Channelwise_OnlyPretone_sigchans_PFC(idx) = Channelwise_OnlyPretone_sigchans_PFC(idx) + 1;
            %                 Channelwise_OnlyPretone_LFPpower_PFC{idx} = [Channelwise_OnlyPretone_LFPpower_PFC{idx}; fd_OnlyPretone.powspctrm(idx, :)];
            %             end
            %         end
            % 
            %         % === OnlyPretone AC ===
            % 
            %         for i = 1:length(OnlyPretone_correct_AC_chans)
            %             chNum = sscanf(OnlyPretone_correct_AC_chans{i}, '*AC_ch%d');
            %             idx = chNum - 2;
            %             if idx >= 1 && idx <= 20
            %                 Channelwise_OnlyPretone_sigchans_AC(idx) = Channelwise_OnlyPretone_sigchans_AC(idx) + 1;
            %                 Channelwise_OnlyPretone_LFPpower_AC{idx} = [Channelwise_OnlyPretone_LFPpower_AC{idx}; fd_OnlyPretone.powspctrm(idx + 20, :)];
            %             end
            %         end
            % 
            %         % === OnlyPrior PFC: Significant Evoked Response ===
            % 
            %         for ch = 1:20
            %             if sig_OnlyPrior_PFC(ch)
            %                 Channelwise_OnlyPrior_Sigevoked_PFC(ch) = Channelwise_OnlyPrior_Sigevoked_PFC(ch) + 1;
            %                 Channelwise_OnlyPrior_LFPevoked_PFC{ch} = [Channelwise_OnlyPrior_LFPevoked_PFC{ch}; fd_OnlyPrior.powspctrm(ch, :)];
            %             end
            %         end
            % 
            %         % === OnlyPrior AC: Significant Evoked Response ===
            % 
            %         for ch = 1:20
            %             if sig_OnlyPrior_AC(ch)
            %                 Channelwise_OnlyPrior_Sigevoked_AC(ch) = Channelwise_OnlyPrior_Sigevoked_AC(ch) + 1;
            %                 Channelwise_OnlyPrior_LFPevoked_AC{ch} = [Channelwise_OnlyPrior_LFPevoked_AC{ch}; fd_OnlyPrior.powspctrm(ch + 20, :)];
            %             end
            %         end
            % 
            %         % === OnlyPretone PFC: Significant Evoked Response ===
            % 
            %         for ch = 1:20
            %             if sig_OnlyPretone_PFC(ch)
            %                 Channelwise_OnlyPretone_Sigevoked_PFC(ch) = Channelwise_OnlyPretone_Sigevoked_PFC(ch) + 1;
            %                 Channelwise_OnlyPretone_LFPevoked_PFC{ch} = [Channelwise_OnlyPretone_LFPevoked_PFC{ch}; fd_OnlyPretone.powspctrm(ch, :)];
            %             end
            %         end
            % 
            %         % === OnlyPretone AC: Significant Evoked Response ===
            % 
            %         for ch = 1:20
            %             if sig_OnlyPretone_AC(ch)
            %                 Channelwise_OnlyPretone_Sigevoked_AC(ch) = Channelwise_OnlyPretone_Sigevoked_AC(ch) + 1;
            %                 Channelwise_OnlyPretone_LFPevoked_AC{ch} = [Channelwise_OnlyPretone_LFPevoked_AC{ch}; fd_OnlyPretone.powspctrm(ch + 20, :)];
            %             end
            %         end

    

        end


 %%   
%         % === PLOTS ===
% 
%         % Only Prior Plot
%         % bin size modification
%         all_vals = [theta_power_OnlyPrior_PFC; theta_power_OnlyPrior_AC];
%         min_val = min(all_vals);
%         max_val = max(all_vals);
%         num_bins = 15;
%         bin_edges = linspace(min_val, max_val, num_bins +1);

% 
%         figure;
%         hold on;
%         histogram(theta_power_OnlyPrior_PFC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPrior-PFC');
%         histogram(theta_power_OnlyPrior_AC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPrior-AC');
%         xlabel('Mean Theta Power');
%         ylabel('Channel Count');
%         legend;
%         title('Theta Power Distribution: OnlyPrior Trials PFC versus AC');
% 
% 
%         % Only Pretone Plot
%         all_vals = [theta_power_OnlyPretone_PFC; theta_power_OnlyPretone_AC];
%         min_val = min(all_vals);
%         max_val = max(all_vals);
%         num_bins = 15;
%         bin_edges = linspace(min_val, max_val, num_bins +1);
% 
%         figure;
%         hold on;
%         histogram(theta_power_OnlyPretone_PFC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPretone-PFC');
%         histogram(theta_power_OnlyPretone_AC, bin_edges, 'FaceAlpha', 0.5, 'DisplayName', 'OnlyPretone-AC');
%         xlabel('Mean Theta Power');
%         ylabel('Channel Count');
%         legend;
%         title('Theta Power Distribution: OnlyPretone Trials PFC versus AC');
% %%
    
        % 
        % % Calling the Function
        % 
        % % === Extract region-specific theta power ===  (enter in any amount
        % % of channels)
        % theta_pfc_prior = theta_power_OnlyPrior(1:20);   % Channels 1–20 are PFC
        % theta_ac_prior  = theta_power_OnlyPrior(21:40);  % Channels 21–40 are AC
        % 
        % % === Correct labels for 3:22 for both regions ===
        % pfc_labels = arrayfun(@(ch) sprintf('*PFC_ch%d', ch), 3:22, 'UniformOutput', false);
        % ac_labels  = arrayfun(@(ch) sprintf('*AC_ch%d', ch), 3:22, 'UniformOutput', false);
        % 
        % % === Run range-based balancing ===
        % [SpecBalanced_Labels_OnlyPrior_PFC, SpecBalanced_Labels_OnlyPrior_AC] = ...
        %     SpecBalancedRangeBased(theta_pfc_prior, theta_ac_prior, pfc_labels, ac_labels);
        % 
        % % === Extract region-specific theta power ===
        % theta_pfc_pretone = theta_power_OnlyPretone(1:20);   % Channels 1–20 are PFC
        % theta_ac_pretone  = theta_power_OnlyPretone(21:40);  % Channels 21–40 are AC
        % 
        % % === Correct labels for 3:22 for both regions ===
        % pfc_labels = arrayfun(@(ch) sprintf('*PFC_ch%d', ch), 3:22, 'UniformOutput', false);
        % ac_labels  = arrayfun(@(ch) sprintf('*AC_ch%d', ch), 3:22, 'UniformOutput', false);
        % 
        % % === Run range-based balancing for Pretone ===
        % [SpecBalanced_Labels_OnlyPretone_PFC, SpecBalanced_Labels_OnlyPretone_AC] = ...
        %     SpecBalancedRangeBased(theta_pfc_pretone, theta_ac_pretone, pfc_labels, ac_labels);
        % 
        % 
        % % CHECK: overlapping fake data
        % 
        % % Fake theta power values for Group A (e.g., PFC)
        % theta_a = [2.1, 2.3, 2.5, 2.7, 2.9, 3.0, 3.2, 3.5, 3.7, 3.9];
        % 
        % % Fake theta power values for Group B (e.g., AC)
        % theta_b = [3.3, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.7, 4.8];
        % 
        % % Corresponding channel labels
        % labels_a = arrayfun(@(x) sprintf('A_ch%d', x), 1:10, 'UniformOutput', false);
        % labels_b = arrayfun(@(x) sprintf('B_ch%d', x), 1:10, 'UniformOutput', false);
        % 
        % % Call the function
        % [labels_a_overlap, labels_b_overlap] = SpecBalancedRangeBased(theta_a, theta_b, labels_a, labels_b);
        % 
        % % Display result
        % disp('Channels from A in overlap range:');
        % disp(labels_a_overlap);
        % 
        % disp('Channels from B in overlap range:');
        % disp(labels_b_overlap);


        % % === SPECTRAL BALANCING === ATTEMPT 1
        % function [labels_small, labels_large] = spectral_balance(theta_small, theta_large, labels_small_all, labels_large_all)
        %     pct = 0.25;
        %     n_small = length(theta_small);
        %     n_large = length(theta_large);
        %     n_top = floor(pct * n_small);
        %     n_bot = floor(pct * n_large);
        % 
        %     [~, idx_sorted_small] = sort(theta_small, 'ascend');
        %     [~, idx_sorted_large] = sort(theta_large, 'ascend');
        % 
        %     top_small_idx = idx_sorted_small(end - n_top + 1:end);  % Top 25% from small
        %     bot_large_idx = idx_sorted_large(1:n_bot);              % Bottom 25% from large
        % 
        %     labels_small = cell(n_small, 1);
        %     labels_large = cell(n_large, 1);
        % 
        %     for i = top_small_idx(:)'
        %         if i <= numel(labels_small_all)
        %             labels_small{i} = labels_small_all{i};
        %         end
        %     end
        % 
        %     for i = bot_large_idx(:)'
        %         if i <= numel(labels_large_all)
        %             labels_large{i} = labels_large_all{i};
        %         end
        %     end
        % end
        % 
        % % === ONLY PRIOR ===
        % mean_pfc_prior = mean(theta_power_OnlyPrior_PFC);
        % mean_ac_prior  = mean(theta_power_OnlyPrior_AC);
        % 
        % if mean_pfc_prior < mean_ac_prior
        %     [SpecBalanced_Labels_OnlyPrior_PFC, SpecBalanced_Labels_OnlyPrior_AC] = ...
        %         spectral_balance(theta_power_OnlyPrior_PFC, theta_power_OnlyPrior_AC, pfc_labels, ac_labels);
        % else
        %     [SpecBalanced_Labels_OnlyPrior_AC, SpecBalanced_Labels_OnlyPrior_PFC] = ...
        %         spectral_balance(theta_power_OnlyPrior_AC, theta_power_OnlyPrior_PFC, ac_labels, pfc_labels);
        % end
        % 
        % % === ONLY PRETONE ===
        % mean_pfc_pretone = mean(theta_power_OnlyPretone_PFC);
        % mean_ac_pretone  = mean(theta_power_OnlyPretone_AC);
        % 
        % if mean_pfc_pretone < mean_ac_pretone
        %     [SpecBalanced_Labels_OnlyPretone_PFC, SpecBalanced_Labels_OnlyPretone_AC] = ...
        %         spectral_balance(theta_power_OnlyPretone_PFC, theta_power_OnlyPretone_AC, pfc_labels, ac_labels);
        % else
        %     [SpecBalanced_Labels_OnlyPretone_AC, SpecBalanced_Labels_OnlyPretone_PFC] = ...
        %         spectral_balance(theta_power_OnlyPretone_AC, theta_power_OnlyPretone_PFC, ac_labels, pfc_labels);
        % end
        % 
        % % Remove empty cells
        % SpecBalanced_Labels_OnlyPrior_PFC   = SpecBalanced_Labels_OnlyPrior_PFC(~cellfun(@isempty, SpecBalanced_Labels_OnlyPrior_PFC));
        % SpecBalanced_Labels_OnlyPrior_AC    = SpecBalanced_Labels_OnlyPrior_AC(~cellfun(@isempty, SpecBalanced_Labels_OnlyPrior_AC));
        % SpecBalanced_Labels_OnlyPretone_PFC = SpecBalanced_Labels_OnlyPretone_PFC(~cellfun(@isempty, SpecBalanced_Labels_OnlyPretone_PFC));
        % SpecBalanced_Labels_OnlyPretone_AC  = SpecBalanced_Labels_OnlyPretone_AC(~cellfun(@isempty, SpecBalanced_Labels_OnlyPretone_AC));
        
        % % Only Prior Spectral Balancing 
        % % Goal --> these do not overlap so get top 25% of PFC channels and
        % % bottom 25% of AC channels, then store channel labels into
        % % SpecBalanced
        % 
        % % Define percentage
        % pct = 0.25;
        % 
        % % --- PFC: Top 25% theta power ---
        % n_pfc_OnlyPrior = length(theta_power_OnlyPrior_PFC);         % Number of PFC channels
        % n_top_OnlyPrior = round(pct *  n_pfc_OnlyPrior);             % Channels to keep
        % 
        % [~, sorted_idx_pfc_OnlyPrior] = sort(theta_power_OnlyPrior_PFC, 'ascend');
        % top_pfc_idx = sorted_idx_pfc_OnlyPrior(end - n_top_OnlyPrior + 1 : end);   % Grabbing top indices
        % 
        % for i = top_pfc_idx(:)'
        %     SpecBalanced_Labels_OnlyPrior_PFC{i} = pfc_labels{i};
        % end
        % 
        % % --- AC: Bottom 25% theta power ---
        % n_ac_OnlyPrior = length(theta_power_OnlyPrior_AC);           % Number of AC channels
        % n_bot_OnlyPrior = round(pct * n_ac_OnlyPrior);               % Channels to keep
        % 
        % [~, sorted_idx_ac_OnlyPrior] = sort(theta_power_OnlyPrior_AC, 'ascend');
        % bot_ac_idx = sorted_idx_ac_OnlyPrior(1 : n_bot_OnlyPrior);
        % 
        % for i = bot_ac_idx(:)'
        %     SpecBalanced_Labels_OnlyPrior_AC{i} = ac_labels{i};
        % end
        % 
        % 
        % 
        % % Only Pretone Spectral Balancing 
        % % Goal --> these do not overlap so get top 25% of PFC channels and
        % % bottom 25% of AC channels, then store channel labels into
        % % SpecBalanced
        % 
        % % Define percentage
        % pct = 0.25;
        % 
        % % --- PFC: Top 25% theta power ---
        % n_pfc_OnlyPretone = length(theta_power_OnlyPretone_PFC);         % Number of PFC channels
        % n_top_OnlyPretone = round(pct *  n_pfc_OnlyPretone);             % Channels to keep
        % 
        % [~, sorted_idx_pfc_OnlyPretone] = sort(theta_power_OnlyPretone_PFC, 'ascend');
        % top_pfc_idx = sorted_idx_pfc_OnlyPretone(end - n_top_OnlyPretone + 1 : end);
        % 
        % for i = top_pfc_idx(:)'
        %     SpecBalanced_Labels_OnlyPretone_PFC{i} = pfc_labels{i};
        % end
        % 
        % % --- AC: Bottom 25% theta power ---
        % n_ac_OnlyPretone = length(theta_power_OnlyPretone_AC);           % Number of AC channels
        % n_bot_OnlyPretone = round(pct * n_ac_OnlyPretone);               % Channels to keep
        % 
        % [~, sorted_idx_ac_OnlyPretone] = sort(theta_power_OnlyPretone_AC, 'ascend');
        % bot_ac_idx = sorted_idx_ac_OnlyPretone(1 : n_bot_OnlyPretone);
        % 
        % for i = bot_ac_idx(:)'
        %     SpecBalanced_Labels_OnlyPretone_AC{i} = ac_labels{i};
        % end
        % 
        % disp('Stored top 25% PFC and bottom 25% AC channel labels.');
        % 






        % % === PFC: Spectral Balancing ===

        % % === 1. Define overlap range ===
        % min_overlap = max([min(theta_prior_pfc), min(theta_pretone_pfc)]);
        % max_overlap = min([max(theta_prior_pfc), max(theta_pretone_pfc)]);
        % 
        % % === 2. Find overlapping channels (those within both ranges) ===
        % overlap_idx_prior     = find(theta_prior_pfc >= min_overlap & theta_prior_pfc <= max_overlap);
        % overlap_idx_pretone   = find(theta_pretone_pfc >= min_overlap & theta_pretone_pfc <= max_overlap);
        % 
        % 
        % % === 4. Store overlapping labels ===
        % for i = overlap_idx_prior(:)'
        %     SpecBalanced_Labels_OnlyPrior_PFC{i} = pfc_labels{i};
        % end
        % for i = overlap_idx_pretone(:)'
        %     SpecBalanced_Labels_OnlyPretone_PFC{i} = pfc_labels{i};
        % end
        % 
        % % === 5. NON-overlapping indices ===
        % nonoverlap_idx_prior   = setdiff(1:20, overlap_idx_prior);
        % nonoverlap_idx_pretone = setdiff(1:20, overlap_idx_pretone);
        % 
        % % === 6. Top & bottom quartiles from non-overlapping ===
        % % PRIOR
        % if ~isempty(nonoverlap_idx_prior)
        %     vals = theta_prior_pfc(nonoverlap_idx_prior);
        %     [~, sort_order] = sort(vals);
        %     n = length(sort_order);
        %     top_q = sort_order(round(0.75*n):end);
        %     bot_q = sort_order(1:round(0.25*n));
        %     for j = [bot_q; top_q]'
        %         idx = nonoverlap_idx_prior(j);
        %         SpecBalanced_Labels_OnlyPrior_PFC{idx} = pfc_labels{idx};
        %     end
        % end
        % 
        % % PRETONE
        % if ~isempty(nonoverlap_idx_pretone)
        %     vals = theta_pretone_pfc(nonoverlap_idx_pretone);
        %     [~, sort_order] = sort(vals);
        %     n = length(sort_order);
        %     top_q = sort_order(round(0.75*n):end);
        %     bot_q = sort_order(1:round(0.25*n));
        %     for j = [bot_q; top_q]'
        %         idx = nonoverlap_idx_pretone(j);
        %         SpecBalanced_Labels_OnlyPretone_PFC{idx} = pfc_labels{idx};
        %     end
        % end
        % 
        % 
        % % === AC: Spectral Balancing ===
        % 
        % % Channel index and label setup
        % ac_idx = 3:22;
        % theta_prior_ac   = theta_power_OnlyPrior(ac_idx);     % [20x1]
        % theta_pretone_ac = theta_power_OnlyPretone(ac_idx);   % [20x1]
        % ac_labels = arrayfun(@(ch) sprintf('*AC_ch%d', ch), ac_idx, 'UniformOutput', false);
        % 
        % % 1. Define overlap range
        % min_overlap_ac = max([min(theta_prior_ac), min(theta_pretone_ac)]);
        % max_overlap_ac = min([max(theta_prior_ac), max(theta_pretone_ac)]);
        % 
        % % 2. Find overlapping channels
        % overlap_idx_prior_ac   = find(theta_prior_ac   >= min_overlap_ac & theta_prior_ac   <= max_overlap_ac);
        % overlap_idx_pretone_ac = find(theta_pretone_ac >= min_overlap_ac & theta_pretone_ac <= max_overlap_ac);
        % 
        % % 4. Store overlapping labels
        % for i = overlap_idx_prior_ac(:)'
        %     SpecBalanced_Labels_OnlyPrior_AC{i} = ac_labels{i};
        % end
        % for i = overlap_idx_pretone_ac(:)'
        %     SpecBalanced_Labels_OnlyPretone_AC{i} = ac_labels{i};
        % end
        % 
        % % 5. Non-overlapping indices
        % nonoverlap_idx_prior_ac   = setdiff(1:20, overlap_idx_prior_ac);
        % nonoverlap_idx_pretone_ac = setdiff(1:20, overlap_idx_pretone_ac);
        % 
        % % 6. Top & bottom quartiles from non-overlapping
        % % PRIOR
        % if ~isempty(nonoverlap_idx_prior_ac)
        %     vals = theta_prior_ac(nonoverlap_idx_prior_ac);
        %     [~, sort_order] = sort(vals);
        %     n = length(sort_order);
        %     top_q = sort_order(round(0.75*n):end);
        %     bot_q = sort_order(1:round(0.25*n));
        %     for j = [bot_q; top_q]'
        %         idx = nonoverlap_idx_prior_ac(j);
        %         SpecBalanced_Labels_OnlyPrior_AC{idx} = ac_labels{idx};
        %     end
        % end
        % 
        % % PRETONE
        % if ~isempty(nonoverlap_idx_pretone_ac)
        %     vals = theta_pretone_ac(nonoverlap_idx_pretone_ac);
        %     [~, sort_order] = sort(vals);
        %     n = length(sort_order);
        %     top_q = sort_order(round(0.75*n):end);
        %     bot_q = sort_order(1:round(0.25*n));
        %     for j = [bot_q; top_q]'
        %         idx = nonoverlap_idx_pretone_ac(j);
        %         SpecBalanced_Labels_OnlyPretone_AC{idx} = ac_labels{idx};
        %     end
        % 
        % end
        % 
        % 

    


    % STANDARD DEV METHOD -- computed differences in theta power between
    % conditions, then use std of these differences as a threshold value 

    % === Theta Power Extraction ===
    % theta_idx_prior   = tfreq_OnlyPrior.freq >= 4 & tfreq_OnlyPrior.freq <= 8;
    % time_idx_prior    = tfreq_OnlyPrior.time >= 0 & tfreq_OnlyPrior.time <= 0.2;
    % 
    % theta_idx_pretone = tfreq_OnlyPretone.freq >= 4 & tfreq_OnlyPretone.freq <= 8;
    % time_idx_pretone  = tfreq_OnlyPretone.time >= 0 & tfreq_OnlyPretone.time <= 0.2;
    % 
    % theta_power_OnlyPrior   = squeeze(mean(mean(tfreq_OnlyPrior.powspctrm(:, theta_idx_prior, time_idx_prior), 3), 2));   % [40x1]
    % theta_power_OnlyPretone = squeeze(mean(mean(tfreq_OnlyPretone.powspctrm(:, theta_idx_pretone, time_idx_pretone), 3), 2)); % [40x1]
    % 
    % % === Setup
    % pfc_idx = 3:22;
    % ac_idx  = 3:22;
    % 
    % pfc_labels = arrayfun(@(ch) sprintf('*PFC_ch%d', ch), pfc_idx, 'UniformOutput', false);
    % ac_labels  = arrayfun(@(ch) sprintf('*AC_ch%d', ch), ac_idx,  'UniformOutput', false);
    % 
    % % PFC
    % % === Prep ===
    % pfc_idx = 3:22;
    % theta_pfc_prior   = theta_power_OnlyPrior_PFC;     % [20x1]
    % theta_pfc_pretone = theta_power_OnlyPretone_PFC;   % [20x1]
    % pfc_labels = arrayfun(@(ch) sprintf('*PFC_ch%d', ch), pfc_idx, 'UniformOutput', false);
    % 
    % % === Difference and Overlap
    % diffs = abs(theta_pfc_prior - theta_pfc_pretone);
    % threshold = std(diffs);
    % overlap_mask = diffs < threshold;
    % 
    % % Initialize
    % SpecBalanced_Labels_OnlyPrior_PFC   = cell(20,1);
    % SpecBalanced_Labels_OnlyPretone_PFC = cell(20,1);
    % 
    % % Overlapping channels
    % for i = find(overlap_mask)'
    %     label = pfc_labels{i};
    %     SpecBalanced_Labels_OnlyPrior_PFC{i}   = label;
    %     SpecBalanced_Labels_OnlyPretone_PFC{i} = label;
    % end
    % 
    % % === Non-overlapping → Top + Bottom 25% ===
    % nonoverlap_idx = find(~overlap_mask);
    % n = numel(nonoverlap_idx);
    % 
    % if n >= 4  % Ensure we can extract top/bottom 25%
    %     % --- PRIOR ---
    %     theta_vals = theta_pfc_prior(nonoverlap_idx);
    %     [~, sorted_idx] = sort(theta_vals);
    %     bot_q = sorted_idx(1:round(n*0.25));
    %     top_q = sorted_idx(round(n*0.75):end);
    %     for idx = [bot_q; top_q]'
    %         real_idx = nonoverlap_idx(idx);
    %         SpecBalanced_Labels_OnlyPrior_PFC{real_idx} = pfc_labels{real_idx};
    %     end
    % 
    %     % --- PRETONE ---
    %     theta_vals = theta_pfc_pretone(nonoverlap_idx);
    %     [~, sorted_idx] = sort(theta_vals);
    %     bot_q = sorted_idx(1:round(n*0.25));
    %     top_q = sorted_idx(round(n*0.75):end);
    %     for idx = [bot_q; top_q]'
    %         real_idx = nonoverlap_idx(idx);
    %         SpecBalanced_Labels_OnlyPretone_PFC{real_idx} = pfc_labels{real_idx};
    %     end
    % end
    % 
    % % AC
    % % === Prep ===
    % ac_idx = 3:22;
    % theta_ac_prior   = theta_power_OnlyPrior_AC;     % [20x1]
    % theta_ac_pretone = theta_power_OnlyPretone_AC;   % [20x1]
    % ac_labels = arrayfun(@(ch) sprintf('*AC_ch%d', ch), ac_idx, 'UniformOutput', false);
    % 
    % % === Difference and Overlap
    % diffs_ac = abs(theta_ac_prior - theta_ac_pretone);
    % threshold_ac = std(diffs_ac);
    % overlap_mask_ac = diffs_ac < threshold_ac;
    % 
    % % Initialize
    % SpecBalanced_Labels_OnlyPrior_AC   = cell(20,1);
    % SpecBalanced_Labels_OnlyPretone_AC = cell(20,1);
    % 
    % % Overlapping channels
    % for i = find(overlap_mask_ac)'
    %     label = ac_labels{i};
    %     SpecBalanced_Labels_OnlyPrior_AC{i}   = label;
    %     SpecBalanced_Labels_OnlyPretone_AC{i} = label;
    % end
    % 
    % % === Non-overlapping → Top + Bottom 25% ===
    % nonoverlap_idx_ac = find(~overlap_mask_ac);
    % n_ac = numel(nonoverlap_idx_ac);
    % 
    % if n_ac >= 4  % Make sure we can slice quartiles
    %     % --- PRIOR ---
    %     theta_vals = theta_ac_prior(nonoverlap_idx_ac);
    %     [~, sorted_idx] = sort(theta_vals);
    %     bot_q = sorted_idx(1:round(n_ac*0.25));
    %     top_q = sorted_idx(round(n_ac*0.75):end);
    %     for idx = [bot_q; top_q]'
    %         real_idx = nonoverlap_idx_ac(idx);
    %         SpecBalanced_Labels_OnlyPrior_AC{real_idx} = ac_labels{real_idx};
    %     end
    % 
    %     % --- PRETONE ---
    %     theta_vals = theta_ac_pretone(nonoverlap_idx_ac);
    %     [~, sorted_idx] = sort(theta_vals);
    %     bot_q = sorted_idx(1:round(n_ac*0.25));
    %     top_q = sorted_idx(round(n_ac*0.75):end);
    %     for idx = [bot_q; top_q]'
    %         real_idx = nonoverlap_idx_ac(idx);
    %         SpecBalanced_Labels_OnlyPretone_AC{real_idx} = ac_labels{real_idx};
    %     end
    % end
    % 
    % 
        
    

      % %% Per channel what percentage of sessions pass phase-shuffle evaluation
      % 
      % propSessions_OnlyPrior_PFC     = (Channelwise_OnlyPrior_sigchans_PFC/rd)   * 100;
      % propSessions_OnlyPrior_AC      = (Channelwise_OnlyPrior_sigchans_AC/rd)    * 100;
      % propSessions_OnlyPretone_PFC   = (Channelwise_OnlyPretone_sigchans_PFC/rd) * 100;
      % propSessions_OnlyPretone_AC    = (Channelwise_OnlyPretone_sigchans_AC/rd)  * 100;
      % 
      % %% Per channel what percentage of sessions pass evoked-response evaluation
      % 
      % propSessions_OnlyPrior_evoked_PFC     = (Channelwise_OnlyPrior_Sigevoked_PFC/rd)   * 100;
      % propSessions_OnlyPrior_evoked_AC      = (Channelwise_OnlyPrior_Sigevoked_AC/rd)    * 100;
      % propSessions_OnlyPretone_evoked_PFC   = (Channelwise_OnlyPretone_Sigevoked_PFC/rd) * 100;
      % propSessions_OnlyPretone_evoked_AC    = (Channelwise_OnlyPretone_Sigevoked_AC/rd)  * 100;

  %% For the channels that pass the first filter, what is the average across sessions. 
        % % 
        %  % === Final averaged power spectra across sessions (per channel) ===
        %  Final_ChannelAvg_OnlyPrior_PFC   = nan(20, 199);
        %  Final_ChannelAvg_OnlyPrior_AC    = nan(20, 199);
        % Final_ChannelAvg_OnlyPretone_PFC = nan(20, 199);
        % Final_ChannelAvg_OnlyPretone_AC  = nan(20, 199);
        % 
        %  for ch = 1:20
        %     if ~isempty(Channelwise_OnlyPrior_LFPpower_PFC{ch})
        %         Final_ChannelAvg_OnlyPrior_PFC(ch, :) = mean(Channelwise_OnlyPrior_LFPpower_PFC{ch}, 1);
        %      end
        %      if ~isempty(Channelwise_OnlyPrior_LFPpower_AC{ch})
        %          Final_ChannelAvg_OnlyPrior_AC(ch, :) = mean(Channelwise_OnlyPrior_LFPpower_AC{ch}, 1);
        %     end
        %      if ~isempty(Channelwise_OnlyPretone_LFPpower_PFC{ch})
        %          Final_ChannelAvg_OnlyPretone_PFC(ch, :) = mean(Channelwise_OnlyPretone_LFPpower_PFC{ch}, 1);
        %      end
        %      if ~isempty(Channelwise_OnlyPretone_LFPpower_AC{ch})
        %         Final_ChannelAvg_OnlyPretone_AC(ch, :) = mean(Channelwise_OnlyPretone_LFPpower_AC{ch}, 1);
        %     end
        %  end
%% Plot Proportions of Sig/Evoked Sig Sessions
%       
% figure ()
% 
%       subplot(1,2,1)
%       plot(propSessions_OnlyPrior_AC);
%       ylim([50 100])
%       subplot(1,2,2)
%       plot(propSessions_OnlyPrior_PFC);
%       ylim([50 100])
% 
% figure ()
% 
%       subplot(1,2,1)
%       plot(propSessions_OnlyPrior_AC);
%       ylim([50 100])
%       subplot(1,2,2)
%       plot(propSessions_OnlyPrior_PFC);
%       ylim([50 100])


%%  Plot Spectrograms 
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
% % 
% % Get baseline and apply z-scoring
% tfreq_baseline = GetSpectrogram(datadir, Animal, RecDate, Epoch, Condition, Behavior);  % Use shorter toi inside function for baseline
% tfreq.zpowspctrm = zscore_spectrogram(tfreq, tfreq_baseline);

% % Plot all PFC channels
% plot_spectrogram(tfreq_OnlyPrior, 'PFC');
% 
% % Plot all AC channels
% plot_spectrogram(tfreq_OnlyPretone, 'AC');
% 

