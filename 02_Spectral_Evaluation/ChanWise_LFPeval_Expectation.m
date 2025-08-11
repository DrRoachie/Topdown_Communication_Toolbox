
%% Congruency_Connectivity_Test

% ALERT: This code carries out an analysis on PriorOnly and PretoneOnly
% trials extracted from the full testTone epoch, so the spectral resolution
% will be higher (.5 Hz) versus (4 Hz) but will not just include the
% testTone (likely some pretones and also some of the decision epoch). The
% reason for this error is that CR was originally interested in filtering
% with evoked potentials. To get a baseline period, you need to used the
% preprocessed data not the epoch cut data. This logic then leaked into the
% thinking about the average power observed on each channel. 
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

% Initialize counters for each of the 20 possible channels that pass the
% phase shuffles
Channelwise_OnlyPrior_sigchans_PFC     = zeros(20,1);
Channelwise_OnlyPrior_sigchans_AC      = zeros(20,1);
Channelwise_OnlyPretone_sigchans_PFC   = zeros(20,1);
Channelwise_OnlyPretone_sigchans_AC    = zeros(20,1);

% Initialize counters for each of the 20 possible channels that exhibit and
% acoustic response 
Channelwise_OnlyPrior_Sigevoked_PFC    = zeros(20, 1);
Channelwise_OnlyPrior_Sigevoked_AC     = zeros(20, 1);
Channelwise_OnlyPretone_Sigevoked_PFC  = zeros(20, 1);
Channelwise_OnlyPretone_Sigevoked_AC   = zeros(20, 1);

Channelwise_OnlyPrior_LFPpower_PFC     = cell(20, 1);
Channelwise_OnlyPrior_LFPpower_AC      = cell(20, 1);
Channelwise_OnlyPretone_LFPpower_PFC   = cell(20, 1);
Channelwise_OnlyPretone_LFPpower_AC    = cell(20, 1);

Channelwise_OnlyPrior_LFPevoked_PFC    = cell(20, 1);
Channelwise_OnlyPrior_LFPevoked_AC     = cell(20, 1);
Channelwise_OnlyPretone_LFPevoked_PFC  = cell(20, 1);
Channelwise_OnlyPretone_LFPevoked_AC   = cell(20, 1);

Channelwise_OnlyPrior_LFPavg_PFC       = cell(20,1);
Channelwise_OnlyPrior_LFPavg_AC        = cell(20,1);
Channelwise_OnlyPretone_LFPavg_PFC     = cell(20,1);
Channelwise_OnlyPretone_LFPavg_AC      = cell(20,1);


%% Session wise generate the meta data 
       

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 
        
     
        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir           = fullfile('D:\02_Preprocessed', Animal, 'testToneOnset', RecDate);                      % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('D:\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);    % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\Spectral_Visualization';                 % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
        
      
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

        OnlyPretone_correct_SigChans            = load(fullfile(chandir, OnlyPretone_correct_SigChans_fn), 'significant_channels');
        OnlyPretone_correct_modifiedCellArray   = regexprep(OnlyPretone_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        OnlyPretone_correct_PFC_chans           =  OnlyPretone_correct_modifiedCellArray(startsWith(OnlyPretone_correct_modifiedCellArray, '*PFC'));
        OnlyPretone_correct_AC_chans            =  OnlyPretone_correct_modifiedCellArray(startsWith(OnlyPretone_correct_modifiedCellArray, '*AC'));
      
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
            
            % Spectral Analysis

            cfg                 = [];
            cfg.method          = 'mtmfft';
            cfg.output          = 'fourier';
            cfg.taper           = 'dpss';                                                           
            cfg.tapsmofrq       = 4;      
            cfg.pad             = 2.0;
            cfg.foilim          = [1 100];
            freq_OnlyPrior      = ft_freqanalysis(cfg, OnlyPrior_trials);
            freq_OnlyPretone    = ft_freqanalysis(cfg, OnlyPretone_trials);
            fd_OnlyPrior        = ft_freqdescriptives(cfg, freq_OnlyPrior);
            fd_OnlyPretone      = ft_freqdescriptives(cfg, freq_OnlyPretone);
            nTrial_OnlyPrior    = numel(OnlyPrior_trials.trial);
            nTrial_OnlyPretone  = numel(OnlyPretone_trials.trial);
            nTaper_OnlyPrior    = size(freq_OnlyPrior.fourierspctrm,1) / nTrial_OnlyPrior;     % number of tapers
            nTaper_OnlyPretone  = size(freq_OnlyPretone.fourierspctrm,1) / nTrial_OnlyPretone; % number of tapers

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

            if strcmp(Frequency_Band, 'theta')  
                    % Get index of baseline period and the evoked period. 
                    t            = -250:10:600;  % in milliseconds
                    t_sec        = t / 1000;  % from milliseconds to seconds
                    baseline_idx = find(t >= -250 & t < 0);     % pre-stimulus
                    poststim_idx = find(t >= 10 & t <= 250);    % post-stimulus window, 10 second lag to account for synaptic latency in A1 
                    freq_idx     = find(tfreq_OnlyPrior.freq   >= 4 & tfreq_OnlyPrior.freq   <= 8); % Frequencies 4–8 Hz
            end

            nChan_OnlyPrior        = size(tfreq_OnlyPrior.powspctrm, 1);
            p_vals_OnlyPrior       = nan(nChan_OnlyPrior, 1);                 % One test per channel, collapsing over Frequency_Band
            nChan_OnlyPretone      = size(tfreq_OnlyPretone.powspctrm, 1);
            p_vals_OnlyPretone     = nan(nChan_OnlyPretone, 1);               % One test per channel, collapsing over Frequency_Band
            evoked_avg_OnlyPrior   = nan(nChan_OnlyPrior, 1);
            evoked_avg_OnlyPretone = nan(nChan_OnlyPretone, 1);


            for ch = 1:nChan_OnlyPrior
            % Average power over Frequency_Band for each time bin
            baseline_power_OnlyPrior = squeeze(tfreq_OnlyPrior.powspctrm(ch, freq_idx, baseline_idx));  % f x t → 1 x t
            poststim_power_OnlyPrior = squeeze(tfreq_OnlyPrior.powspctrm(ch, freq_idx, poststim_idx));

            % pool all frequency values over the entire 200 ms epoch (25,
            % 10 ms bins)
            baseline_vals_OnlyPrior = mean(baseline_power_OnlyPrior, 2);  % mean power at each frequency, n = number of ms in condition
            poststim_vals_OnlyPrior = mean(poststim_power_OnlyPrior, 2);
            poststim_avg_OnlyPrior  = mean(poststim_power_OnlyPrior(:));

            % Non-parametric paired test (Wilcoxon signed-rank test)
            p = ranksum(poststim_vals_OnlyPrior, baseline_vals_OnlyPrior);
            p_vals_OnlyPrior(ch) = p;
            
            evoked_avg_OnlyPrior(ch) = poststim_avg_OnlyPrior;
            end


            for ch = 1:nChan_OnlyPretone
            % Average power over 1–35 Hz for each time bin
            baseline_power_OnlyPretone = squeeze(tfreq_OnlyPretone.powspctrm(ch, freq_idx, baseline_idx));  % chan x f x t → f x t
            poststim_power_OnlyPretone = squeeze(tfreq_OnlyPretone.powspctrm(ch, freq_idx, poststim_idx));
            
            % pool all frequency values over the entire 200 ms epoch (25,
            % 10 ms bins)
            baseline_vals_OnlyPretone = mean(baseline_power_OnlyPretone, 2);  % f
            poststim_vals_OnlyPretone = mean(poststim_power_OnlyPretone, 2);
            poststim_avg_OnlyPretone  = mean(poststim_power_OnlyPretone(:));

            % Non-parametric paired test (Wilcoxon signed-rank test)
            p = ranksum(poststim_vals_OnlyPretone, baseline_vals_OnlyPretone);
            p_vals_OnlyPretone(ch) = p;

            evoked_avg_OnlyPretone(ch) = poststim_avg_OnlyPretone;
            end

            % --- Split into PFC (first 20 channels) and AC (last 20 channels) ---

            p_vals_OnlyPrior_PFC         = p_vals_OnlyPrior(1:20);
            p_vals_OnlyPrior_AC          = p_vals_OnlyPrior(21:40);
            p_vals_OnlyPretone_PFC       = p_vals_OnlyPretone(1:20);
            p_vals_OnlyPretone_AC        = p_vals_OnlyPretone(21:40);
            evoked_avg_OnlyPrior_PFC     = evoked_avg_OnlyPrior(1:20);
            evoked_avg_OnlyPrior_AC      = evoked_avg_OnlyPrior(21:40);
            evoked_avg_OnlyPretone_PFC   = evoked_avg_OnlyPretone(1:20);
            evoked_avg_OnlyPretone_AC    = evoked_avg_OnlyPretone(21:40);

            % --- Logical indices of significant channels (p < 0.05) ---

            sig_OnlyPrior_PFC    = p_vals_OnlyPrior_PFC   < 0.05;
            sig_OnlyPrior_AC     = p_vals_OnlyPrior_AC    < 0.05;

            sig_OnlyPretone_PFC  = p_vals_OnlyPretone_PFC < 0.05;
            sig_OnlyPretone_AC   = p_vals_OnlyPretone_AC  < 0.05;

                    % This section takes the phase shuffled sig channels, and
                    % extracts the powerspectra  

                    % === OnlyPrior PFC ===

                    for i = 1:length(OnlyPrior_correct_PFC_chans)
                        chNum = sscanf(OnlyPrior_correct_PFC_chans{i}, '*PFC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_OnlyPrior_sigchans_PFC(idx) = Channelwise_OnlyPrior_sigchans_PFC(idx) + 1; % keeps count
                            Channelwise_OnlyPrior_LFPpower_PFC{idx} = [Channelwise_OnlyPrior_LFPpowerAVG_PFC{idx}; fd_OnlyPrior.powspctrm(idx, :)];
                        end
                    end

                    % === OnlyPrior AC ===

                    for i = 1:length(OnlyPrior_correct_AC_chans)
                        chNum = sscanf(OnlyPrior_correct_AC_chans{i}, '*AC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_OnlyPrior_sigchans_AC(idx) = Channelwise_OnlyPrior_sigchans_AC(idx) + 1;
                            Channelwise_OnlyPrior_LFPpower_AC{idx} = [Channelwise_OnlyPrior_LFPpower_AC{idx}; fd_OnlyPrior.powspctrm(idx + 20, :)];
                        end
                    end

                    % === OnlyPretone PFC ===

                    for i = 1:length(OnlyPretone_correct_PFC_chans)
                        chNum = sscanf(OnlyPretone_correct_PFC_chans{i}, '*PFC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_OnlyPretone_sigchans_PFC(idx) = Channelwise_OnlyPretone_sigchans_PFC(idx) + 1;
                            Channelwise_OnlyPretone_LFPpower_PFC{idx} = [Channelwise_OnlyPretone_LFPpower_PFC{idx}; fd_OnlyPretone.powspctrm(idx, :)];
                        end
                    end

                    % === OnlyPretone AC ===

                    for i = 1:length(OnlyPretone_correct_AC_chans)
                        chNum = sscanf(OnlyPretone_correct_AC_chans{i}, '*AC_ch%d');
                        idx = chNum - 2;
                        if idx >= 1 && idx <= 20
                            Channelwise_OnlyPretone_sigchans_AC(idx) = Channelwise_OnlyPretone_sigchans_AC(idx) + 1;
                            Channelwise_OnlyPretone_LFPpower_AC{idx} = [Channelwise_OnlyPretone_LFPpower_AC{idx}; fd_OnlyPretone.powspctrm(idx + 20, :)];
                        end
                    end

                    % This section takes the channels with sig evoked channels, and
                    % extracts the powerspectra  

                    % === OnlyPrior PFC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPrior_PFC(ch)
                            Channelwise_OnlyPrior_Sigevoked_PFC(ch) = Channelwise_OnlyPrior_Sigevoked_PFC(ch) + 1;
                            Channelwise_OnlyPrior_LFPevoked_PFC{ch} = [Channelwise_OnlyPrior_LFPevoked_PFC{ch}; fd_OnlyPrior.powspctrm(ch, :)];
                            Channelwise_OnlyPrior_LFPavg_PFC{ch}    = [Channelwise_OnlyPrior_LFPavg_PFC{ch}; evoked_avg_OnlyPrior_PFC(ch, :)];
                        end
                    end

                    % === OnlyPrior AC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPrior_AC(ch)
                            Channelwise_OnlyPrior_Sigevoked_AC(ch) = Channelwise_OnlyPrior_Sigevoked_AC(ch) + 1;
                            Channelwise_OnlyPrior_LFPevoked_AC{ch} = [Channelwise_OnlyPrior_LFPevoked_AC{ch}; fd_OnlyPrior.powspctrm(ch + 20, :)];
                            Channelwise_OnlyPrior_LFPavg_AC{ch}    = [Channelwise_OnlyPrior_LFPavg_AC{ch}; evoked_avg_OnlyPrior_AC(ch, :)];
                        end
                    end

                    % === OnlyPretone PFC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPretone_PFC(ch)
                            Channelwise_OnlyPretone_Sigevoked_PFC(ch) = Channelwise_OnlyPretone_Sigevoked_PFC(ch) + 1;
                            Channelwise_OnlyPretone_LFPevoked_PFC{ch} = [Channelwise_OnlyPretone_LFPevoked_PFC{ch}; fd_OnlyPretone.powspctrm(ch, :)];
                            Channelwise_OnlyPretone_LFPavg_PFC{ch}    = [Channelwise_OnlyPretone_LFPavg_PFC{ch}; evoked_avg_OnlyPretone_PFC(ch, :)];
                        end
                    end

                    % === OnlyPretone AC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPretone_AC(ch)
                            Channelwise_OnlyPretone_Sigevoked_AC(ch) = Channelwise_OnlyPretone_Sigevoked_AC(ch) + 1;
                            Channelwise_OnlyPretone_LFPevoked_AC{ch} = [Channelwise_OnlyPretone_LFPevoked_AC{ch}; fd_OnlyPretone.powspctrm(ch + 20, :)];
                            Channelwise_OnlyPretone_LFPavg_AC{ch}    = [Channelwise_OnlyPretone_LFPavg_AC{ch}; evoked_avg_OnlyPretone_AC(ch, :)];
                        end
                    end


      end

      
      %% Per channel what percentage of sessions pass phase-shuffle evaluation

      propSessions_OnlyPrior_PFC     = (Channelwise_OnlyPrior_sigchans_PFC/rd)   * 100;
      propSessions_OnlyPrior_AC      = (Channelwise_OnlyPrior_sigchans_AC/rd)    * 100;
      propSessions_OnlyPretone_PFC   = (Channelwise_OnlyPretone_sigchans_PFC/rd) * 100;
      propSessions_OnlyPretone_AC    = (Channelwise_OnlyPretone_sigchans_AC/rd)  * 100;

    %% OnlyPrior_AC AVG power spectra
    % average across frequencies and sessions in power spectra from fft 
       
    
    % Preallocate 20x1 arrays
    Channelwise_OnlyPrior_LFPpowerAVG_AC     = zeros(20, 1);
    Channelwise_OnlyPrior_LFPpowerAVG_CI_AC  = zeros(20, 1);
    Channelwise_OnlyPrior_LFPpowerSUM_AC     = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = 7:17;
    
    % Loop through each channel
    for ch = 1:20
        data = Channelwise_OnlyPrior_LFPpower_AC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        flat_vals = selected_data(:);                   % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(flat_vals);
        se = std(flat_vals) / sqrt(length(flat_vals));
        ci = 1.96 * se;
    
        % Compute sum of all extracted values
        total_power = sum(flat_vals);
    
        % Store results
        Channelwise_OnlyPrior_LFPpowerAVG_AC(ch)    = mu;
        Channelwise_OnlyPrior_LFPpowerAVG_CI_AC(ch) = ci;
        Channelwise_OnlyPrior_LFPpowerSUM_AC(ch)    = total_power;
    end

    %% OnlyPretone_AC AVG power spectra
    % average across frequencies and sessions in power spectra from fft 
       
    
    % Preallocate 20x1 arrays
    Channelwise_OnlyPretone_LFPpowerAVG_AC     = zeros(20, 1);
    Channelwise_OnlyPretone_LFPpowerAVG_CI_AC  = zeros(20, 1);
    Channelwise_OnlyPretone_LFPpowerSUM_AC     = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = 7:17;
    
    % Loop through each channel
    for ch = 1:20
        data = Channelwise_OnlyPretone_LFPpower_AC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        flat_vals = selected_data(:);                   % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(flat_vals);
        se = std(flat_vals) / sqrt(length(flat_vals));
        ci = 1.96 * se;
    
        % Compute sum of all extracted values
        total_power = sum(flat_vals);
    
        % Store results
        Channelwise_OnlyPretone_LFPpowerAVG_AC(ch)    = mu;
        Channelwise_OnlyPretone_LFPpowerAVG_CI_AC(ch) = ci;
        Channelwise_OnlyPretone_LFPpowerSUM_AC(ch)    = total_power;
    end


    %% OnlyPrior_PFC AVG power spectra
    % average across frequencies and sessions in power spectra from fft 
       
    
    % Preallocate 20x1 arrays
    Channelwise_OnlyPrior_LFPpowerAVG_PFC     = zeros(20, 1);
    Channelwise_OnlyPrior_LFPpowerAVG_CI_PFC  = zeros(20, 1);
    Channelwise_OnlyPrior_LFPpowerSUM_PFC     = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = 7:17;
    
    % Loop through each channel
    for ch = 1:20
        data = Channelwise_OnlyPrior_LFPpower_PFC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        flat_vals = selected_data(:);                   % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(flat_vals);
        se = std(flat_vals) / sqrt(length(flat_vals));
        ci = 1.96 * se;
    
        % Compute sum of all extracted values
        total_power = sum(flat_vals);
    
        % Store results
        Channelwise_OnlyPrior_LFPpowerAVG_PFC(ch)    = mu;
        Channelwise_OnlyPrior_LFPpowerAVG_CI_PFC(ch) = ci;
        Channelwise_OnlyPrior_LFPpowerSUM_PFC(ch)    = total_power;
    end

    %% OnlyPretone_PFC AVG power spectra
    % average across frequencies and sessions in power spectra from fft 
       
    
    % Preallocate 20x1 arrays
    Channelwise_OnlyPretone_LFPpowerAVG_PFC     = zeros(20, 1);
    Channelwise_OnlyPretone_LFPpowerAVG_CI_PFC  = zeros(20, 1);
    Channelwise_OnlyPretone_LFPpowerSUM_PFC     = zeros(20, 1);
    
    % Define the column range for 4–8 Hz (7 to 17 for 0.5 Hz bins)
    freq_cols = 7:17;
    
    % Loop through each channel
    for ch = 1:20
        data = Channelwise_OnlyPretone_LFPpower_PFC{ch};   % n x 199
        selected_data = data(:, freq_cols);             % n x 11
        flat_vals = selected_data(:);                   % all values as a vector
    
        % Compute mean and 95% CI
        mu = mean(flat_vals);
        se = std(flat_vals) / sqrt(length(flat_vals));
        ci = 1.96 * se;
    
        % Compute sum of all extracted values
        total_power = sum(flat_vals);
    
        % Store results
        Channelwise_OnlyPretone_LFPpowerAVG_PFC(ch)    = mu;
        Channelwise_OnlyPretone_LFPpowerAVG_CI_PFC(ch) = ci;
        Channelwise_OnlyPretone_LFPpowerSUM_PFC(ch)    = total_power;
    end


%%
 
% Define custom RGB colors
blue   = [0.2 0.4 1];      % OnlyPrior
orange = [1 0.5 0];        % OnlyPretone

% ----- AC -----
figure();
hold on;
plot(propSessions_OnlyPrior_AC,  'Color', blue,   'LineWidth', 2);
plot(propSessions_OnlyPretone_AC,'Color', orange, 'LineWidth', 2);
ylim([50 100]);
xlabel('Channel');
ylabel('Proportion of Sessions (%)');
title('Auditory Cortex (AC)');
%legend({'OnlyPrior', 'OnlyPretone'});
hold off;

% ----- PFC -----
figure();
hold on;
plot(propSessions_OnlyPrior_PFC,  'Color', blue,   'LineWidth', 2);
plot(propSessions_OnlyPretone_PFC,'Color', orange, 'LineWidth', 2);
ylim([50 100]);
xlabel('Channel');
ylabel('Proportion of Sessions (%)');
title('Prefrontal Cortex (PFC)');
%legend({'OnlyPrior', 'OnlyPretone'});
hold off;



%%
% Define custom RGB colors
blue   = [0.2 0.4 1];      % OnlyPrior
orange = [1 0.5 0];        % OnlyPretone


% --- Compute AC bounds ---
Channelwise_OnlyPrior_LFPpowerAVG_ACupper    = Channelwise_OnlyPrior_LFPpowerAVG_AC  + Channelwise_OnlyPrior_LFPpowerAVG_CI_AC;
Channelwise_OnlyPrior_LFPpowerAVG_AClower    = Channelwise_OnlyPrior_LFPpowerAVG_AC  - Channelwise_OnlyPrior_LFPpowerAVG_CI_AC;
Channelwise_OnlyPretone_LFPpowerAVG_ACupper  = Channelwise_OnlyPretone_LFPpowerAVG_AC + Channelwise_OnlyPretone_LFPpowerAVG_CI_AC;
Channelwise_OnlyPretone_LFPpowerAVG_AClower  = Channelwise_OnlyPretone_LFPpowerAVG_AC - Channelwise_OnlyPretone_LFPpowerAVG_CI_AC;

% --- Compute PFC bounds ---
Channelwise_OnlyPrior_LFPpowerAVG_PFCupper   = Channelwise_OnlyPrior_LFPpowerAVG_PFC  + Channelwise_OnlyPrior_LFPpowerAVG_CI_PFC;
Channelwise_OnlyPrior_LFPpowerAVG_PFClower   = Channelwise_OnlyPrior_LFPpowerAVG_PFC  - Channelwise_OnlyPrior_LFPpowerAVG_CI_PFC;
Channelwise_OnlyPretone_LFPpowerAVG_PFCupper = Channelwise_OnlyPretone_LFPpowerAVG_PFC + Channelwise_OnlyPretone_LFPpowerAVG_CI_PFC;
Channelwise_OnlyPretone_LFPpowerAVG_PFClower = Channelwise_OnlyPretone_LFPpowerAVG_PFC - Channelwise_OnlyPretone_LFPpowerAVG_CI_PFC;

% --- Plot: AC ---
figure();
hold on;
fill([1:20, fliplr(1:20)], [Channelwise_OnlyPrior_LFPpowerAVG_ACupper',   fliplr(Channelwise_OnlyPrior_LFPpowerAVG_AClower')],   blue,   'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([1:20, fliplr(1:20)], [Channelwise_OnlyPretone_LFPpowerAVG_ACupper', fliplr(Channelwise_OnlyPretone_LFPpowerAVG_AClower')], orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Channelwise_OnlyPrior_LFPpowerAVG_AC,   'Color', blue,   'LineWidth', 2);
plot(Channelwise_OnlyPretone_LFPpowerAVG_AC, 'Color', orange, 'LineWidth', 2);
xlabel('Channel');
ylabel('Mean LFP Power (4–8 Hz)');
title('Auditory Cortex (AC)');
ylim([175 1300]);
%legend({'OnlyPrior CI', 'OnlyPretone CI', 'OnlyPrior', 'OnlyPretone'});
hold off;

% --- Plot: PFC ---
figure();
hold on;
fill([1:20, fliplr(1:20)], [Channelwise_OnlyPrior_LFPpowerAVG_PFCupper',   fliplr(Channelwise_OnlyPrior_LFPpowerAVG_PFClower')],   blue,   'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([1:20, fliplr(1:20)], [Channelwise_OnlyPretone_LFPpowerAVG_PFCupper', fliplr(Channelwise_OnlyPretone_LFPpowerAVG_PFClower')], orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(Channelwise_OnlyPrior_LFPpowerAVG_PFC,   'Color', blue,   'LineWidth', 2);
plot(Channelwise_OnlyPretone_LFPpowerAVG_PFC, 'Color', orange, 'LineWidth', 2);
xlabel('Channel');
ylabel('Mean LFP Power (4–8 Hz)');
title('Prefrontal Cortex (PFC)');
ylim([175 1300]);
%legend({'OnlyPrior CI', 'OnlyPretone CI', 'OnlyPrior', 'OnlyPretone'});
hold off;


% Define custom RGB colors (use same as previous figures)
blue   = [0.2 0.4 1];      % OnlyPrior
orange = [1 0.5 0];        % OnlyPretone

% ----- AC: OnlyPrior vs. OnlyPretone -----
figure();
hold on;
scatter(propSessions_OnlyPrior_AC,   Channelwise_OnlyPrior_LFPpowerAVG_AC,   60, 'filled', 'MarkerFaceColor', blue);
scatter(propSessions_OnlyPretone_AC, Channelwise_OnlyPretone_LFPpowerAVG_AC, 60, 'filled', 'MarkerFaceColor', orange);
xlabel('Proportion of Sessions (%)');
ylabel('Mean LFP Power (4–8 Hz)');
title('Auditory Cortex (AC)');
%legend({'OnlyPrior', 'OnlyPretone'}, 'Location', 'best');
grid on;
xlim([50 100]);
ylim([175 1300]);
hold off;

% ----- PFC: OnlyPrior vs. OnlyPretone -----
figure();
hold on;
scatter(propSessions_OnlyPrior_PFC,   Channelwise_OnlyPrior_LFPpowerAVG_PFC,   60, 'filled', 'MarkerFaceColor', blue);
scatter(propSessions_OnlyPretone_PFC, Channelwise_OnlyPretone_LFPpowerAVG_PFC, 60, 'filled', 'MarkerFaceColor', orange);
xlabel('Proportion of Sessions (%)');
ylabel('Mean LFP Power (4–8 Hz)');
title('Prefrontal Cortex (PFC)');
%legend({'OnlyPrior', 'OnlyPretone'}, 'Location', 'best');
grid on;
xlim([50 100]);
ylim([175 400]);
hold off;


%% Statistics

% AC - OnlyPrior
[rho_AC_Prior, p_AC_Prior] = corr(propSessions_OnlyPrior_AC, Channelwise_OnlyPrior_LFPpowerAVG_AC, 'Type', 'Spearman');

% AC - OnlyPretone
[rho_AC_Pretone, p_AC_Pretone] = corr(propSessions_OnlyPretone_AC, Channelwise_OnlyPretone_LFPpowerAVG_AC, 'Type', 'Spearman');

% PFC - OnlyPrior
[rho_PFC_Prior, p_PFC_Prior] = corr(propSessions_OnlyPrior_PFC, Channelwise_OnlyPrior_LFPpowerAVG_PFC, 'Type', 'Spearman');

% PFC - OnlyPretone
[rho_PFC_Pretone, p_PFC_Pretone] = corr(propSessions_OnlyPretone_PFC, Channelwise_OnlyPretone_LFPpowerAVG_PFC, 'Type', 'Spearman');




%  %% Per channel what percentage of sessions pass evoked-response evaluation
% 
%       propSessions_OnlyPrior_evoked_PFC     = (Channelwise_OnlyPrior_Sigevoked_PFC/rd)   * 100;
%       propSessions_OnlyPrior_evoked_AC      = (Channelwise_OnlyPrior_Sigevoked_AC/rd)    * 100;
%       propSessions_OnlyPretone_evoked_PFC   = (Channelwise_OnlyPretone_Sigevoked_PFC/rd) * 100;
%       propSessions_OnlyPretone_evoked_AC    = (Channelwise_OnlyPretone_Sigevoked_AC/rd)  * 100;
% 
% 
% %%%PICK UP HERE: The same logic can be used for evoked sessions, 
% 
% 
% %% Plot average powerspectra of significant phase-evokedsessions that have a signficant evoked response as a function of condition/area
% 
%      figure ()
%       hold on 
%       plot(propSessions_OnlyPrior_evoked_AC);
%       ylim([50 100])
%       plot(propSessions_OnlyPretone_evoked_AC);
%       ylim([50 100])
%       hold off
% 
% 
%  figure()
%      hold on 
%      plot(propSessions_OnlyPrior_evoked_PFC);
%      ylim([50 100])
%      plot(propSessions_OnlyPretone_evoked_PFC);
%      ylim([50 100])
%      hold off


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
