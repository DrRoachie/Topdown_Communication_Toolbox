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
                        'MrCassius', '190725';}

    end



%% Initiate arrays for meta data collection

% Channel-wise count for phase shuffle pass

Channelwise_OnlyPrior_Sigevoked_PFC   = zeros(20, 1);
Channelwise_OnlyPrior_Sigevoked_AC    = zeros(20, 1);
Channelwise_OnlyPretone_Sigevoked_PFC = zeros(20, 1);
Channelwise_OnlyPretone_Sigevoked_AC  = zeros(20, 1);

Channelwise_OnlyPrior_LFPevoked_PFC   = cell(20, 1);
Channelwise_OnlyPrior_LFPevoked_AC    = cell(20, 1);
Channelwise_OnlyPretone_LFPevoked_PFC = cell(20, 1);
Channelwise_OnlyPretone_LFPevoked_AC  = cell(20, 1);

% === Containers to retain TF results per-session ===
nSessions = size(session_info, 1);

% Per-session FieldTrip outputs
tfreq_OnlyPrior_all   = cell(nSessions, 1);
tfreq_OnlyPretone_all = cell(nSessions, 1);

% --- Session-level metadata containers ---
Session_OnlyPrior_SigCount_PFC   = zeros(nSessions,1);
Session_OnlyPrior_SigCount_AC    = zeros(nSessions,1);
Session_OnlyPrior_SigCount_All   = zeros(nSessions,1);

Session_OnlyPretone_SigCount_PFC = zeros(nSessions,1);
Session_OnlyPretone_SigCount_AC  = zeros(nSessions,1);
Session_OnlyPretone_SigCount_All = zeros(nSessions,1);

% Evoked magnitude (post - pre), averaged within session
%   *_meanSig = mean over significant channels only
%   *_meanAll = mean over all channels (useful sanity check)
Session_OnlyPrior_Effect_PFC_meanSig    = nan(nSessions,1);
Session_OnlyPrior_Effect_AC_meanSig     = nan(nSessions,1);
Session_OnlyPrior_Effect_All_meanSig    = nan(nSessions,1);
Session_OnlyPrior_Effect_PFC_meanAll    = nan(nSessions,1);
Session_OnlyPrior_Effect_AC_meanAll     = nan(nSessions,1);
Session_OnlyPrior_Effect_All_meanAll    = nan(nSessions,1);

Session_OnlyPretone_Effect_PFC_meanSig  = nan(nSessions,1);
Session_OnlyPretone_Effect_AC_meanSig   = nan(nSessions,1);
Session_OnlyPretone_Effect_All_meanSig  = nan(nSessions,1);
Session_OnlyPretone_Effect_PFC_meanAll  = nan(nSessions,1);
Session_OnlyPretone_Effect_AC_meanAll   = nan(nSessions,1);
Session_OnlyPretone_Effect_All_meanAll  = nan(nSessions,1);




%% Session wise generate the meta data 
       

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 
        
     
        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir           = fullfile('\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\02_Preprocessed_bipolar', Animal, 'testToneOnset', RecDate);             % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_bipolar', Animal, 'testTone', RecDate);    % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Desktop\Figure_Spectral_Evaluation';                                                           % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
        
        
        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention; only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')  % make folders if they don't exist already
           mkdir(savedir);
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
            iSelect.congruency   = [];                           % Select congruent trials
            OnlyPrior_trials     = selectData(OnlyPrior_data.data, OnlyPrior_params, iSelect);     % 

            % Select the OnlyPretone trials 
            iSelect              = setStimulusCondition('OnlyPretone');
            iSelect.err          = 'c';                                    % Select  correct trials
            iSelect.SNR          = SNR;                                    % Select all SNR values
            iSelect.congruency   = [];                            % Select congruent trials
            OnlyPretone_trials   = selectData(OnlyPretone_data.data, OnlyPretone_params,iSelect);

            % fourier analysis

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
            t_win        = 0.15;               % 150 ms 
            t_slidingwin = -0.25:0.01:0.6;     % 10-ms sliding window
        
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

            % Save per-session TF results
            tfreq_OnlyPrior_all{rd}   = tfreq_OnlyPrior;
            tfreq_OnlyPretone_all{rd} = tfreq_OnlyPretone;


            % ======= Window + band for STATISTICAL TESTS ONLY =======
            useNonParam = false;   % false = ttest2 (parametric), true = ranksum (non-parametric)
            
            % Time axis in ms (your trials are on a fixed -250:10:600 grid)
            t     = -250:10:600;             % ms
            t_sec = t / 1000;                % sec (kept here if you need it elsewhere)
            
            % Pre/post windows for the test
            baseline_idx = find(t >= -100 & t <  0);     % -100..-10 ms (10 bins)
            poststim_idx = find(t >=   10 & t <= 250);   % 10..250 ms  (25 bins)
            
            % Frequency range for the test (keep original 1–35 Hz)
            freq_idx_prior   = find(tfreq_OnlyPrior  .freq >= 1 & tfreq_OnlyPrior  .freq <= 35);
            freq_idx_pretone = find(tfreq_OnlyPretone.freq >= 1 & tfreq_OnlyPretone.freq <= 35);
            
            % Allocate
            nChan_OnlyPrior    = size(tfreq_OnlyPrior.powspctrm,   1);
            nChan_OnlyPretone  = size(tfreq_OnlyPretone.powspctrm, 1);
            p_vals_OnlyPrior   = nan(nChan_OnlyPrior,   1);
            p_vals_OnlyPretone = nan(nChan_OnlyPretone, 1);
            
            % ---------- OnlyPrior: per-channel tests ----------
            for ch = 1:nChan_OnlyPrior
                % Average across 1–35 Hz → time series of 10-ms bins
                pre_ts  = squeeze(mean(tfreq_OnlyPrior.powspctrm (ch, freq_idx_prior,   baseline_idx), 2)); % [1 x nPre]
                post_ts = squeeze(mean(tfreq_OnlyPrior.powspctrm (ch, freq_idx_prior,   poststim_idx), 2)); % [1 x nPost]
            
                if useNonParam
                    % Wilcoxon rank-sum (non-parametric)
                    p = ranksum(post_ts(:), pre_ts(:));
                else
                    % Two-sample t-test (parametric; unequal n, bins treated as independent)
                    [~, p] = ttest2(post_ts(:), pre_ts(:));  % add 'Vartype','unequal' if desired
                end
                p_vals_OnlyPrior(ch) = p;
            end
            
            % ---------- OnlyPretone: per-channel tests ----------
            for ch = 1:nChan_OnlyPretone
                pre_ts  = squeeze(mean(tfreq_OnlyPretone.powspctrm(ch, freq_idx_pretone, baseline_idx), 2));
                post_ts = squeeze(mean(tfreq_OnlyPretone.powspctrm(ch, freq_idx_pretone, poststim_idx), 2));
            
                if useNonParam
                    p = ranksum(post_ts(:), pre_ts(:));
                else
                    [~, p] = ttest2(post_ts(:), pre_ts(:));
                end
                p_vals_OnlyPretone(ch) = p;
            end
            
            % --- Split by region (unchanged) ---
            p_vals_OnlyPrior_PFC    = p_vals_OnlyPrior(1:20);
            p_vals_OnlyPrior_AC     = p_vals_OnlyPrior(21:40);
            p_vals_OnlyPretone_PFC  = p_vals_OnlyPretone(1:20);
            p_vals_OnlyPretone_AC   = p_vals_OnlyPretone(21:40);

            % --- Logical indices of significant channels (p < 0.05) ---

            sig_OnlyPrior_PFC    = p_vals_OnlyPrior_PFC   < 0.05;
            sig_OnlyPrior_AC     = p_vals_OnlyPrior_AC    < 0.05;
            sig_OnlyPretone_PFC  = p_vals_OnlyPretone_PFC < 0.05;
            sig_OnlyPretone_AC   = p_vals_OnlyPretone_AC  < 0.05;

                    % === OnlyPrior PFC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPrior_PFC(ch)
                            Channelwise_OnlyPrior_Sigevoked_PFC(ch) = Channelwise_OnlyPrior_Sigevoked_PFC(ch) + 1;
                            Channelwise_OnlyPrior_LFPevoked_PFC{ch} = [Channelwise_OnlyPrior_LFPevoked_PFC{ch}; fd_OnlyPrior.powspctrm(ch, :)];
                        end
                    end

                    % === OnlyPrior AC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPrior_AC(ch)
                            Channelwise_OnlyPrior_Sigevoked_AC(ch) = Channelwise_OnlyPrior_Sigevoked_AC(ch) + 1;
                            Channelwise_OnlyPrior_LFPevoked_AC{ch} = [Channelwise_OnlyPrior_LFPevoked_AC{ch}; fd_OnlyPrior.powspctrm(ch + 20, :)];
                        end
                    end

                    % === OnlyPretone PFC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPretone_PFC(ch)
                            Channelwise_OnlyPretone_Sigevoked_PFC(ch) = Channelwise_OnlyPretone_Sigevoked_PFC(ch) + 1;
                            Channelwise_OnlyPretone_LFPevoked_PFC{ch} = [Channelwise_OnlyPretone_LFPevoked_PFC{ch}; fd_OnlyPretone.powspctrm(ch, :)];
                        end
                    end

                    % === OnlyPretone AC: Significant Evoked Response ===

                    for ch = 1:20
                        if sig_OnlyPretone_AC(ch)
                            Channelwise_OnlyPretone_Sigevoked_AC(ch) = Channelwise_OnlyPretone_Sigevoked_AC(ch) + 1;
                            Channelwise_OnlyPretone_LFPevoked_AC{ch} = [Channelwise_OnlyPretone_LFPevoked_AC{ch}; fd_OnlyPretone.powspctrm(ch + 20, :)];
                        end
                    end

                    % ---------- Per-session SIGNIFICANT COUNTS ----------
            Session_OnlyPrior_SigCount_PFC(rd)   = sum(sig_OnlyPrior_PFC);
            Session_OnlyPrior_SigCount_AC(rd)    = sum(sig_OnlyPrior_AC);
            Session_OnlyPrior_SigCount_All(rd)   = Session_OnlyPrior_SigCount_PFC(rd) + Session_OnlyPrior_SigCount_AC(rd);
            
            Session_OnlyPretone_SigCount_PFC(rd) = sum(sig_OnlyPretone_PFC);
            Session_OnlyPretone_SigCount_AC(rd)  = sum(sig_OnlyPretone_AC);
            Session_OnlyPretone_SigCount_All(rd) = Session_OnlyPretone_SigCount_PFC(rd) + Session_OnlyPretone_SigCount_AC(rd);
            
            % ---------- Per-session EVOKED MAGNITUDE (post - pre) ----------
            % Use SAME windows/band you used for the tests: baseline -100..0 ms, post 10..250 ms, 1–35 Hz
            % Build per-channel effects for each condition
            eff_prior   = nan(40,1);
            eff_pretone = nan(40,1);
            
            for ch = 1:40
                % OnlyPrior effect
                pre_ts  = squeeze(mean(tfreq_OnlyPrior .powspctrm(ch, freq_idx_prior,   baseline_idx), 2)); % average over 1–35 Hz → time series
                post_ts = squeeze(mean(tfreq_OnlyPrior .powspctrm(ch, freq_idx_prior,   poststim_idx), 2));
                eff_prior(ch) = mean(post_ts) - mean(pre_ts);   % post - pre
            
                % OnlyPretone effect
                pre_ts2  = squeeze(mean(tfreq_OnlyPretone.powspctrm(ch, freq_idx_pretone, baseline_idx), 2));
                post_ts2 = squeeze(mean(tfreq_OnlyPretone.powspctrm(ch, freq_idx_pretone, poststim_idx), 2));
                eff_pretone(ch) = mean(post_ts2) - mean(pre_ts2);
            end
            
            % Split to regions
            pfc_eff_prior   = eff_prior(1:20);
            ac_eff_prior    = eff_prior(21:40);
            pfc_eff_pret    = eff_pretone(1:20);
            ac_eff_pret     = eff_pretone(21:40);
            
            % Means over SIGNIFICANT channels only (NaN if none significant)
            Session_OnlyPrior_Effect_PFC_meanSig(rd)   = mean(pfc_eff_prior(sig_OnlyPrior_PFC),   'omitnan');
            Session_OnlyPrior_Effect_AC_meanSig(rd)    = mean(ac_eff_prior(sig_OnlyPrior_AC),     'omitnan');
            Session_OnlyPrior_Effect_All_meanSig(rd)   = mean([pfc_eff_prior(sig_OnlyPrior_PFC); ac_eff_prior(sig_OnlyPrior_AC)], 'omitnan');
            
            Session_OnlyPretone_Effect_PFC_meanSig(rd) = mean(pfc_eff_pret(sig_OnlyPretone_PFC),  'omitnan');
            Session_OnlyPretone_Effect_AC_meanSig(rd)  = mean(ac_eff_pret(sig_OnlyPretone_AC),    'omitnan');
            Session_OnlyPretone_Effect_All_meanSig(rd) = mean([pfc_eff_pret(sig_OnlyPretone_PFC); ac_eff_pret(sig_OnlyPretone_AC)], 'omitnan');
            
            % (Optional) Means over ALL channels (regardless of significance)
            Session_OnlyPrior_Effect_PFC_meanAll(rd)   = mean(pfc_eff_prior, 'omitnan');
            Session_OnlyPrior_Effect_AC_meanAll(rd)    = mean(ac_eff_prior,  'omitnan');
            Session_OnlyPrior_Effect_All_meanAll(rd)   = mean(eff_prior,     'omitnan');
            
            Session_OnlyPretone_Effect_PFC_meanAll(rd) = mean(pfc_eff_pret,  'omitnan');
            Session_OnlyPretone_Effect_AC_meanAll(rd)  = mean(ac_eff_pret,   'omitnan');
            Session_OnlyPretone_Effect_All_meanAll(rd) = mean(eff_pretone,   'omitnan');



        end


%% Grand-average across sessions WITHOUT collapsing channels

firstIdx = find(~cellfun(@isempty, tfreq_OnlyPretone_all), 1, 'first');
assert(~isempty(firstIdx), 'No TF data found.');
ref    = tfreq_OnlyPretone_all{firstIdx};
sumPow = zeros(size(ref.powspctrm), 'like', ref.powspctrm);
nValid = 0;

for k = 1:numel(tfreq_OnlyPretone_all)
    if isempty(tfreq_OnlyPretone_all{k}), continue; end
    A = tfreq_OnlyPretone_all{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow))
        sumPow = sumPow + A;
        nValid = nValid + 1;
    else
        warning('Skipping session %d (size mismatch)', k);
    end
end
assert(nValid>0, 'No sessions matched expected size.');
grandMean_pow_pretone = sumPow / nValid;              % [40 x F x T]

f      = ref.freq(:);
tms    = ref.time(:) * 1000;
labels = ref.label;                           % 40x1 cell (FieldTrip)


%% Grand-average across sessions WITHOUT collapsing channels (OnlyPrior)
firstIdx_prior = find(~cellfun(@isempty, tfreq_OnlyPrior_all), 1, 'first');
assert(~isempty(firstIdx_prior), 'No OnlyPrior TF data found.');
ref_prior    = tfreq_OnlyPrior_all{firstIdx_prior};
sumPow_prior = zeros(size(ref_prior.powspctrm), 'like', ref_prior.powspctrm);
nValid_prior = 0;

for k = 1:numel(tfreq_OnlyPrior_all)
    if isempty(tfreq_OnlyPrior_all{k}), continue; end
    A = tfreq_OnlyPrior_all{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow_prior))
        sumPow_prior = sumPow_prior + A;
        nValid_prior = nValid_prior + 1;
    else
        warning('OnlyPrior: Skipping session %d (size mismatch)', k);
    end
end
assert(nValid_prior > 0, 'No OnlyPrior sessions matched expected size.');
grandMean_pow_prior = sumPow_prior / nValid_prior;      % [40 x F x T]

f_prior      = ref_prior.freq(:);
tms_prior    = ref_prior.time(:) * 1000;
labels_prior = ref_prior.label;                          % 40x1 cell (FieldTrip)

%% Plot 20 per-channel spectrograms for PFC (1–20) and AC (21–40) — OnlyPrior

plot_region_grid(grandMean_pow_prior, 1:20,  f_prior, tms_prior, labels_prior, 'OnlyPrior · PFC (ch 1–20)', nValid_prior);
plot_region_grid(grandMean_pow_prior, 21:40, f_prior, tms_prior, labels_prior, 'OnlyPrior · AC  (ch 21–40)', nValid_prior);

%% Plot 20 per-channel spectrograms for PFC (1–20) and AC (21–40)

plot_region_grid(grandMean_pow_pretone, 1:20,  f, tms, labels, 'OnlyPretone · PFC (ch 1–20)', nValid);
plot_region_grid(grandMean_pow_pretone, 21:40, f, tms, labels, 'OnlyPretone · AC  (ch 21–40)', nValid);

%% === Across-session averages ===

% (A) On average, how many channels per session were significant?
AvgSigCh_OnlyPrior_PFC   = mean(Session_OnlyPrior_SigCount_PFC);
AvgSigCh_OnlyPrior_AC    = mean(Session_OnlyPrior_SigCount_AC);
AvgSigCh_OnlyPrior_All   = mean(Session_OnlyPrior_SigCount_All);

AvgSigCh_OnlyPretone_PFC = mean(Session_OnlyPretone_SigCount_PFC);
AvgSigCh_OnlyPretone_AC  = mean(Session_OnlyPretone_SigCount_AC);
AvgSigCh_OnlyPretone_All = mean(Session_OnlyPretone_SigCount_All);

% (B) On average, what was the magnitude of the evoked potential (post - pre)?
%     Using significant channels only (swap to *_meanAll if you prefer all channels)
AvgEff_OnlyPrior_PFC_meanSig   = mean(Session_OnlyPrior_Effect_PFC_meanSig,   'omitnan');
AvgEff_OnlyPrior_AC_meanSig    = mean(Session_OnlyPrior_Effect_AC_meanSig,    'omitnan');
AvgEff_OnlyPrior_All_meanSig   = mean(Session_OnlyPrior_Effect_All_meanSig,   'omitnan');

AvgEff_OnlyPretone_PFC_meanSig = mean(Session_OnlyPretone_Effect_PFC_meanSig, 'omitnan');
AvgEff_OnlyPretone_AC_meanSig  = mean(Session_OnlyPretone_Effect_AC_meanSig,  'omitnan');
AvgEff_OnlyPretone_All_meanSig = mean(Session_OnlyPretone_Effect_All_meanSig, 'omitnan');

% (Optional) quick printout
fprintf('\n=== Avg # significant channels per session ===\n');
fprintf('OnlyPrior   PFC: %.2f,  AC: %.2f,  All: %.2f\n',  AvgSigCh_OnlyPrior_PFC,   AvgSigCh_OnlyPrior_AC,   AvgSigCh_OnlyPrior_All);
fprintf('OnlyPretone PFC: %.2f,  AC: %.2f,  All: %.2f\n',  AvgSigCh_OnlyPretone_PFC, AvgSigCh_OnlyPretone_AC, AvgSigCh_OnlyPretone_All);

fprintf('\n=== Avg evoked magnitude (post - pre), significant chans only ===\n');
fprintf('OnlyPrior   PFC: %.4g, AC: %.4g, All: %.4g\n',     AvgEff_OnlyPrior_PFC_meanSig,   AvgEff_OnlyPrior_AC_meanSig,   AvgEff_OnlyPrior_All_meanSig);
fprintf('OnlyPretone PFC: %.4g, AC: %.4g, All: %.4g\n',     AvgEff_OnlyPretone_PFC_meanSig, AvgEff_OnlyPretone_AC_meanSig, AvgEff_OnlyPretone_All_meanSig);

%% === Across-session statistical comparisons (Expectation x Area) ===
useNonParamAcross = true;  % false = paired t-tests, true = Wilcoxon signed-rank

% ---------- (1) Expectation main effect, collapsed across area: COUNTS ----------
X_prior_all   = Session_OnlyPrior_SigCount_All;
X_pretone_all = Session_OnlyPretone_SigCount_All;

if useNonParamAcross
    p_counts_main = signrank(X_pretone_all, X_prior_all);
    test_counts   = 'Wilcoxon signed-rank';
else
    [~, p_counts_main, ~, stats_counts] = ttest(X_pretone_all, X_prior_all);
    % paired effect size (Cohen''s dz)
    diff_counts = X_pretone_all - X_prior_all;
    dz_counts   = mean(diff_counts, 'omitnan') / std(diff_counts, 'omitnan');
    test_counts = sprintf('paired t-test (dz=%.3f, t=%.2f, df=%d)', dz_counts, stats_counts.tstat, stats_counts.df);
end

fprintf('\n[COUNTS] Expectation main effect (OnlyPretone vs OnlyPrior, ALL channels): p=%.4g, %s\n', ...
    p_counts_main, test_counts);

% ---------- (2) Expectation main effect BY REGION: COUNTS ----------
X_prior_pfc   = Session_OnlyPrior_SigCount_PFC;
X_pretone_pfc = Session_OnlyPretone_SigCount_PFC;

X_prior_ac    = Session_OnlyPrior_SigCount_AC;
X_pretone_ac  = Session_OnlyPretone_SigCount_AC;

if useNonParamAcross
    p_counts_pfc = signrank(X_pretone_pfc, X_prior_pfc);
    p_counts_ac  = signrank(X_pretone_ac,  X_prior_ac);
    note_byreg   = 'Wilcoxon signed-rank';
else
    [~, p_counts_pfc, ~, stats_pfc] = ttest(X_pretone_pfc, X_prior_pfc);
    [~, p_counts_ac,  ~, stats_ac ] = ttest(X_pretone_ac,  X_prior_ac);
    diff_pfc = X_pretone_pfc - X_prior_pfc; dz_pfc = mean(diff_pfc,'omitnan')/std(diff_pfc,'omitnan');
    diff_ac  = X_pretone_ac  - X_prior_ac;  dz_ac  = mean(diff_ac, 'omitnan')/std(diff_ac, 'omitnan');
    note_byreg = sprintf('paired t-tests (PFC dz=%.3f, t=%.2f, df=%d; AC dz=%.3f, t=%.2f, df=%d)', ...
        dz_pfc, stats_pfc.tstat, stats_pfc.df, dz_ac, stats_ac.tstat, stats_ac.df);
end

fprintf('[COUNTS] Expectation effect BY REGION: PFC p=%.4g, AC p=%.4g, %s\n', ...
    p_counts_pfc, p_counts_ac, note_byreg);

% ---------- (3) 2×2 RM-ANOVA on COUNTS (Expectation x Area) ----------
T_counts = table( ...
    X_prior_pfc, X_prior_ac, X_pretone_pfc, X_pretone_ac, ...
    'VariableNames', {'Prior_PFC','Prior_AC','Pretone_PFC','Pretone_AC'});

within = table( ...
    categorical({'Prior';   'Prior';    'Pretone';   'Pretone'}), ...
    categorical({'PFC';     'AC';       'PFC';       'AC'}), ...
    'VariableNames', {'Expectation','Area'}, ...
    'RowNames', {'Prior_PFC','Prior_AC','Pretone_PFC','Pretone_AC'});

rm_counts = fitrm(T_counts, 'Prior_PFC-Pretone_AC ~ 1', 'WithinDesign', within);
ran_counts = ranova(rm_counts, 'WithinModel', 'Expectation*Area');

disp('RM-ANOVA (COUNTS):'); disp(ran_counts);

% ---------- (4) Expectation main effect, collapsed across area: EFFECT MAGNITUDE ----------
Y_prior_all_eff   = Session_OnlyPrior_Effect_All_meanSig;
Y_pretone_all_eff = Session_OnlyPretone_Effect_All_meanSig;

valid_eff_all = ~isnan(Y_prior_all_eff) & ~isnan(Y_pretone_all_eff);
Ya = Y_prior_all_eff(valid_eff_all);
Yb = Y_pretone_all_eff(valid_eff_all);

if useNonParamAcross
    p_eff_main = signrank(Yb, Ya);
    test_eff   = 'Wilcoxon signed-rank';
else
    [~, p_eff_main, ~, stats_eff] = ttest(Yb, Ya);
    diff_eff = Yb - Ya;
    dz_eff   = mean(diff_eff,'omitnan')/std(diff_eff,'omitnan');
    test_eff = sprintf('paired t-test (dz=%.3f, t=%.2f, df=%d)', dz_eff, stats_eff.tstat, stats_eff.df);
end

fprintf('\n[EFFECT] Expectation main effect (OnlyPretone vs OnlyPrior, ALL channels; significant chans only): p=%.4g, %s (n=%d sessions)\n', ...
    p_eff_main, test_eff, sum(valid_eff_all));

% ---------- (5) Expectation main effect BY REGION: EFFECT MAGNITUDE ----------
Y_prior_pfc = Session_OnlyPrior_Effect_PFC_meanSig;
Y_prior_ac  = Session_OnlyPrior_Effect_AC_meanSig;
Y_pret_pfc  = Session_OnlyPretone_Effect_PFC_meanSig;
Y_pret_ac   = Session_OnlyPretone_Effect_AC_meanSig;

valid_eff_pfc = ~isnan(Y_prior_pfc) & ~isnan(Y_pret_pfc);
valid_eff_ac  = ~isnan(Y_prior_ac)  & ~isnan(Y_pret_ac);

if useNonParamAcross
    p_eff_pfc = signrank(Y_pret_pfc(valid_eff_pfc), Y_prior_pfc(valid_eff_pfc));
    p_eff_ac  = signrank(Y_pret_ac(valid_eff_ac),   Y_prior_ac(valid_eff_ac));
    note_eff_byreg = 'Wilcoxon signed-rank';
else
    [~, p_eff_pfc, ~, stats_eff_pfc] = ttest(Y_pret_pfc(valid_eff_pfc), Y_prior_pfc(valid_eff_pfc));
    [~, p_eff_ac,  ~, stats_eff_ac ] = ttest(Y_pret_ac(valid_eff_ac),   Y_prior_ac(valid_eff_ac));
    diff_eff_pfc = Y_pret_pfc(valid_eff_pfc) - Y_prior_pfc(valid_eff_pfc);
    diff_eff_ac  = Y_pret_ac(valid_eff_ac)  - Y_prior_ac(valid_eff_ac);
    dz_eff_pfc   = mean(diff_eff_pfc,'omitnan')/std(diff_eff_pfc,'omitnan');
    dz_eff_ac    = mean(diff_eff_ac,'omitnan') /std(diff_eff_ac,'omitnan');
    note_eff_byreg = sprintf('paired t-tests (PFC dz=%.3f, t=%.2f, df=%d; AC dz=%.3f, t=%.2f, df=%d)', ...
        dz_eff_pfc, stats_eff_pfc.tstat, stats_eff_pfc.df, dz_eff_ac, stats_eff_ac.tstat, stats_eff_ac.df);
end

fprintf('[EFFECT] Expectation effect BY REGION (significant chans only): PFC p=%.4g (n=%d), AC p=%.4g (n=%d), %s\n', ...
    p_eff_pfc, sum(valid_eff_pfc), p_eff_ac, sum(valid_eff_ac), note_eff_byreg);

% ---------- (6) 2×2 RM-ANOVA on EFFECT MAGNITUDE (Expectation x Area) ----------
T_eff = table(Y_prior_pfc, Y_prior_ac, Y_pret_pfc, Y_pret_ac, ...
    'VariableNames', {'Prior_PFC','Prior_AC','Pretone_PFC','Pretone_AC'});

% Drop sessions with NaN in any of the four cells
keep = all(~isnan(T_eff{:,:}), 2);
T_eff = T_eff(keep, :);

within_eff = within;  % same design matrix
rm_eff = fitrm(T_eff, 'Prior_PFC-Pretone_AC ~ 1', 'WithinDesign', within_eff);
ran_eff = ranova(rm_eff, 'WithinModel', 'Expectation*Area');

disp('RM-ANOVA (EFFECT magnitude; significant chans only):'); disp(ran_eff);


%% ===== Summary tables: means and SDs across sessions =====

rows = {'PFC'; 'AC'; 'All'};

% ---------- COUNTS ----------
% Session arrays expected:
%   Session_OnlyPrior_SigCount_PFC,  Session_OnlyPrior_SigCount_AC,  Session_OnlyPrior_SigCount_All
%   Session_OnlyPretone_SigCount_PFC,Session_OnlyPretone_SigCount_AC,Session_OnlyPretone_SigCount_All

% Means
m_prior_pfc = mean(Session_OnlyPrior_SigCount_PFC,  'omitnan');
m_prior_ac  = mean(Session_OnlyPrior_SigCount_AC,   'omitnan');
m_prior_all = mean(Session_OnlyPrior_SigCount_All,  'omitnan');

m_pret_pfc  = mean(Session_OnlyPretone_SigCount_PFC,'omitnan');
m_pret_ac   = mean(Session_OnlyPretone_SigCount_AC, 'omitnan');
m_pret_all  = mean(Session_OnlyPretone_SigCount_All,'omitnan');

% SDs
sd_prior_pfc = std(Session_OnlyPrior_SigCount_PFC, 0, 'omitnan');
sd_prior_ac  = std(Session_OnlyPrior_SigCount_AC,  0, 'omitnan');
sd_prior_all = std(Session_OnlyPrior_SigCount_All, 0, 'omitnan');

sd_pret_pfc  = std(Session_OnlyPretone_SigCount_PFC,0,'omitnan');
sd_pret_ac   = std(Session_OnlyPretone_SigCount_AC, 0,'omitnan');
sd_pret_all  = std(Session_OnlyPretone_SigCount_All,0,'omitnan');

% Ns (non-NaN)
n_prior_pfc = sum(~isnan(Session_OnlyPrior_SigCount_PFC));
n_prior_ac  = sum(~isnan(Session_OnlyPrior_SigCount_AC));
n_prior_all = sum(~isnan(Session_OnlyPrior_SigCount_All));

n_pret_pfc  = sum(~isnan(Session_OnlyPretone_SigCount_PFC));
n_pret_ac   = sum(~isnan(Session_OnlyPretone_SigCount_AC));
n_pret_all  = sum(~isnan(Session_OnlyPretone_SigCount_All));

T_counts = table( ...
    [m_prior_pfc; m_prior_ac; m_prior_all], [sd_prior_pfc; sd_prior_ac; sd_prior_all], [n_prior_pfc; n_prior_ac; n_prior_all], ...
    [m_pret_pfc;  m_pret_ac;  m_pret_all],  [sd_pret_pfc;  sd_pret_ac;  sd_pret_all],  [n_pret_pfc;  n_pret_ac;  n_pret_all], ...
    'VariableNames', {'OnlyPrior_mean','OnlyPrior_sd','OnlyPrior_n','OnlyPretone_mean','OnlyPretone_sd','OnlyPretone_n'}, ...
    'RowNames', rows);

disp('--- Counts per session (means ± SD) ---');
disp(T_counts);

% ---------- EFFECT MAGNITUDE (significant channels only) ----------
% Session arrays expected:
%   Session_OnlyPrior_Effect_PFC_meanSig,  Session_OnlyPrior_Effect_AC_meanSig,  Session_OnlyPrior_Effect_All_meanSig
%   Session_OnlyPretone_Effect_PFC_meanSig,Session_OnlyPretone_Effect_AC_meanSig,Session_OnlyPretone_Effect_All_meanSig

m_prior_pfc_e = mean(Session_OnlyPrior_Effect_PFC_meanSig,  'omitnan');
m_prior_ac_e  = mean(Session_OnlyPrior_Effect_AC_meanSig,   'omitnan');
m_prior_all_e = mean(Session_OnlyPrior_Effect_All_meanSig,  'omitnan');

m_pret_pfc_e  = mean(Session_OnlyPretone_Effect_PFC_meanSig,'omitnan');
m_pret_ac_e   = mean(Session_OnlyPretone_Effect_AC_meanSig, 'omitnan');
m_pret_all_e  = mean(Session_OnlyPretone_Effect_All_meanSig,'omitnan');

sd_prior_pfc_e = std(Session_OnlyPrior_Effect_PFC_meanSig, 0, 'omitnan');
sd_prior_ac_e  = std(Session_OnlyPrior_Effect_AC_meanSig,  0, 'omitnan');
sd_prior_all_e = std(Session_OnlyPrior_Effect_All_meanSig, 0, 'omitnan');

sd_pret_pfc_e  = std(Session_OnlyPretone_Effect_PFC_meanSig,0,'omitnan');
sd_pret_ac_e   = std(Session_OnlyPretone_Effect_AC_meanSig, 0,'omitnan');
sd_pret_all_e  = std(Session_OnlyPretone_Effect_All_meanSig,0,'omitnan');

n_prior_pfc_e = sum(~isnan(Session_OnlyPrior_Effect_PFC_meanSig));
n_prior_ac_e  = sum(~isnan(Session_OnlyPrior_Effect_AC_meanSig));
n_prior_all_e = sum(~isnan(Session_OnlyPrior_Effect_All_meanSig));

n_pret_pfc_e  = sum(~isnan(Session_OnlyPretone_Effect_PFC_meanSig));
n_pret_ac_e   = sum(~isnan(Session_OnlyPretone_Effect_AC_meanSig));
n_pret_all_e  = sum(~isnan(Session_OnlyPretone_Effect_All_meanSig));

T_effect = table( ...
    [m_prior_pfc_e; m_prior_ac_e; m_prior_all_e], [sd_prior_pfc_e; sd_prior_ac_e; sd_prior_all_e], [n_prior_pfc_e; n_prior_ac_e; n_prior_all_e], ...
    [m_pret_pfc_e;  m_pret_ac_e;  m_pret_all_e],  [sd_pret_pfc_e;  sd_pret_ac_e;  sd_pret_all_e],  [n_pret_pfc_e;  n_pret_ac_e;  n_pret_all_e], ...
    'VariableNames', {'OnlyPrior_mean','OnlyPrior_sd','OnlyPrior_n','OnlyPretone_mean','OnlyPretone_sd','OnlyPretone_n'}, ...
    'RowNames', rows);

disp('--- Evoked magnitude on significant channels (means ± SD) ---');
disp(T_effect);

% ---------- Optional: write to CSV ----------
% writetable(T_counts, 'Summary_Counts.csv', 'WriteRowNames', true);
% writetable(T_effect, 'Summary_EffectMagnitude_SigOnly.csv', 'WriteRowNames', true);


%% Counts (bar) + Magnitude (box) — no mean diamonds
set(groot,'defaultFigureRenderer','painters');  % safer 2D rendering

% ---------- Session vectors ----------
% COUNTS
AC_prior_counts  = Session_OnlyPrior_SigCount_AC;
AC_pret_counts   = Session_OnlyPretone_SigCount_AC;
PFC_prior_counts = Session_OnlyPrior_SigCount_PFC;
PFC_pret_counts  = Session_OnlyPretone_SigCount_PFC;

% MAGNITUDE (significant channels only)
AC_prior_mag  = Session_OnlyPrior_Effect_AC_meanSig;
AC_pret_mag   = Session_OnlyPretone_Effect_AC_meanSig;
PFC_prior_mag = Session_OnlyPrior_Effect_PFC_meanSig;
PFC_pret_mag  = Session_OnlyPretone_Effect_PFC_meanSig;

% Colors
c_prior = [0.23 0.49 0.85];   % Informative (OnlyPrior)
c_pret  = [0.85 0.45 0.23];   % Non-informative (OnlyPretone)

fig = figure('Color','w','Renderer','painters');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ================= LEFT: COUNTS — grouped bar (mean ± SD) =================
nexttile(1);
data_counts = [ mean(AC_prior_counts,'omitnan')   mean(AC_pret_counts,'omitnan');
                mean(PFC_prior_counts,'omitnan')  mean(PFC_pret_counts,'omitnan') ];
err_counts  = [ std(AC_prior_counts, 0,'omitnan') std(AC_pret_counts, 0,'omitnan');
                std(PFC_prior_counts,0,'omitnan') std(PFC_pret_counts,0,'omitnan') ];

b1 = bar(data_counts, 'grouped'); hold on;
b1(1).FaceColor = c_prior;      % informative
b1(2).FaceColor = c_pret;       % non-informative

for s = 1:size(data_counts,2)
    x = b1(s).XEndPoints;
    errorbar(x, data_counts(:,s), err_counts(:,s), 'k', 'linestyle','none', 'LineWidth',1);
end
set(gca,'XTickLabel',{'AC','vlPFC'});
ylabel('Evoked channels / session');
title('Counts (mean \pm SD)');
legend({'Informative Priors','Stimulus-Based Cues'}, 'Location','northeastoutside');
grid off; box off;

% ================= RIGHT: MAGNITUDE — grouped boxplot (no diamonds) =================
nexttile(2);

yM = [AC_prior_mag; AC_pret_mag; PFC_prior_mag; PFC_pret_mag];
areaM = [repmat("AC",  numel(AC_prior_mag),1);  repmat("AC",  numel(AC_pret_mag),1); ...
         repmat("PFC", numel(PFC_prior_mag),1); repmat("PFC", numel(PFC_pret_mag),1)];
condM = [repmat("Informative",      numel(AC_prior_mag),1);  repmat("Non-informative", numel(AC_pret_mag),1); ...
         repmat("Informative",      numel(PFC_prior_mag),1); repmat("Non-informative", numel(PFC_pret_mag),1)];

% Boxplot grouped by Area with Condition as ColorGroup
boxplot(yM, {areaM,condM}, 'FactorSeparator',1, ...
        'LabelVerbosity','majorminor', 'Symbol','', ...
        'ColorGroup',condM, 'Colors',[c_prior; c_pret]);

ylabel('Evoked magnitude (a.u.)');
title('Magnitude (post - pre)');
grid off; box off;

% Thicken boxplot lines for readability
set(findobj(gca,'Tag','Box'),     'LineWidth',1.4);
set(findobj(gca,'Tag','Median'),  'LineWidth',1.4);
set(findobj(gca,'Tag','Whisker'), 'LineWidth',1.2);

% Legend to indicate color ↔ condition
hold on;
hInf = plot(nan,nan,'-','Color',c_prior,'LineWidth',2);
hNon = plot(nan,nan,'-','Color',c_pret, 'LineWidth',2);
legend([hInf hNon], {'Informative Priors','Stimulus-Based Cues'}, 'Location','northeastoutside');

% Optional: if magnitude spans a huge range, try log scale
% set(gca,'YScale','log');

% Optional: save
% exportgraphics(fig,'counts_bar__magnitude_box_noMeanMarkers.png','Resolution',300);

%% === Pairwise comparisons (AC & PFC only) ===

Pairwise = table;
Pairwise.Metric = [ ...
    "Evoked channels / session"; ...
    "Evoked channels / session"; ...
    "Evoked magnitude (uV^2/Hz; sig chans)"; ...
    "Evoked magnitude (uV^2/Hz; sig chans)"];
Pairwise.Region = ["AC"; "PFC"; "AC"; "PFC"];

% Means & SDs (from your summaries)
Pairwise.Informative_mean = [15.8; 11.6; 1651.3; 144.4];
Pairwise.Informative_sd   = [5.4922; 4.571; 1586.6; 186.21];
Pairwise.NonInform_mean   = [11.0; 8.8857; 189.67; 102.6];
Pairwise.NonInform_sd     = [5.2971; 3.3499; 920.86; 128.75];

% Sample sizes, t(df), p, dz
Pairwise.n_pairs = [35; 35; 34; 35];
Pairwise.t_df    = ["-3.53 (34)"; "-2.97 (34)"; "-4.87 (33)"; "-1.96 (34)"];
Pairwise.p_value = [0.001204; 0.005498; 2.73e-05; 0.05798];
Pairwise.dz      = [-0.597; -0.501; -0.835; -0.332];

disp('Table 1a — Pairwise by region (AC & PFC only):'); disp(Pairwise);

%% === 2x2 RM-ANOVA (Expectation x Area) ===

ANOVA_Counts = table(["Expectation";"Area";"Expectation×Area"], ...
                     [23.57; 13.924; 1.465], ...
                     [2.6466e-05; 6.9408e-04; 0.23456], ...
                     'VariableNames', {'Effect','F','p'});
ANOVA_Mag    = table(["Expectation";"Area";"Expectation×Area"], ...
                     [24.822; 21.413; 22.325], ...
                     [1.9444e-05; 5.5054e-05; 4.1395e-05], ...
                     'VariableNames', {'Effect','F','p'});

disp('Table 1b — RM-ANOVA (Counts):');     disp(ANOVA_Counts);
disp('Table 1b — RM-ANOVA (Magnitude):');  disp(ANOVA_Mag);

% Optional: write CSVs for manuscript source
% writetable(Pairwise,    'Table1a_Pairwise_AC_PFC.csv');
% writetable(ANOVA_Counts,'Table1b_RMAnova_Counts_AC_PFC.csv');
% writetable(ANOVA_Mag,   'Table1b_RMAnova_Magnitude_AC_PFC.csv');


%% Helper

function plot_region_grid(data40xFxT, idx, f, tms, labels, figTitle, nValid)
    D = data40xFxT(idx,:,:);  % [20 x F x T]

    % --- Restrict to -200:300 ms and 1–50 Hz ---
    tmask = (tms >= -50 & tms <= 300);
    fmask = (f >= 4   & f <= 30);
    tms   = tms(tmask);
    f     = f(fmask);
    D     = D(:,fmask,tmask);   % [20 x F x T subset]

    % Shared color scale across the 20 plots
    vals = D(:); vals = vals(~isnan(vals));
    if isempty(vals), vals = 0; end
    clim = [min(vals), max(vals)];

    % --- Plot grid ---
    figure('Color','w');
    tl = tiledlayout(4,5,'Padding','compact','TileSpacing','compact');
    for ch = 1:20
        nexttile;
        imagesc(tms, f, squeeze(D(ch,:,:)));
        axis xy tight;
        set(gca,'CLim',clim);
        xlabel('Time (ms)'); ylabel('Hz');
        if ~isempty(labels) && numel(labels) >= idx(ch)
            title(labels{idx(ch)}, 'Interpreter','none');
        else
            title(sprintf('Ch %d', idx(ch)));
        end
    end
    colormap jet;  % jet color palette
    cb = colorbar; cb.Layout.Tile = 'east';
    title(tl, sprintf('%s (n=%d sessions)', figTitle, nValid));
end

