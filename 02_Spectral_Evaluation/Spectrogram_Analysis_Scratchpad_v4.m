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

        session_info = {'MrCassius', '190330'; % aligned 
                        %'MrCassius', '190404'; % not aligned, excitatory
                        'MrCassius', '190413'; % aligned
                        'MrCassius', '190416'; % aligned
                        'MrM', '190417';       % aligned inhibitory
                        'MrCassius', '190418'; % aligned inhibitory
                        'MrCassius', '190419'; % aligned inhibitory 
                        'MrCassius', '190421'; % aligned inhibitory 
                        'MrM', '190422';       % aligned inhibitory
                        'MrCassius', '190423'; % aligned excitatory
                        'MrM', '190425';       % aligned excitatory 
                        'MrM', '190427';       % aligned excitatory
                        'MrM', '190502';       % aligned excitatory
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


        % if strcmp(Frequency_Band, 'theta') || strcmp(Frequency_Band, 'alpha') % theta and alpha have the same shared pairs according to LFP_spec_eval, so they can use the same set/
        % 
        % session_info = {
        %                 'MrCassius', '190330';
        %                };
        % 
        % end

%% Initiate arrays for meta data collection

% === Containers to retain TF results per-session ===
nSessions = size(session_info, 1);

% Per-session FieldTrip outputs
tfreq_OnlyPrior_high_sessions     = cell(nSessions, 1);
tfreq_OnlyPrior_low_sessions      = cell(nSessions, 1);
tfreq_OnlyPrior_all_sessions      = cell(nSessions, 1);
tfreq_OnlyPretone_high_sessions   = cell(nSessions, 1);
tfreq_OnlyPretone_low_sessions    = cell(nSessions, 1);
tfreq_OnlyPretone_all_sessions    = cell(nSessions, 1);

%% Session wise generate the meta data 
       

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 
        
     
        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir           = fullfile('\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\02_Preprocessed_bipolar', Animal, 'testToneOnset', RecDate);             % epoch data, with the time chunk that you want isolated
        % chandir           = fullfile('\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_bipolar', Animal, 'testTone', RecDate);    % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Desktop\Figure_Spectral_Evaluation';                                                           % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
        
        
        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention; only ever one session at a time
        addpath(genpath(datadir));
        % addpath(genpath(chandir));

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

         % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

         OnlyPrior_data   = ConvertStim(OnlyPrior_data, RecDate);
         OnlyPretone_data = ConvertStim(OnlyPretone_data, RecDate);

         % set parameters required for fieldtrip functions and data selection
         OnlyPrior_params.choice        = OnlyPrior_data.choice;
         OnlyPrior_params.err           = OnlyPrior_data.err;
         OnlyPrior_params.pretone       = OnlyPrior_data.pretone;
         OnlyPrior_params.pretoneLength = OnlyPrior_data.pretoneLength;
         OnlyPrior_params.prior         = cell2char(OnlyPrior_data.prior); 
         OnlyPrior_params.target        = OnlyPrior_data.target;
         OnlyPrior_params.SNR           = OnlyPrior_data.SNR;

         % set parameters required for fieldtrip functions and data selection
         OnlyPretone_params.choice        =  OnlyPretone_data.choice;
         OnlyPretone_params.err           =  OnlyPretone_data.err;
         OnlyPretone_params.pretone       =  OnlyPretone_data.pretone;
         OnlyPretone_params.pretoneLength =  OnlyPretone_data.pretoneLength;
         OnlyPretone_params.prior         =  cell2char(OnlyPretone_data.prior); 
         OnlyPretone_params.target        =  OnlyPretone_data.target;
         OnlyPretone_params.SNR           =  OnlyPretone_data.SNR;

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

            % Select the OnlyPrior high tone trials 
            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = 'c';                                   % Select correct trials
            iSelect.SNR          = SNR;                                   % Select all SNR values
            iSelect.congruency   = [];                                    % Select congruent trials
            iSelect.target       = 'H';                                     % 'H' or 'L' target tone identity 
            OnlyPrior_trials_high    = selectData(OnlyPrior_data.data, OnlyPrior_params, iSelect);     % 

            % Select the OnlyPrior low tone trials 
            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = 'c';                                   % Select correct trials
            iSelect.SNR          = SNR;                                   % Select all SNR values
            iSelect.congruency   = [];                                    % Select congruent trials
            iSelect.target       = 'L';                                     % 'H' or 'L' target tone identity 
            OnlyPrior_trials_low    = selectData(OnlyPrior_data.data, OnlyPrior_params, iSelect);     % 

            % Select all  OnlyPrior trials 
            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = 'c';                                   % Select correct trials
            iSelect.SNR          = SNR;                                   % Select all SNR values
            iSelect.congruency   = [];                                    % Select congruent trials
            iSelect.target       = [];                                     % 'H' or 'L' target tone identity 
            OnlyPrior_trials_all     = selectData(OnlyPrior_data.data, OnlyPrior_params, iSelect);     % 

            % Select the OnlyPretone high tone trials 
            iSelect              = setStimulusCondition('OnlyPretone');
            iSelect.err          = 'c';                                   % Select correct trials
            iSelect.SNR          = SNR;                                   % Select all SNR values
            iSelect.congruency   = [];                                    % Select congruent trials
            iSelect.target       = 'H';                                     % 'H' or 'L' target tone identity 
            OnlyPretone_trials_high    = selectData(OnlyPretone_data.data, OnlyPretone_params, iSelect);     % 

            % Select the OnlyPretone low tone trials 
            iSelect              = setStimulusCondition('OnlyPretone');
            iSelect.err          = 'c';                                   % Select correct trials
            iSelect.SNR          = SNR;                                   % Select all SNR values
            iSelect.congruency   = [];                                    % Select congruent trials
            iSelect.target       = 'L';                                     % 'H' or 'L' target tone identity 
            OnlyPretone_trials_low    = selectData(OnlyPretone_data.data, OnlyPretone_params, iSelect);     % 

            % Select the OnlyPretone trials 
            iSelect              = setStimulusCondition('OnlyPretone');
            iSelect.err          = 'c';                                    % Select  correct trials
            iSelect.SNR          = SNR;                                    % Select all SNR values
            iSelect.congruency   = [];                                     % Select congruent trials
            iSelect.target       = [];                                          
            OnlyPretone_trials_all   = selectData(OnlyPretone_data.data, OnlyPretone_params,iSelect);

            % fourier analysis

            % cfg                 = [];
            % cfg.method          = 'mtmfft';
            % cfg.output          = 'fourier';
            % cfg.taper           = 'dpss';                                                           
            % cfg.tapsmofrq       = 4;      
            % cfg.pad             = 2.0;
            % cfg.foilim          = [1 100];
            % freq_OnlyPrior      = ft_freqanalysis(cfg, OnlyPrior_trials_both);
            % freq_OnlyPretone    = ft_freqanalysis(cfg, OnlyPretone_trials_both);
            % fd_OnlyPrior        = ft_freqdescriptives(cfg, freq_OnlyPrior);
            % fd_OnlyPretone      = ft_freqdescriptives(cfg, freq_OnlyPretone);
            % nTrial_OnlyPrior    = numel(OnlyPrior_trials_both.trial);
            % nTrial_OnlyPretone  = numel(OnlyPretone_trials_both.trial);
            % nTaper_OnlyPrior    = size(freq_OnlyPrior.fourierspctrm,1) / nTrial_OnlyPrior;     % number of tapers
            % nTaper_OnlyPretone  = size(freq_OnlyPretone.fourierspctrm,1) / nTrial_OnlyPretone; % number of tapers

            % Time frequency analysis

            t_length     = 2;                  % make the trial length 2 sec with zero-padding
            t_win        = 0.15;               % 150 ms 
            t_slidingwin = -0.8:0.01:0.6;     % 10-ms sliding window
        
            cfg               = [];
            cfg.method        = 'mtmconvol';
            cfg.output        = 'powandcsd';                          % Options:'powandcsd'; 'fourier'
            cfg.foi           = 1:1/(t_length/2):100;
            cfg.taper         = 'hanning';                            % Options: 'hanning'; 'dpss'
            cfg.t_ftimwin     = ones(1,length(cfg.foi)) .* t_win;     % 250-ms sliding window
            cfg.toi           = t_slidingwin; 
            cfg.keeptrials    = 'no';                                 % Options: 'no'; 'yes'
            cfg.pad           = t_length;
        
            tfreq_OnlyPrior_high     = ft_freqanalysis(cfg, OnlyPrior_trials_high);
            tfreq_OnlyPrior_low      = ft_freqanalysis(cfg, OnlyPrior_trials_low);
            tfreq_OnlyPrior_all      = ft_freqanalysis(cfg, OnlyPrior_trials_all);
            
            tfreq_OnlyPretone_high   = ft_freqanalysis(cfg, OnlyPretone_trials_high);
            tfreq_OnlyPretone_low    = ft_freqanalysis(cfg, OnlyPretone_trials_low);
            tfreq_OnlyPretone_all    = ft_freqanalysis(cfg, OnlyPretone_trials_all);

            % Save per-session TF results
            tfreq_OnlyPrior_high_sessions{rd}   = tfreq_OnlyPrior_high;
            tfreq_OnlyPrior_low_sessions{rd}    = tfreq_OnlyPrior_low;
            tfreq_OnlyPrior_all_sessions{rd}    = tfreq_OnlyPrior_all;
            
            tfreq_OnlyPretone_high_sessions{rd} = tfreq_OnlyPretone_high;
            tfreq_OnlyPretone_low_sessions{rd} = tfreq_OnlyPretone_low;
            tfreq_OnlyPretone_all_sessions{rd} = tfreq_OnlyPretone_all;


      end



%% Grand-average across sessions  collapsing channels (OnlyPrior, ALL trials)
firstIdx_prior = find(~cellfun(@isempty, tfreq_OnlyPrior_all_sessions), 1, 'first');
assert(~isempty(firstIdx_prior), 'No OnlyPrior TF data found.');
ref_prior    = tfreq_OnlyPrior_all_sessions{firstIdx_prior};
sumPow_prior = zeros(size(ref_prior.powspctrm), 'like', ref_prior.powspctrm);
nValid_prior = 0;

for k = 1:numel(tfreq_OnlyPrior_all_sessions)
    if isempty(tfreq_OnlyPrior_all_sessions{k}), continue; end
    A = tfreq_OnlyPrior_all_sessions{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow_prior))
        sumPow_prior = sumPow_prior + A;
        nValid_prior = nValid_prior + 1;
    else
        warning('OnlyPrior (ALL): Skipping session %d (size mismatch)', k);
    end
end
assert(nValid_prior > 0, 'No OnlyPrior sessions matched expected size.');
grandMean_pow_prior = sumPow_prior / nValid_prior;      % [40 x F x T]

f_prior      = ref_prior.freq(:);
tms_prior    = ref_prior.time(:) * 1000;                % ms
labels_prior = ref_prior.label;                         % 40x1 cell


%% Grand-average across sessions  collapsing channels (OnlyPretone, ALL trials)

firstIdx_pret = find(~cellfun(@isempty, tfreq_OnlyPretone_all_sessions), 1, 'first');
assert(~isempty(firstIdx_pret), 'No OnlyPretone TF data found.');
ref    = tfreq_OnlyPretone_all_sessions{firstIdx_pret};
sumPow = zeros(size(ref.powspctrm), 'like', ref.powspctrm);
nValid = 0;

for k = 1:numel(tfreq_OnlyPretone_all_sessions)
    if isempty(tfreq_OnlyPretone_all_sessions{k}), continue; end
    A = tfreq_OnlyPretone_all_sessions{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow))
        sumPow = sumPow + A;
        nValid = nValid + 1;
    else
        warning('OnlyPretone (ALL): Skipping session %d (size mismatch)', k);
    end
end
assert(nValid>0, 'No OnlyPretone sessions matched expected size.');
grandMean_pow_pretone = sumPow / nValid;                % [40 x F x T]

f      = ref.freq(:);
tms    = ref.time(:) * 1000;                            % ms
labels = ref.label;


%% Grand-average across sessions BY TARGET (OnlyPrior: high vs low)

% ---------- High target ----------
firstIdx_prior_high = find(~cellfun(@isempty, tfreq_OnlyPrior_high_sessions), 1, 'first');
assert(~isempty(firstIdx_prior_high), 'No OnlyPrior HIGH TF data found.');
ref_prior_high    = tfreq_OnlyPrior_high_sessions{firstIdx_prior_high};
sumPow_prior_high = zeros(size(ref_prior_high.powspctrm), 'like', ref_prior_high.powspctrm);
nValid_prior_high = 0;

for k = 1:numel(tfreq_OnlyPrior_high_sessions)
    if isempty(tfreq_OnlyPrior_high_sessions{k}), continue; end
    A = tfreq_OnlyPrior_high_sessions{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow_prior_high))
        sumPow_prior_high = sumPow_prior_high + A;
        nValid_prior_high = nValid_prior_high + 1;
    else
        warning('OnlyPrior HIGH: Skipping session %d (size mismatch)', k);
    end
end
assert(nValid_prior_high > 0, 'No OnlyPrior HIGH sessions matched expected size.');
grandMean_pow_prior_high = sumPow_prior_high / nValid_prior_high;  % [40 x F x T]

% ---------- Low target ----------
firstIdx_prior_low = find(~cellfun(@isempty, tfreq_OnlyPrior_low_sessions), 1, 'first');
assert(~isempty(firstIdx_prior_low), 'No OnlyPrior LOW TF data found.');
ref_prior_low    = tfreq_OnlyPrior_low_sessions{firstIdx_prior_low};
sumPow_prior_low = zeros(size(ref_prior_low.powspctrm), 'like', ref_prior_low.powspctrm);
nValid_prior_low = 0;

for k = 1:numel(tfreq_OnlyPrior_low_sessions)
    if isempty(tfreq_OnlyPrior_low_sessions{k}), continue; end
    A = tfreq_OnlyPrior_low_sessions{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow_prior_low))
        sumPow_prior_low = sumPow_prior_low + A;
        nValid_prior_low = nValid_prior_low + 1;
    else
        warning('OnlyPrior LOW: Skipping session %d (size mismatch)', k);
    end
end
assert(nValid_prior_low > 0, 'No OnlyPrior LOW sessions matched expected size.');
grandMean_pow_prior_low = sumPow_prior_low / nValid_prior_low;    % [40 x F x T]


%% Grand-average across sessions BY TARGET (OnlyPretone: high vs low)

% ---------- High target ----------
firstIdx_pret_high = find(~cellfun(@isempty, tfreq_OnlyPretone_high_sessions), 1, 'first');
assert(~isempty(firstIdx_pret_high), 'No OnlyPretone HIGH TF data found.');
ref_pret_high    = tfreq_OnlyPretone_high_sessions{firstIdx_pret_high};
sumPow_pret_high = zeros(size(ref_pret_high.powspctrm), 'like', ref_pret_high.powspctrm);
nValid_pret_high = 0;

for k = 1:numel(tfreq_OnlyPretone_high_sessions)
    if isempty(tfreq_OnlyPretone_high_sessions{k}), continue; end
    A = tfreq_OnlyPretone_high_sessions{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow_pret_high))
        sumPow_pret_high = sumPow_pret_high + A;
        nValid_pret_high = nValid_pret_high + 1;
    else
        warning('OnlyPretone HIGH: Skipping session %d (size mismatch)', k);
    end
end
assert(nValid_pret_high > 0, 'No OnlyPretone HIGH sessions matched expected size.');
grandMean_pow_pretone_high = sumPow_pret_high / nValid_pret_high; % [40 x F x T]

% ---------- Low target ----------
firstIdx_pret_low = find(~cellfun(@isempty, tfreq_OnlyPretone_low_sessions), 1, 'first');
assert(~isempty(firstIdx_pret_low), 'No OnlyPretone LOW TF data found.');
ref_pret_low    = tfreq_OnlyPretone_low_sessions{firstIdx_pret_low};
sumPow_pret_low = zeros(size(ref_pret_low.powspctrm), 'like', ref_pret_low.powspctrm);
nValid_pret_low = 0;

for k = 1:numel(tfreq_OnlyPretone_low_sessions)
    if isempty(tfreq_OnlyPretone_low_sessions{k}), continue; end
    A = tfreq_OnlyPretone_low_sessions{k}.powspctrm;   % [40 x F x T]
    if isequal(size(A), size(sumPow_pret_low))
        sumPow_pret_low = sumPow_pret_low + A;
        nValid_pret_low = nValid_pret_low + 1;
    else
        warning('OnlyPretone LOW: Skipping session %d (size mismatch)', k);
    end
end
assert(nValid_pret_low > 0, 'No OnlyPretone LOW sessions matched expected size.');
grandMean_pow_pretone_low = sumPow_pret_low / nValid_pret_low;    % [40 x F x T]


%% Collapse across frequency to get band-averaged power over time

% Map the Frequency_Band string to a numeric range
switch lower(Frequency_Band)
    case 'theta'
        fwin = [4 8];
    case 'alpha'
        fwin = [9 14];
    case 'beta'
        fwin = [15 30];
    otherwise
        warning('Unknown Frequency_Band "%s", defaulting to theta (4–8 Hz).', Frequency_Band);
        fwin = [4 8];
end

% ---------- ALL TRIALS: OnlyPrior ----------
fmask_prior_all = (f_prior >= fwin(1) & f_prior <= fwin(2));

% band-averaged power per channel over time: [40 x T]
band_prior_ch = squeeze(mean(grandMean_pow_prior(:, fmask_prior_all, :), 2));  

% region averages (PFC = 1:20, AC = 21:40)
PFC_prior = mean(band_prior_ch(1:20, :), 1);   % [1 x T]
AC_prior  = mean(band_prior_ch(21:40, :), 1);  % [1 x T]

% ---------- ALL TRIALS: OnlyPretone ----------
fmask_pret_all = (f >= fwin(1) & f <= fwin(2));

band_pret_ch = squeeze(mean(grandMean_pow_pretone(:, fmask_pret_all, :), 2));  % [40 x T]

PFC_pret = mean(band_pret_ch(1:20, :), 1);    % [1 x T]
AC_pret  = mean(band_pret_ch(21:40, :), 1);   % [1 x T]


%% Collapse across frequency for HIGH vs LOW targets (4–8 Hz band-averaged power over time)

% ---------- OnlyPrior HIGH vs LOW ----------
fmask_prior = (f_prior >= fwin(1) & f_prior <= fwin(2));

band_prior_ch_high = squeeze(mean(grandMean_pow_prior_high(:, fmask_prior, :), 2));  % [40 x T]
band_prior_ch_low  = squeeze(mean(grandMean_pow_prior_low(:,  fmask_prior, :), 2));  % [40 x T]

PFC_prior_high = mean(band_prior_ch_high(1:20, :), 1);   % [1 x T]
AC_prior_high  = mean(band_prior_ch_high(21:40, :), 1);  % [1 x T]

PFC_prior_low  = mean(band_prior_ch_low(1:20, :), 1);    % [1 x T]
AC_prior_low   = mean(band_prior_ch_low(21:40, :), 1);   % [1 x T]

% ---------- OnlyPretone HIGH vs LOW ----------
fmask_pret = (f >= fwin(1) & f <= fwin(2));

band_pret_ch_high = squeeze(mean(grandMean_pow_pretone_high(:, fmask_pret, :), 2)); % [40 x T]
band_pret_ch_low  = squeeze(mean(grandMean_pow_pretone_low(:,  fmask_pret, :), 2)); % [40 x T]

PFC_pret_high = mean(band_pret_ch_high(1:20, :), 1);   % [1 x T]
AC_pret_high  = mean(band_pret_ch_high(21:40, :), 1);  % [1 x T]

PFC_pret_low  = mean(band_pret_ch_low(1:20, :), 1);    % [1 x T]
AC_pret_low   = mean(band_pret_ch_low(21:40, :), 1);   % [1 x T]


%% Compute theta response index per channel (stim vs prestim, ALL trials)

% Time windows (ms)
prestim_win = [-100 0];   % baseline
stim_win    = [0 250];    % response

% masks on tms_prior / tms (in ms)
preMask_prior  = (tms_prior >= prestim_win(1) & tms_prior <  prestim_win(2));
stimMask_prior = (tms_prior >= stim_win(1)    & tms_prior <= stim_win(2));

preMask_pret   = (tms       >= prestim_win(1) & tms       <  prestim_win(2));
stimMask_pret  = (tms       >= stim_win(1)    & tms       <= stim_win(2));

% mean theta power per channel in each window: [40 x 1]
pre_prior  = mean(band_prior_ch(:, preMask_prior),  2);
stim_prior = mean(band_prior_ch(:, stimMask_prior), 2);

pre_pret   = mean(band_pret_ch(:,  preMask_pret),   2);
stim_pret  = mean(band_pret_ch(:,  stimMask_pret),  2);

% Response index (modulation index): (stim - pre) / (stim + pre)
epsVal = 1e-12;  % avoid divide-by-zero

respIdx_prior = (stim_prior - pre_prior) ./ (stim_prior + pre_prior + epsVal);
respIdx_pret  = (stim_pret  - pre_pret)  ./ (stim_pret  + pre_pret  + epsVal);


%% Tone selectivity index (baseline-corrected high vs low) per channel

% ---------- OnlyPrior ----------
% baseline & stim for HIGH
pre_prior_high  = mean(band_prior_ch_high(:, preMask_prior),  2);   % [40 x 1]
stim_prior_high = mean(band_prior_ch_high(:, stimMask_prior), 2);   % [40 x 1]

% baseline & stim for LOW
pre_prior_low   = mean(band_prior_ch_low(:,  preMask_prior),  2);   % [40 x 1]
stim_prior_low  = mean(band_prior_ch_low(:,  stimMask_prior), 2);   % [40 x 1]

% baseline-corrected responses
resp_prior_high = stim_prior_high - pre_prior_high;   % [40 x 1]
resp_prior_low  = stim_prior_low  - pre_prior_low;    % [40 x 1]

% selectivity index: (High - Low) / (High + Low)
selIdx_prior = (resp_prior_high - resp_prior_low) ./ ...
               (resp_prior_high + resp_prior_low + epsVal);   % [40 x 1]

% ---------- OnlyPretone ----------
pre_pret_high  = mean(band_pret_ch_high(:, preMask_pret),  2);  % [40 x 1]
stim_pret_high = mean(band_pret_ch_high(:, stimMask_pret), 2);  % [40 x 1]

pre_pret_low   = mean(band_pret_ch_low(:,  preMask_pret),  2);  % [40 x 1]
stim_pret_low  = mean(band_pret_ch_low(:,  stimMask_pret), 2);  % [40 x 1]

resp_pret_high = stim_pret_high - pre_pret_high;   % [40 x 1]
resp_pret_low  = stim_pret_low  - pre_pret_low;    % [40 x 1]

selIdx_pret = (resp_pret_high - resp_pret_low) ./ ...
              (resp_pret_high + resp_pret_low + epsVal);       % [40 x 1]


%% Scatter plot: Relationship between expectations (response index)
% X: OnlyPretone response index (0–250 vs -100–0)
% Y: OnlyPrior   response index (0–250 vs -100–0)
% vlPFC = channels 1–20 (triangles, filled)
% AC   = channels 21–40 (squares, filled)

figure('Color','w'); hold on;

% vlPFC (1:20) – filled triangles
plot(respIdx_pret(1:20), respIdx_prior(1:20), '^', ...
    'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerFaceColor','auto');

% AC (21:40) – filled squares
plot(respIdx_pret(21:40), respIdx_prior(21:40), 's', ...
    'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerFaceColor','auto');

% ----- Data-driven axis limits with small padding -----
allX = respIdx_pret(:);
allY = respIdx_prior(:);

minX = min(allX);  maxX = max(allX);
minY = min(allY);  maxY = max(allY);

padX = 0.05 * max(maxX - minX, 1e-3);
padY = 0.05 * max(maxY - minY, 1e-3);

xlims = [minX - padX, maxX + padX];
ylims = [minY - padY, maxY + padY];

xlim(xlims);
ylim(ylims);

% Unity line (OnlyPrior = OnlyPretone)
lineMin = min(xlims(1), ylims(1));
lineMax = max(xlims(2), ylims(2));
plot([lineMin lineMax], [lineMin lineMax], 'k--', 'LineWidth', 1);

axis square;

xlabel('Theta response index: OnlyPretone (0–250 ms vs -100–0 ms)', ...
    'FontName','Arial','FontSize',18);
ylabel('Theta response index: OnlyPrior (0–250 ms vs -100–0 ms)', ...
    'FontName','Arial','FontSize',18);
title('Channel-wise theta response: OnlyPretone vs OnlyPrior', ...
    'FontName','Arial','FontSize',18);

leg = legend({'vlPFC (ch 1–20)', 'AC (ch 21–40)'}, 'Location', 'best');
set(leg, 'FontName','Arial', 'FontSize',18);

grid off;
set(gca, 'FontName','Arial', 'FontSize',18);


%% Line plots: High vs Low targets, per region and expectation

% OnlyPrior · PFC
figure('Color','w'); hold on;
p1 = plot(tms_prior, PFC_prior_high, 'LineWidth', 2);
p2 = plot(tms_prior, PFC_prior_low,  'LineWidth', 2);

xlabel('Time (ms)', 'FontName','Arial','FontSize',18);
ylabel('4–8 Hz power (a.u.)', 'FontName','Arial','FontSize',18);
title('OnlyPrior · PFC (ch 1–20): High vs Low targets', ...
      'FontName','Arial','FontSize',10);
xlim([-100 250]);
grid off;

yl = ylim;
plot([0 0], [yl(1) yl(2)], 'k--', 'LineWidth', 1.5);
ylim(yl);

legend([p1 p2], {'High target', 'Low target'}, ...
       'FontName','Arial','FontSize',10, 'Location','best');

set(gca, 'FontName','Arial','FontSize',18);


% OnlyPrior · AC
figure('Color','w'); hold on;
p1 = plot(tms_prior, AC_prior_high, 'LineWidth', 2);
p2 = plot(tms_prior, AC_prior_low,  'LineWidth', 2);

xlabel('Time (ms)', 'FontName','Arial','FontSize',18);
ylabel('4–8 Hz power (a.u.)', 'FontName','Arial','FontSize',18);
title('OnlyPrior · AC (ch 21–40): High vs Low targets', ...
      'FontName','Arial','FontSize',10);
xlim([-100 250]);
grid off;

yl = ylim;
plot([0 0], [yl(1) yl(2)], 'k--', 'LineWidth', 1.5);
ylim(yl);

% legend([p1 p2], {'High target', 'Low target'}, ...
%        'FontName','Arial','FontSize',18, 'Location','best');

set(gca, 'FontName','Arial','FontSize',18);


% OnlyPretone · PFC
figure('Color','w'); hold on;
p1 = plot(tms, PFC_pret_high, 'LineWidth', 2);
p2 = plot(tms, PFC_pret_low,  'LineWidth', 2);

xlabel('Time (ms)', 'FontName','Arial','FontSize',18);
ylabel('4–8 Hz power (a.u.)', 'FontName','Arial','FontSize',18);
title('OnlyPretone · PFC (ch 1–20): High vs Low targets', ...
      'FontName','Arial','FontSize',10);
xlim([-100 250]);
grid off;

yl = ylim;
plot([0 0], [yl(1) yl(2)], 'k--', 'LineWidth', 1.5);
ylim(yl);

% legend([p1 p2], {'High target', 'Low target'}, ...
%        'FontName','Arial','FontSize',18, 'Location','best');

set(gca, 'FontName','Arial','FontSize',18);


% OnlyPretone · AC
figure('Color','w'); hold on;
p1 = plot(tms, AC_pret_high, 'LineWidth', 2);
p2 = plot(tms, AC_pret_low,  'LineWidth', 2);

xlabel('Time (ms)', 'FontName','Arial','FontSize',18);
ylabel('4–8 Hz power (a.u.)', 'FontName','Arial','FontSize',18);
title('OnlyPretone · AC (ch 21–40): High vs Low targets', ...
      'FontName','Arial','FontSize',10);
xlim([-100 250]);
grid off;

yl = ylim;
plot([0 0], [yl(1) yl(2)], 'k--', 'LineWidth', 1.5);
ylim(yl);

% legend([p1 p2], {'High target', 'Low target'}, ...
%        'FontName','Arial','FontSize',18, 'Location','best');

set(gca, 'FontName','Arial','FontSize',18);


%% Scatter plot: Tone selectivity index (OnlyPretone vs OnlyPrior)
% X: OnlyPretone selectivity (High vs Low)
% Y: OnlyPrior   selectivity (High vs Low)

figure('Color','w'); hold on;

% vlPFC (1:20)
plot(selIdx_pret(1:20), selIdx_prior(1:20), '^', ...
    'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerFaceColor','auto');

% AC (21:40)
plot(selIdx_pret(21:40), selIdx_prior(21:40), 's', ...
    'MarkerSize', 12, 'LineWidth', 1.5, ...
    'MarkerFaceColor','auto');

allX = selIdx_pret(:);
allY = selIdx_prior(:);

minX = min(allX);  maxX = max(allX);
minY = min(allY);  maxY = max(allY);

padX = 0.05 * max(maxX - minX, 1e-3);
padY = 0.05 * max(maxY - minY, 1e-3);

xlims = [minX - padX, maxX + padX];
ylims = [minY - padY, maxY + padY];

xlim(xlims);
ylim(ylims);

lineMin = min(xlims(1), ylims(1));
lineMax = max(xlims(2), ylims(2));
plot([lineMin lineMax], [lineMin lineMax], 'k--', 'LineWidth', 1);

axis square;

xlabel('Theta selectivity index: OnlyPretone (High vs Low targets)', ...
    'FontName','Arial','FontSize',18);
ylabel('Theta selectivity index: OnlyPrior (High vs Low targets)', ...
    'FontName','Arial','FontSize',18);
title('Channel-wise theta selectivity: OnlyPretone vs OnlyPrior', ...
    'FontName','Arial','FontSize',18);

leg = legend({'vlPFC (ch 1–20)', 'AC (ch 21–40)'}, 'Location', 'best');
set(leg, 'FontName','Arial','FontSize',18);

grid off;
set(gca, 'FontName','Arial','FontSize',18);


%% Heatmaps: OnlyPrior (4–8 Hz only, ALL trials)

plot_region_grid_4to8(grandMean_pow_prior, 1:20,  f_prior, tms_prior, labels_prior, ...
    'OnlyPrior · PFC (ch 1–20)', nValid_prior);

plot_region_grid_4to8(grandMean_pow_prior, 21:40, f_prior, tms_prior, labels_prior, ...
    'OnlyPrior · AC  (ch 21–40)', nValid_prior);


%% Heatmaps: OnlyPretone (4–8 Hz only, ALL trials)

plot_region_grid_4to8(grandMean_pow_pretone, 1:20,  f, tms, labels, ...
    'OnlyPretone · PFC (ch 1–20)', nValid);

plot_region_grid_4to8(grandMean_pow_pretone, 21:40, f, tms, labels, ...
    'OnlyPretone · AC  (ch 21–40)', nValid);


%% Helper function

function plot_region_grid_4to8(data40xFxT, idx, f, tms, labels, figTitle, nValid)
    % data40xFxT : [40 x F x T]
    % idx        : indices for this region (e.g., 1:20 or 21:40)
    % f, tms     : full freq/time vectors for this condition

    D = data40xFxT(idx,:,:);  % [20 x F x T]

    % --- Restrict to -100:250 ms and 4–8 Hz ---
    tmask = (tms >= -100 & tms <= 250);
    fmask = (f   >= 4    & f   <= 8);
    tms   = tms(tmask);
    f     = f(fmask);
    D     = D(:, fmask, tmask);   % [20 x Fsubset x Tsubset]

    % Local color scale for this region/condition
    vals = D(:);
    vals = vals(~isnan(vals));
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

        hold on;
        plot([0 0], [f(1) f(end)], 'k--', 'LineWidth', 1);
        hold off;
    end

    colormap jet;
    cb = colorbar;
    cb.Layout.Tile = 'east';
    title(tl, sprintf('%s (4–8 Hz, n=%d sessions)', figTitle, nValid), ...
        'FontName','Arial','FontSize',18);
end




