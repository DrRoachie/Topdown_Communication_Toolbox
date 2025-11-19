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

    %     session_info = {'MrCassius', '190330'; % aligned 
    %                     'MrCassius', '190404'; % not aligned
    %                     'MrCassius', '190413'; % aligned
    %                     'MrCassius', '190416'; % aligned
    %                     'MrM', '190417';       % aligned inhibitory
    %                     'MrCassius', '190418'; % aligned inhibitory
    %                     'MrCassius', '190419'; % aligned inhibitory 
    %                     'MrCassius', '190421'; % aligned inhibitory 
    %                     'MrM', '190422';       % aligned inhibitory
    %                     'MrCassius', '190423'; % aligned excitatory
    %                     'MrM', '190425';       % aligned excitatory 
    %                     'MrM', '190427';       % aligned excitatory
    %                     'MrM', '190502';       % aligned excitatory
    %                     'MrM', '190514';
    %                     'MrCassius', '190515';
    %                     'MrCassius', '190517';
    %                     'MrM', '190525';
    %                     'MrM', '190527';
    %                     'MrM', '190530';
    %                     'MrCassius', '190531';
    %                     'MrM', '190601';
    %                     'MrCassius', '190603';
    %                     'MrM', '190604';
    %                     'MrCassius', '190605';
    %                     'MrCassius', '190703';
    %                     'MrM', '190704';
    %                     'MrM', '190709';
    %                     'MrCassius', '190711';
    %                     'MrCassius', '190713';
    %                     'MrM', '190717';
    %                     'MrM', '190719';
    %                     'MrCassius', '190720';
    %                     'MrM', '190722';
    %                     'MrCassius', '190723';
    %                     'MrCassius', '190725';};
    % 
    % end


        session_info = {
                        'MrCassius', '190404';
                        
                       };

    end
%% Initiate arrays for meta data collection

% === Containers to retain TF results per-session ===
nSessions = size(session_info, 1);

% Per-session FieldTrip outputs
tfreq_OnlyPrior_all   = cell(nSessions, 1);
tfreq_OnlyPretone_all = cell(nSessions, 1);

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
            iSelect.congruency   = [];                                     % Select congruent trials
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
        
            tfreq_OnlyPrior     = ft_freqanalysis(cfg, OnlyPrior_trials);
            tfreq_OnlyPretone   = ft_freqanalysis(cfg, OnlyPretone_trials);

            % Save per-session TF results
            tfreq_OnlyPrior_all{rd}   = tfreq_OnlyPrior;
            tfreq_OnlyPretone_all{rd} = tfreq_OnlyPretone;


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

%% Compute global color limits across OnlyPrior and OnlyPretone (theta, -250 to 250 ms)

% Define the window you are actually plotting
twin = [-250 250];   % ms
fwin = [4 8];        % Hz (theta band)

% Masks for OnlyPrior
tmask_prior = (tms_prior >= twin(1) & tms_prior <= twin(2));
fmask_prior = (f_prior  >= fwin(1)  & f_prior  <= fwin(2));

% Masks for OnlyPretone
tmask_pret  = (tms       >= twin(1) & tms      <= twin(2));
fmask_pret  = (f         >= fwin(1) & f        <= fwin(2));

% Subset the data to the plotting window
D_prior   = grandMean_pow_prior(:, fmask_prior, tmask_prior);     % [40 x F x T]
D_pretone = grandMean_pow_pretone(:, fmask_pret, tmask_pret);     % [40 x F x T]

% Flatten and concatenate both conditions
vals_prior   = D_prior(:);
vals_pretone = D_pretone(:);

% Remove NaNs
vals_all = [vals_prior(:); vals_pretone(:)];
vals_all = vals_all(~isnan(vals_all));
if isempty(vals_all)
    vals_all = 0;   % fallback
end

% Global color limits
clim_global = [min(vals_all), max(vals_all)];

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

% ---------- OnlyPrior ----------
% band mask
fmask_prior = (f_prior >= fwin(1) & f_prior <= fwin(2));

% band-averaged power per channel over time: [40 x T]
band_prior_ch = squeeze(mean(grandMean_pow_prior(:, fmask_prior, :), 2));  

% region averages (PFC = 1:20, AC = 21:40)
PFC_prior = mean(band_prior_ch(1:20, :), 1);   % [1 x T]
AC_prior  = mean(band_prior_ch(21:40, :), 1);  % [1 x T]

% ---------- OnlyPretone ----------
fmask_pret = (f >= fwin(1) & f <= fwin(2));

band_pret_ch = squeeze(mean(grandMean_pow_pretone(:, fmask_pret, :), 2));  % [40 x T]

PFC_pret = mean(band_pret_ch(1:20, :), 1);    % [1 x T]
AC_pret  = mean(band_pret_ch(21:40, :), 1);   % [1 x T]


%% Plot 20 per-channel spectrograms for PFC (1–20) and AC (21–40) — OnlyPrior

plot_region_grid(grandMean_pow_prior, 1:20,  f_prior, tms_prior, labels_prior, ...
    'OnlyPrior · PFC (ch 1–20)', nValid_prior, clim_global);

plot_region_grid(grandMean_pow_prior, 21:40, f_prior, tms_prior, labels_prior, ...
    'OnlyPrior · AC  (ch 21–40)', nValid_prior, clim_global);

%% Plot 20 per-channel spectrograms for PFC (1–20) and AC (21–40) - OnlyPretone

plot_region_grid(grandMean_pow_pretone, 1:20,  f, tms, labels, ...
    'OnlyPretone · PFC (ch 1–20)', nValid, clim_global);

plot_region_grid(grandMean_pow_pretone, 21:40, f, tms, labels, ...
    'OnlyPretone · AC  (ch 21–40)', nValid, clim_global);

%% Plot OnlyPrior: PFC vs AC band-averaged power over time

figure('Color','w'); hold on;
plot(tms_prior, PFC_prior, 'LineWidth', 2);
plot(tms_prior, AC_prior,  'LineWidth', 2);
xlabel('Time (ms)');
ylabel(sprintf('%s-band power (a.u.)', Frequency_Band));
title(sprintf('OnlyPrior · mean %s-band power (PFC vs AC)', Frequency_Band));
legend({'PFC','AC'}, 'Location','best');
xlim([-250 250]);  % optional, adjust to taste
grid on;

%% Plot OnlyPretone: PFC vs AC band-averaged power over time

figure('Color','w'); hold on;
plot(tms, PFC_pret, 'LineWidth', 2);
plot(tms, AC_pret,  'LineWidth', 2);
xlabel('Time (ms)');
ylabel(sprintf('%s-band power (a.u.)', Frequency_Band));
title(sprintf('OnlyPretone · mean %s-band power (PFC vs AC)', Frequency_Band));
legend({'PFC','AC'}, 'Location','best');
xlim([-250 250]);  % optional
grid on;


%% Helper

function plot_region_grid(data40xFxT, idx, f, tms, labels, figTitle, nValid, clim)
    % data40xFxT : [40 x F x T]
    % idx        : indices for this region (e.g., 1:20 or 21:40)
    % f, tms     : full freq/time vectors for this condition
    % clim       : 1x2 global [cmin cmax] (precomputed)

    D = data40xFxT(idx,:,:);  % [20 x F x T]

    % --- Restrict to -250:250 ms and 4–8 Hz (theta) ---
    tmask = (tms >= -800 & tms <= 800);
    fmask = (f   >= 4   & f   <= 8);
    tms   = tms(tmask);
    f     = f(fmask);
    D     = D(:, fmask, tmask);   % [20 x Fsubset x Tsubset]

    % If (for some reason) clim wasn’t provided, fall back to local scaling
    if nargin < 8 || isempty(clim)
        vals = D(:); 
        vals = vals(~isnan(vals));
        if isempty(vals), vals = 0; end
        clim = [min(vals), max(vals)];
    end

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

    colormap jet;
    cb = colorbar;
    cb.Layout.Tile = 'east';
    title(tl, sprintf('%s (n=%d sessions)', figTitle, nValid));
end

