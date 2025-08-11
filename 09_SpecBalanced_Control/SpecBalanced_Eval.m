
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


%% saving top-cut LFP labels

% Initialize tfreq collection for OnlyPrior and OnlyPretone
SpecBalanced_Labels_OnlyPrior_PFC   = cell(20, 1);     % each cell = matrix for one PFC channel
SpecBalanced_Labels_OnlyPretone_PFC = cell(20, 1);

SpecBalanced_Labels_OnlyPrior_AC    = cell(20, 1);      % each cell = matrix for one AC channel
SpecBalanced_Labels_OnlyPretone_AC  = cell(20, 1);




%% Session wise generate the meta data 
       

     for rd = 1:length(session_info(:,1))

         % Initialize outputs for SpecBalancedRangeBased

        SpecBalanced_theta.OnlyPrior.PFC = {};
        SpecBalanced_theta.OnlyPrior.AC  = {};
        SpecBalanced_theta.OnlyPretone.PFC = {};
        SpecBalanced_theta.OnlyPretone.AC  = {};


        % Establish the subset of data that you want to process. 

        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 
        
     
        % Set directory and get a list of all files in the folder with the desired file name pattern.

        datadir           = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\02_Preprocessed', Animal, 'testToneOnset', RecDate);                      % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);             % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\auditory research la\OneDrive\Desktop\Sonia_JunkPile';                                                                    % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
        
        
        sessions = dir(fullfile(datadir,'*.mat'));              % bad naming convention; only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory

        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')                               % make folders if they don't exist already
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
                        OnlyPrior_data.congruency{i} = 'congruent';     % If both are 'H'
                    elseif strcmp(OnlyPrior_data.prior{i}, 'L') && strcmp(OnlyPrior_data.target{i}, 'L')
                        OnlyPrior_data.congruency{i} = 'congruent';     % If both are 'L'
                    elseif strcmp(OnlyPrior_data.prior{i}, 'N')
                       OnlyPrior_data.congruency{i} = 'neutral';        % If prior is 'N'
                    elseif (strcmp(OnlyPrior_data.prior{i}, 'H') && strcmp(OnlyPrior_data.target{i}, 'L')) || (strcmp(OnlyPrior_data.prior{i}, 'L') && strcmp(OnlyPrior_data.target{i}, 'H'))
                       OnlyPrior_data.congruency{i} = 'incongruent';    % If prior and target are different
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
                        OnlyPretone_data.congruency{i} = 'congruent';       % If both are 'H'
                    elseif strcmp(OnlyPretone_data.pretone{i}, 'L') && strcmp(OnlyPretone_data.target{i}, 'L')
                        OnlyPretone_data.congruency{i} = 'congruent';       % If both are 'L'
                    elseif strcmp(OnlyPretone_data.pretone{i}, 'N')
                        OnlyPretone_data.congruency{i} = 'neutral';         % If prior is 'N'
                    elseif (strcmp(OnlyPretone_data.pretone{i}, 'H') && strcmp(OnlyPretone_data.target{i}, 'L')) || (strcmp(OnlyPretone_data.pretone{i}, 'L') && strcmp(OnlyPretone_data.target{i}, 'H'))
                        OnlyPretone_data.congruency{i} = 'incongruent';     % If prior and target are different
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
                % clear SpecBalanced_theta_OnlyPrior SpecBalanced_theta_OnlyPretone

                % SpecBalanced_theta = {};
     

        end


