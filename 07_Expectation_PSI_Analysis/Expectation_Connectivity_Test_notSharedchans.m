
%% Expectation_Connectivity_Test

% This code tests how connectivity between PFC and AC changes as a function
% of epoch, expectation, and congruence. Specifically, this code tests how the PFC-AC connectivity
% during the testTone epoch changes as a function of expectation. 
% 
% Expectation Condition 1: PriorOnly Trials (informative LED followed by
% silence)
% 
% Expectation Condition 2: PretoneOnly Trials (neutral LED followed by
% pretones).

% This analysis will be carried out session-wise on Correct-Congruent Trials, 
% and on all SNRs.

% Both coherence and PSI will be calculated. 

% In order to circumvent monte carlo estimation, we decided to subsample
% the data, and just deal with the 40-60 percent loss in trials in the
% condition with more trials. 
 
% Statistics

% Because all heuristics are trials balanced, the magnitude of the two
% conditins are directly comparable via nonparametric testing. This also
% allows us to subtract the values between conditions which is useful for
% assessing how the two conditions change across sessions. 

%% Set parameters for the testTone Condition

Epoch             = 'testToneOnset';
datadir_root      = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\04b_Epoc_Cut';
chandir_root      = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05b_Signficant_Channels_epoch_fullArray';
savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\205_07_28_PSI_Test';
Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];   % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 

%% Calculate PSI Spectra during the testTone epoch; PriorOnly Trials, Collapsed Across Behavior (Correct and Wrong) and SNR (High and Low)
% ahead of the development of this code, we found that there was not a
% strong influence of SNR on PSI so we collapse across all levels of SNR

    % Step 1: Get list of modulated channels 

    % Because theta/alpha & beta have the chance to use different channels,
    % you have to make the code Frequency_Band dependent. 

   
        % session_info = {'MrCassius', '190330';
        %                 'MrCassius', '190404';
        %                 'MrCassius', '190413';
        %                 'MrCassius', '190416';
        %                 'MrM', '190417';
        %                 'MrCassius', '190418';
        %                 'MrCassius', '190419';
        %                 'MrCassius', '190421';
        %                 'MrM', '190422';
        %                 'MrCassius', '190423';
        %                 'MrM', '190425';
        %                 'MrM', '190427';
        %                 'MrM', '190502';
        %                 'MrM', '190514';
        %                 'MrCassius', '190515';
        %                 'MrCassius', '190517';
        %                 'MrM', '190525';
        %                 'MrM', '190527';
        %                 'MrM', '190530';
        %                 'MrCassius', '190531';
        %                 'MrM', '190601';
        %                 'MrCassius', '190603';
        %                 'MrM', '190604';
        %                 'MrCassius', '190605';
        %                 'MrCassius', '190703';
        %                 'MrM', '190704';
        %                 'MrM', '190709';
        %                 'MrCassius', '190711';
        %                 'MrCassius', '190713';
        %                 'MrM', '190717';
        %                 'MrM', '190719';
        %                 'MrCassius', '190720';
        %                 'MrM', '190722';
        %                 'MrCassius', '190723';
        %                 'MrCassius', '190725';};

                    session_info = {
                            'MrCassius', '190330';
                            'MrCassius', '190404';
                            'MrCassius', '190413';
                            'MrCassius', '190416';
                            'MrCassius', '190418';
                            'MrCassius', '190419';
                            'MrCassius', '190421';
                            'MrCassius', '190423';
                            'MrCassius', '190515';
                            'MrCassius', '190517';
                            'MrCassius', '190531';
                            'MrCassius', '190603';
                            'MrCassius', '190605';
                            'MrCassius', '190703';
                            'MrCassius', '190711';
                            'MrCassius', '190713';
                            'MrCassius', '190720';
                            'MrCassius', '190723';
                            'MrCassius', '190725';};
 

                    % session_info = {
                    %             'MrM', '190417';
                    %             'MrM', '190422';
                    %             'MrM', '190425';
                    %             'MrM', '190427';
                    %             'MrM', '190502';
                    %             'MrM', '190514';
                    %             'MrM', '190525';
                    %             'MrM', '190527';
                    %             'MrM', '190530';
                    %             'MrM', '190601';
                    %             'MrM', '190604';
                    %             'MrM', '190704';
                    %             'MrM', '190709';
                    %             'MrM', '190717';
                    %             'MrM', '190719';
                    %             'MrM', '190722';};

        for rd = 1:length(session_info(:,1))

                % Establish the subset of data that you want to process. 
                Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
                RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 
        
                % Set directory and get a list of all files in the folder with the desired file name pattern.
                datadir           = fullfile(datadir_root, Animal, Epoch, RecDate);                                 % epoch data, with the time chunk that you want isolated
                chandir           = fullfile(chandir_root, Animal, Epoch, RecDate);     % the sig channels that are output LFP_Spectral_Analysis
                sessions          = dir(fullfile(datadir,'*.mat')); % bad naming convention; only ever one session at a time
                addpath(genpath(datadir));
                addpath(genpath(chandir));
        
                % make save folders and directory
                savedir = fullfile(savedir, RecDate, Animal);
                if ~exist(savedir, 'dir')  % make folders if they don't exist already
                   mkdir(savedir);
                end   
        
        % —— Step 1: Load sig‑channel files and build channel‑pair lists ——
        % — OnlyPrior —

        if strcmp(Epoch,'testToneOnset') == 1
            OnlyPrior_correct_SigChans_fn = sprintf('SigChannels_%s_testToneOnset_OnlyPrior_Correct_%s.mat', RecDate, Frequency_Band);
            OnlyPrior_correct_files       = dir(fullfile(chandir, OnlyPrior_correct_SigChans_fn));
            
            if isempty(OnlyPrior_correct_files)
                fprintf('Skipping %s %s: no OnlyPrior sig‑channel file.\n', Animal, RecDate);
                continue;
            end
        end

        if strcmp(Epoch,'preCueOnset') == 1
            OnlyPrior_correct_SigChans_fn = sprintf('SigChannels_%s_preCueOnset_OnlyPrior_Correct_%s.mat', RecDate, Frequency_Band);
            OnlyPrior_correct_files       = dir(fullfile(chandir, OnlyPrior_correct_SigChans_fn));
            
            if isempty(OnlyPrior_correct_files)
                fprintf('Skipping %s %s: no OnlyPrior sig‑channel file.\n', Animal, RecDate);
                continue;
            end
        end


        % load & split
        OnlyPrior_correct_SigChans          = load(fullfile(chandir, OnlyPrior_correct_SigChans_fn), 'significant_channels');
        OnlyPrior_correct_modifiedCellArray   = regexprep(OnlyPrior_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        OnlyPrior_correct_PFC_chans  =  OnlyPrior_correct_modifiedCellArray(startsWith(OnlyPrior_correct_modifiedCellArray, '*PFC'));
        OnlyPrior_correct_AC_chans   =  OnlyPrior_correct_modifiedCellArray(startsWith(OnlyPrior_correct_modifiedCellArray, '*AC'));


        % — OnlyPretone —

         if strcmp(Epoch,'testToneOnset') == 1
            OnlyPretone_correct_SigChans_fn = sprintf('SigChannels_%s_testToneOnset_OnlyPretone_Correct_%s.mat', RecDate, Frequency_Band);
            OnlyPretone_correct_files       = dir(fullfile(chandir, OnlyPretone_correct_SigChans_fn));
            
            if isempty(OnlyPretone_correct_files)
                fprintf('Skipping %s %s: no OnlyPretone sig‑channel file.\n', Animal, RecDate);
                continue;
            end
         end

        if strcmp(Epoch,'preCueOnset') == 1
            OnlyPretone_correct_SigChans_fn = sprintf('SigChannels_%s_preCueOnset_OnlyPretone_Correct_%s.mat', RecDate, Frequency_Band);
            OnlyPretone_correct_files       = dir(fullfile(chandir, OnlyPretone_correct_SigChans_fn));
            
            if isempty(OnlyPretone_correct_files)
                fprintf('Skipping %s %s: no OnlyPretone sig‑channel file.\n', Animal, RecDate);
                continue;
            end
        end

        % load & split
        OnlyPretone_correct_SigChans = load(fullfile(chandir, OnlyPretone_correct_SigChans_fn), 'significant_channels');
        OnlyPretone_correct_modifiedCellArray   = regexprep(OnlyPretone_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        OnlyPretone_correct_PFC_chans  =  OnlyPretone_correct_modifiedCellArray(startsWith(OnlyPretone_correct_modifiedCellArray, '*PFC'));
        OnlyPretone_correct_AC_chans   =  OnlyPretone_correct_modifiedCellArray(startsWith(OnlyPretone_correct_modifiedCellArray, '*AC'));
        
        % — build channel‑pair lists —
        [OnlyPrior_PFC_comb,   OnlyPrior_AC_comb]     = ndgrid(OnlyPrior_correct_PFC_chans,   OnlyPrior_correct_AC_chans);
        OnlyPrior_ChannelPairs                        = [OnlyPrior_PFC_comb(:),   OnlyPrior_AC_comb(:)];
        
        [OnlyPretone_PFC_comb, OnlyPretone_AC_comb]   = ndgrid(OnlyPretone_correct_PFC_chans, OnlyPretone_correct_AC_chans);
        OnlyPretone_ChannelPairs                      = [OnlyPretone_PFC_comb(:), OnlyPretone_AC_comb(:)];

         % clear variables used to generate shared channel pairs
         clear_list = {'OnlyPrior_AC_comb', 'OnlyPrior_PFC_comb', 'OnlyPretone_AC_comb', 'OnlyPretone_PFC_comb', 'OnlyPrior_correct_files', 'OnlyPretone_correct_files','OnlyPrior_correct_PFC_chans', 'OnlyPretone_correct_PFC_chans', ...
                       'OnlyPrior_correct_AC_chans', 'OnlyPretone_correct_AC_chans', 'OnlyPrior_correct_modifiedCellArray', 'OnlyPretone_correct_modifiedCellArray', ...
                       'OnlyPrior_correct_SigChans_fn', 'OnlyPretone_correct_SigChans_fn', 'OnlyPretone_correct_SigChans_fn', 'OnlyPrior_correct_SigChans', 'OnlyPretone_correct_SigChans'};

         clear(clear_list{:});
         clear clear_list;   

            % ——— after generating OnlyPrior_ChannelPairs & OnlyPretone_ChannelPairs ———
            if isempty(OnlyPrior_ChannelPairs) || isempty(OnlyPretone_ChannelPairs)
                fprintf('Skipping %s %s: no channel pairs to process.\n', Animal, RecDate);
                continue;   % go straight to the next rd in the for‑loop
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
         OnlyPretone_params.prior         = cell2char(OnlyPretone_data.prior); 
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

            % Number of bootstrap iterations
            nIterations = 50;
            
            % Preallocate arrays to store PSI calues 
            % Assume size is [nPairs x nFreqs], e.g., 190 channel pairs by 101 frequencies
            OnlyPrior_nPairs   = size(OnlyPrior_ChannelPairs, 1);
            OnlyPretone_nPairs = size(OnlyPretone_ChannelPairs, 1);
            nFreqs = 101;
           
            PSI_OnlyPrior_all   = zeros(OnlyPrior_nPairs, nFreqs, nIterations);
            PSI_OnlyPretone_all = zeros(OnlyPretone_nPairs, nFreqs, nIterations);
            
            for b = 1:nIterations

                % Determine the number of trials in each condition
                numPriorTrials = numel(OnlyPrior_trials.trial);
                numPretoneTrials = numel(OnlyPretone_trials.trial);
                
                % Find the smaller trial count and which dataset is larger
                minTrials = min(numPriorTrials, numPretoneTrials);
                
                % Randomly subsample the larger dataset to match the smaller one,
                % and reconstruct a viable data array for field trip functions. 
                if numPriorTrials > numPretoneTrials
                    idx = randperm(numPriorTrials, minTrials);
                    OnlyPrior_trials_subset.label      = OnlyPrior_trials.label;
                    OnlyPrior_trials_subset.time       = OnlyPrior_trials.time(1:minTrials);
                    OnlyPrior_trials_subset.trial      = OnlyPrior_trials.trial(idx);
                    OnlyPrior_trials_subset.fsample    = OnlyPrior_trials.fsample;
                    OnlyPrior_trials_subset.sampleinfo = OnlyPrior_trials.sampleinfo(1:minTrials);
                    
                    OnlyPretone_trials_subset          = OnlyPretone_trials;
        
                else
        
                    idx = randperm(numPretoneTrials, minTrials);
                    OnlyPretone_trials_subset.label      = OnlyPretone_trials.label;
                    OnlyPretone_trials_subset.time       = OnlyPretone_trials.time(1:minTrials);
                    OnlyPretone_trials_subset.trial      = OnlyPretone_trials.trial(idx);
                    OnlyPretone_trials_subset.fsample    = OnlyPretone_trials.fsample;
                    OnlyPretone_trials_subset.sampleinfo = OnlyPretone_trials.sampleinfo(1:minTrials); 
        
                    OnlyPrior_trials_subset               = OnlyPrior_trials;
        
                end
        
                % Calculate PSI for both Conditions 
        
                    cfg            = [];
                    cfg.method     = 'mtmfft';
                    cfg.taper      = 'dpss';
                    cfg.output     = 'fourier';
                    cfg.tapsmofrq  = 4; 
                    cfg.pad        = 1;
                    cfg.foilim     = [0 100];
        
                    freq_OnlyPrior      = ft_freqanalysis(cfg, OnlyPrior_trials_subset);
                    freq_OnlyPretone    = ft_freqanalysis(cfg, OnlyPretone_trials_subset);
        
                    cfg                   = [];
                    cfg.method            = 'psi';
                    cfg.bandwidth         = 2;
                    cfg.channelcmb        = OnlyPrior_ChannelPairs; % representing the all shared channel pairs for the 
                    PSI_OnlyPrior         = ft_connectivityanalysis(cfg, freq_OnlyPrior);
                    PSI_OnlyPrior.dof     = length(OnlyPrior_trials_subset.trial);
                    
                    cfg                     = [];
                    cfg.method              = 'psi';
                    cfg.bandwidth           = 2;
                    cfg.channelcmb          = OnlyPretone_ChannelPairs; % representing the all shared channel pairs for the 
                    PSI_OnlyPretone         = ft_connectivityanalysis(cfg, freq_OnlyPretone);
                    PSI_OnlyPretone.dof     = length(OnlyPretone_trials_subset.trial);
                    
                    % Store the PSI values
                    PSI_OnlyPrior_all(:, :, b)   = PSI_OnlyPrior.psispctrm;
                    PSI_OnlyPretone_all(:, :, b) = PSI_OnlyPretone.psispctrm;


            end

                    % average the boostrapped subsample 
                    PSI_OnlyPrior.psispctrm    = mean(PSI_OnlyPrior_all, 3);
                    PSI_OnlyPretone.psispctrm  = mean(PSI_OnlyPretone_all, 3);

                    % Now we need to generate the null distribution. 

                    Only
              
                    % save session value 
                    save_file_name = sprintf('%s_%s_%s_%s_%s_%s_PSI.mat', Animal, RecDate, 'OnlyPrior', 'correct', 'congruent', Frequency_Band);
                    save(fullfile(savedir, save_file_name), 'PSI_OnlyPrior');
                    save_file_name = sprintf('%s_%s_%s_%s_%s_%s_PSI.mat', Animal, RecDate, 'OnlyPretone', 'correct', 'congruent', Frequency_Band);
                    save(fullfile(savedir, save_file_name), 'PSI_OnlyPretone')
             
                    % save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Coh.mat', Animal, RecDate, 'OnlyPrior', 'correct', 'congruent', Frequency_Band);
                    % save(fullfile(savedir, save_file_name), 'Coh_OnlyPrior');
                    % save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Coh.mat', Animal, RecDate, 'OnlyPretone', 'correct', 'congruent', Frequency_Band);
                    % save(fullfile(savedir, save_file_name), 'Coh_OnlyPretone');

                    % clear_list = {'freq_OnlyPrior', 'freq_OnlyPretone', 'cfg','PSI_OnlyPrior', 'PSI_OnlyPretone', 'Coh_OnlyPrior', 'Coh_OnlyPretone'};
                    clear_list = {'freq_OnlyPrior', 'freq_OnlyPretone', 'cfg','PSI_OnlyPrior', 'PSI_OnlyPretone'};

                    clear(clear_list{:});
                    clear clear_list;


end

     