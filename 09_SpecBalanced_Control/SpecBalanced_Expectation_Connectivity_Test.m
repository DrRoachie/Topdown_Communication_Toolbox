
%% Expectation_Connectivity_Test

% This code tests how connectivity between PFC and AC changes as a function
% of expectation. Specifically, this code tests how the PFC-AC connectivity
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

%% SpecBalanced Version

% This version of the code is outfitted for Sonia's first mission in which
% she is rerunning all analyses on 'spec-balanced' channels. The pipeline
% is basically the same as the first version of the analysis except rather
% than pulling sig channels from the phase shuffle analysis (those
% modulated over noise), she is running the analysis on channels that share
% spectral overlap between PFC and AC. This is a control that will enable
% us to visualize how matched magnitude influences connectivity metric.
% Simply, do channels that are operating within the same voltage range
% communicate differently than when considering all channels. Following the
% successful implimentation of the code for the theta band, we will revise
% the code for alpha and beta, and also repeat the anlysis with a double
% filter (first filtering for phase shuffle, THEN filtering for overlap). 

% To do this, she must extract sig channels from her output files titled
% 'SpecBalanced_theta'. and use them to assemble cfg.label. 

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


    % session_info = {'MrCassius', '190330';
    %                     'MrCassius', '190404';
    %                     'MrCassius', '190413';
    %                     'MrCassius', '190416';
    %                     'MrM', '190417';
    %                     'MrCassius', '190418';};
    % 

     end


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
    % 
    % 
    % 
    %  end

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 

        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir           = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\04_Epoc_Cut', Animal, 'testTone', RecDate);                       % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);     % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\auditory research la\OneDrive\Desktop\Connectivity_Results';

        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention; only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')  % make folders if they don't exist already
           mkdir(savedir);
        end   

        
        % === Load SpecBalanced channel pairs (instead of SigChannels) ===
        specbalanced_dir = fullfile('D:\SpecBalanced_Eval', RecDate, Animal); % SoniaJunkPile --> RecDate --> Animal
        specbalanced_file = fullfile(specbalanced_dir, sprintf('SpecBalancedTheta_%s_%s.mat', RecDate, Animal));       % Builds the full filename       
        
        if ~exist(specbalanced_file, 'file')
            warning('SpecBalanced file NOT FOUND for %s %s — skipping.\n', Animal, RecDate);                          % if file is missing create empty NULL file 
            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_data.txt', RecDate, Animal, Frequency_Band));
            fid4 = fopen(save_file_name, 'w');
            fclose(fid4);
            continue;
        end
        
        % === Load the file and  extract the variable ===
        loaded = load(specbalanced_file);          %  if file exists, load into a struct
        varname = fieldnames(loaded);              %  get the name of the variable inside the file
        SpecBalanced_Theta = loaded.(varname{1});  %  assign variable to SpecBalanced_Theta
        
        % === Separate OnlyPrior and OnlyPretone channel pair lists, with ZERO-PADDING ===
        % Build OnlyPrior pairwise combinations
        
        prior_PFC_chans = SpecBalanced_Theta.OnlyPrior.PFC;
        prior_AC_chans  = SpecBalanced_Theta.OnlyPrior.AC;
        
        % Clean Prior labels (zero-pad chXX)
        prior_PFC_chans_clean = cellfun(@(x) regexprep(x, 'ch(\d+)$', '${sprintf(''ch%02d'', str2double($1))}'), prior_PFC_chans, 'UniformOutput', false);
        prior_AC_chans_clean  = cellfun(@(x) regexprep(x, 'ch(\d+)$', '${sprintf(''ch%02d'', str2double($1))}'), prior_AC_chans,  'UniformOutput', false);
        
        % Build grid of pairs
        [PFC_prior_grid, AC_prior_grid] = ndgrid(prior_PFC_chans_clean, prior_AC_chans_clean);
        OnlyPrior_Balanced_ChannelPairs = [PFC_prior_grid(:), AC_prior_grid(:)];
        
        % Build OnlyPretone pairwise combinations
        
        pretone_PFC_chans = SpecBalanced_Theta.OnlyPretone.PFC;
        pretone_AC_chans  = SpecBalanced_Theta.OnlyPretone.AC;
        
        % Clean Pretone labels (zero-pad chXX)
        pretone_PFC_chans_clean = cellfun(@(x) regexprep(x, 'ch(\d+)$', '${sprintf(''ch%02d'', str2double($1))}'), pretone_PFC_chans, 'UniformOutput', false);
        pretone_AC_chans_clean  = cellfun(@(x) regexprep(x, 'ch(\d+)$', '${sprintf(''ch%02d'', str2double($1))}'), pretone_AC_chans,  'UniformOutput', false);

        
        % Build grid of pairs
        [PFC_pretone_grid, AC_pretone_grid] = ndgrid(pretone_PFC_chans_clean, pretone_AC_chans_clean);
        OnlyPretone_Balanced_ChannelPairs = [PFC_pretone_grid(:), AC_pretone_grid(:)];
        
        % === Debug print ===
        fprintf('\n--- Session %s %s ---\n', Animal, RecDate);
        fprintf('[OnlyPrior] PFC: %s\n', strjoin(prior_PFC_chans_clean, ', '));
        fprintf('[OnlyPrior] AC: %s\n', strjoin(prior_AC_chans_clean, ', '));
        fprintf('[OnlyPrior] # of pairs: %d\n', size(OnlyPrior_Balanced_ChannelPairs, 1));
        fprintf('[OnlyPretone] PFC: %s\n', strjoin(pretone_PFC_chans_clean, ', '));
        fprintf('[OnlyPretone] AC: %s\n', strjoin(pretone_AC_chans_clean, ', '));
        fprintf('[OnlyPretone] # of pairs: %d\n\n', size(OnlyPretone_Balanced_ChannelPairs, 1));
        

        % extract session information
                 baseFileName = sessions.name;
                 fullFileName = fullfile(sessions.folder, baseFileName);
                 fprintf(1, 'Now reading %s/n', fullFileName);
                 OnlyPrior_data   = load(fullfile(fullFileName));
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

                % --- START OF BOOTSTRAP / CONNECTIVITY TEST ---

            nIterations = 5;
            nFreqs = 101;
            
            % Preallocate result arrays with fixed dimensions
            nPairs_Prior   = size(OnlyPrior_Balanced_ChannelPairs, 1);
            nPairs_Pretone = size(OnlyPretone_Balanced_ChannelPairs, 1);
            
            % % Preallocate numeric arrays
            % PSI_OnlyPrior_all   = NaN(nPairs_Prior, nFreqs, nIterations);    % [pair index x frequency x iteration] 
            % PSI_OnlyPretone_all = NaN(nPairs_Pretone, nFreqs, nIterations);  % use NaN initally so if some pairs are missing we can omit NaN at the end 
            % Coh_OnlyPrior_all   = NaN(nPairs_Prior, nFreqs, nIterations);
            % Coh_OnlyPretone_all = NaN(nPairs_Pretone, nFreqs, nIterations);
            
            for b = 1:nIterations
                % Subsample trials to match size
                numPriorTrials   = numel(OnlyPrior_trials.trial);             % determine how many trials are avialable in each condition
                numPretoneTrials = numel(OnlyPretone_trials.trial);
                minTrials = min(numPriorTrials, numPretoneTrials);
            
                if numPriorTrials > numPretoneTrials
                    idx = randperm(numPriorTrials, minTrials);
                    OnlyPrior_trials_subset.label      = OnlyPrior_trials.label;
                    % OnlyPrior_trials_subset.time       = OnlyPrior_trials.time(1:minTrials);
                    OnlyPrior_trials_subset.time       = OnlyPrior_trials.time(idx);
                    OnlyPrior_trials_subset.trial      = OnlyPrior_trials.trial(idx);
                    OnlyPrior_trials_subset.fsample    = OnlyPrior_trials.fsample;
                    % OnlyPrior_trials_subset.sampleinfo = OnlyPrior_trials.sampleinfo(1:minTrials);
                    OnlyPrior_trials_subset.sampleinfo = OnlyPrior_trials.sampleinfo(idx, :);
                    OnlyPretone_trials_subset          = OnlyPretone_trials;
                else
                    idx = randperm(numPretoneTrials, minTrials);
                    OnlyPretone_trials_subset.label      = OnlyPretone_trials.label;
                    % OnlyPretone_trials_subset.time       = OnlyPretone_trials.time(1:minTrials);
                    OnlyPretone_trials_subset.time       = OnlyPretone_trials.time(idx);
                    OnlyPretone_trials_subset.trial      = OnlyPretone_trials.trial(idx);
                    OnlyPretone_trials_subset.fsample    = OnlyPretone_trials.fsample;
                    % OnlyPretone_trials_subset.sampleinfo = OnlyPretone_trials.sampleinfo(1:minTrials);
                    OnlyPretone_trials_subset.sampleinfo = OnlyPretone_trials.sampleinfo(idx, :);
                    OnlyPrior_trials_subset              = OnlyPrior_trials;
                end
            
                % --- Frequency analysis ---
            
                cfg            = [];
                cfg.method     = 'mtmfft';
                cfg.taper      = 'dpss';
                cfg.output     = 'fourier';
                cfg.tapsmofrq  = 4;
                cfg.pad        = 1;
                cfg.foilim     = [0 100];

                freq_OnlyPrior   = ft_freqanalysis(cfg, OnlyPrior_trials_subset);
                freq_OnlyPretone = ft_freqanalysis(cfg, OnlyPretone_trials_subset);

                % --- PSI ---
                cfg = [];
                cfg.method     = 'psi';
                cfg.bandwidth  = 2;
                cfg.channelcmb = OnlyPrior_Balanced_ChannelPairs;
                PSI_OnlyPrior = ft_connectivityanalysis(cfg, freq_OnlyPrior);
                PSI_OnlyPrior.dof       = length(OnlyPrior_trials_subset.trial);
                Nfreq = length(PSI_OnlyPrior.freq);
                PSI_OnlyPrior_all(:,:,b) = PSI_OnlyPrior.psispctrm;
            


                cfg.channelcmb          = OnlyPretone_Balanced_ChannelPairs;
                PSI_OnlyPretone         = ft_connectivityanalysis(cfg, freq_OnlyPretone);
                PSI_OnlyPretone.dof     = length(OnlyPretone_trials_subset.trial);
                Nfreq = length(PSI_OnlyPretone.freq);
                PSI_OnlyPretone_all(:,:,b) = PSI_OnlyPretone.psispctrm;
               

                % --- Coherence ---
                cfg = [];
                cfg.method     = 'coh';

                cfg.channelcmb = OnlyPrior_Balanced_ChannelPairs;
                Coh_OnlyPrior = ft_connectivityanalysis(cfg, freq_OnlyPrior);
                Coh_OnlyPrior_all(:,:,b) = Coh_OnlyPrior.cohspctrm;
             

                cfg.channelcmb = OnlyPretone_Balanced_ChannelPairs;
                Coh_OnlyPretone = ft_connectivityanalysis(cfg, freq_OnlyPretone);
                Coh_OnlyPretone_all(:,:,b) = Coh_OnlyPretone.cohspctrm;
         
            % --- Average across iterations ---
            PSI_OnlyPrior.psispctrm    = mean(PSI_OnlyPrior_all, 3);
            PSI_OnlyPrior.labelcmb     = OnlyPrior_Balanced_ChannelPairs;
            
            PSI_OnlyPretone.psispctrm  = mean(PSI_OnlyPretone_all, 3);
            PSI_OnlyPretone.labelcmb   = OnlyPretone_Balanced_ChannelPairs;
            
            Coh_OnlyPrior.cohspctrm    = mean(Coh_OnlyPrior_all, 3);
            Coh_OnlyPrior.labelcmb     = OnlyPrior_Balanced_ChannelPairs;
            
            Coh_OnlyPretone.cohspctrm  = mean(Coh_OnlyPretone_all, 3);
            Coh_OnlyPretone.labelcmb   = OnlyPretone_Balanced_ChannelPairs;
            
            % --- Save results ---
            save(fullfile(savedir, sprintf('%s_%s_OnlyPrior_correct_congruent_%s_PSI.mat', Animal, RecDate, Frequency_Band)), 'PSI_OnlyPrior');
            save(fullfile(savedir, sprintf('%s_%s_OnlyPretone_correct_congruent_%s_PSI.mat', Animal, RecDate, Frequency_Band)), 'PSI_OnlyPretone');
            
            save(fullfile(savedir, sprintf('%s_%s_OnlyPrior_correct_congruent_%s_Coh.mat', Animal, RecDate, Frequency_Band)), 'Coh_OnlyPrior');
            save(fullfile(savedir, sprintf('%s_%s_OnlyPretone_correct_congruent_%s_Coh.mat', Animal, RecDate, Frequency_Band)), 'Coh_OnlyPretone');
            
            % --- Clear between sessions ---
            clear_list = {
                'freq_OnlyPrior', 'freq_OnlyPretone', ...
                'OnlyPrior_trials_subset', 'OnlyPretone_trials_subset', ...
                'PSI_OnlyPrior', 'PSI_OnlyPretone', ...
                'Coh_OnlyPrior', 'Coh_OnlyPretone', ...
                'PSI_OnlyPrior_all', 'PSI_OnlyPretone_all', ...
                'Coh_OnlyPrior_all', 'Coh_OnlyPretone_all'};
            clear(clear_list{:});
            clear clear_list;


        
            end
     end 


%% 

