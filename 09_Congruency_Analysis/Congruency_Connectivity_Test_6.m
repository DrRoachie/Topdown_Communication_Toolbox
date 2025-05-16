
%% Congruency_Connectivity_Test

% This code tests how connectivity between PFC and AC changes as a function
% of congrency (whether or not the target matches the expectation).  
% When an informative LED light (ones that indicate a 75% chance of high or low), 
% the target can either be congruent (match) or incongruent (mismatch) 
% the tone predicted by the LED. When an noninformative LED is illuminated, 
% the target is compared to the subsequency pretones. 

% This analysis will be carried out session-wise on Correct Trials and on
% all SNRs.

% Both coherence and PSI will be calculated. 

% In order to circumvent monte carlo estimation, we decided to subsample
% the data, and just deal with the 40-60 percent loss in trials in the
% condition with more trials. 
 
% Statistics

% Because all heuristics are trials balanced, the magnitude of the two
% conditins are directly comparable via nonparametric testing. This also
% allows us to subtract the values between conditions which is useful for
% assessing how the two conditions change across sessions. 

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

        session_info = {
                        'MrM', '190719';
                        'MrCassius', '190720';
                        'MrM', '190722';
                        'MrCassius', '190723';
                        'MrCassius', '190725';};
             
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
        datadir           = fullfile('D:\04_Epoc_Cut', Animal, 'testTone', RecDate);                      % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('D:\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);    % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\Congruency_Analysis';                 % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

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

        % now we build a matrix of all of the shared channel pairs between
        % correct and wrong trials

        % Find common elements
        combined_PFC_chans =  unique([OnlyPrior_correct_PFC_chans; OnlyPretone_correct_PFC_chans], 'stable');
        combined_AC_chans  =  unique([OnlyPrior_correct_AC_chans;  OnlyPretone_correct_AC_chans],  'stable');

        % Get all pairwise combinations
        [PFC_comb, AC_comb] = ndgrid(combined_PFC_chans, combined_AC_chans);

         % Reshape into Nx2 cell array
         Shared_ChannelPairs = [PFC_comb(:), AC_comb(:)]; % these are the channel pairs that we are going to the analysis on, these will change as we loop through sessions.

         % clear variables used to generate shared channel pairs
         clear_list = {'AC_comb', 'combined_AC_chans', 'combined_PFC_chans', 'OnlyPrior_correct_files', 'OnlyPretone_correct_files','PFC_comb', 'OnlyPrior_correct_PFC_chans', 'OnlyPretone_correct_PFC_chans', ...
                       'OnlyPrior_correct_AC_chans', 'OnlyPretone_correct_AC_chans', 'OnlyPrior_correct_modifiedCellArray', 'OnlyPretone_correct_modifiedCellArray', ...
                       'OnlyPrior_correct_SigChans_fn', 'OnlyPretone_correct_SigChans_fn', 'OnlyPretone_correct_SigChans_fn', 'OnlyPrior_correct_SigChans', 'OnlyPretone_correct_SigChans'};

         clear(clear_list{:});
         clear clear_list;

         end        
     
         % STEP 2: Load in the data, and record numerical values to 'H' and
         % 'L' designation, and code each trials as congruent and
         % incongruent. 

         % extract session information
         baseFileName = sessions.name;
         fullFileName = fullfile(sessions.folder, baseFileName);
         fprintf(1, 'Now reading %s/n', fullFileName);
         V   = load(fullfile(fullFileName));
        

         % set parameters required for fieldtrip functions and data selection
         V_params.choice        = V.choice;
         V_params.err           = V.err;
         V_params.pretone       = V.pretone;
         V_params.pretoneLength = V.pretoneLength;
         V_params.prior         = cell2char(V.prior); 
         V_params.SNR           = V.SNR;

         % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

         V   = ConvertStim(V, RecDate);
         
         % Adds new field in data structure for OnlyPrior Trials that
         % indicates congruence between target and LED

             % Initialize OnlyPrior_data.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.prior));  % Preallocate as a cell array
             V.pretone    = cellstr(V.pretone);

             % Loop through each element in V.prior and V.target

                 for i = 1:length(V.prior)
                    if strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'H')
                        V.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'L')
                        V.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(V.prior{i}, 'N')
                       V.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'H'))
                       V.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

             V_params.congruency = V.congruency;

            % Select the OnlyPrior congruent trials 
            iSelect              = setStimulusCondition('Both');
            iSelect.err          = 'c';                                             
            iSelect.SNR          = SNR;                                                             % Select all trials regardless of SNR value
            iSelect.congruency   = 'congruent';                                 
            Congruent_trials     = selectData(V.data, V_params, iSelect);     % 


            %Select the OnlyPrior incongruent trials 
            iSelect              = setStimulusCondition('Both');
            iSelect.err          = 'c';                                   % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = 'incongruent';                                   % Select all trials regarless of congruency
            Incongruent_trials     = selectData(V.data, V_params, iSelect);     % 

            % Number of bootstrap iterations
            nIterations = 50;
            
            % Preallocate arrays to store coherence spectra
            % Assume size is [nPairs x nFreqs], e.g., 190 channel pairs by 101 frequencies
            nPairs = size(Shared_ChannelPairs, 1);
            nFreqs = 101;
            
            Coh_congruent_all   = zeros(nPairs, nFreqs, nIterations);
            Coh_incongruent_all = zeros(nPairs, nFreqs, nIterations);
            
            PSI_congruent_all   = zeros(nPairs, nFreqs, nIterations);
            PSI_incongruent_all = zeros(nPairs, nFreqs, nIterations);
            
            for b = 1:nIterations

                % Determine the number of trials in each condition
                numCongruentTrials   = numel(Congruent_trials.trial);
                numIncongruentTrials = numel(Incongruent_trials.trial);
                
                % Find the smaller trial count and which dataset is larger
                minTrials = min(numCongruentTrials, numIncongruentTrials);
                
                % Randomly subsample the larger dataset to match the smaller one,
                % and reconstruct a viable data array for field trip functions. 
                if numCongruentTrials > numIncongruentTrials
                    idx = randperm(numCongruentTrials, minTrials);
                    Congruent_trials_subset.label      = Congruent_trials.label;
                    Congruent_trials_subset.time       = Congruent_trials.time(1:minTrials);
                    Congruent_trials_subset.trial      = Congruent_trials.trial(idx);
                    Congruent_trials_subset.fsample    = Congruent_trials.fsample;
                    Congruent_trials_subset.sampleinfo = Congruent_trials.sampleinfo(1:minTrials);
                    
                    Incongruent_trials_subset          = Incongruent_trials;
        
                else
        
                    idx = randperm(numIncongruentTrials, minTrials);
                    Incongruent_trials_subset.label      = Incongruent_trials.label;
                    Incongruent_trials_subset.time       = Incongruent_trials.time(1:minTrials);
                    Incongruent_trials_subset.trial      = Incongruent_trials.trial(idx);
                    Incongruent_trials_subset.fsample    = Incongruent_trials.fsample;
                    Incongruent_trials_subset.sampleinfo = Incongruent_trials.sampleinfo(1:minTrials); 
        
                    Congruent_trials_subset               = Congruent_trials;
        
                end
        
                % Calculate PSI for both Conditions 
        
                    cfg            = [];
                    cfg.method     = 'mtmfft';
                    cfg.taper      = 'dpss';
                    cfg.output     = 'fourier';
                    cfg.tapsmofrq  = 4; 
                    cfg.pad        = 1;
                    cfg.foilim     = [0 100];
        
                    freq_congruent      = ft_freqanalysis(cfg, Congruent_trials_subset);
                    freq_incongruent    = ft_freqanalysis(cfg, Incongruent_trials_subset);
        
                    cfg             = [];
                    cfg.method      = 'psi';
                    cfg.bandwidth   = 2;
                    cfg.channelcmb  = Shared_ChannelPairs; % representing the all shared channel pairs for the 
        
                    PSI_congruent         = ft_connectivityanalysis(cfg, freq_congruent);
                    PSI_congruent.dof     = length(Congruent_trials_subset.trial);
                    PSI_incongruent       = ft_connectivityanalysis(cfg, freq_incongruent);
                    PSI_incongruent.dof   = length(Incongruent_trials_subset.trial);
                    
                    % Store the cohspectrm values
                    PSI_congruent_all(:, :, b)   = PSI_congruent.psispctrm;
                    PSI_incongruent_all(:, :, b) = PSI_incongruent.psispctrm;

                    % Calculate Coherence for both Conditions 
                    % non-parametric computation of the cross-spectral density matrix 
        
                    cfg            = [];
                    cfg.method     = 'coh';
                    cfg.channelcmb = Shared_ChannelPairs;   
        
                    Coh_congruent       = ft_connectivityanalysis(cfg, freq_congruent);
                    Coh_incongruent     = ft_connectivityanalysis(cfg, freq_incongruent);

                    % Store the cohspectrm values
                    Coh_congruent_all(:, :, b)   = Coh_congruent.cohspctrm;
                    Coh_incongruent_all(:, :, b) = Coh_incongruent.cohspctrm;

            end

                    % average the boostrapped subsample 
                    PSI_congruent.psispctrm    = mean(PSI_congruent_all, 3);
                    PSI_incongruent.psispctrm  = mean(PSI_incongruent_all, 3);
                    Coh_congruent.cohspctrm    = mean(Coh_congruent_all, 3);
                    Coh_incongruent.cohspctrm  = mean(Coh_incongruent_all, 3);

                    % save session value 
                    save_file_name = sprintf('%s_%s_%s_%s_%s_PSI.mat', Animal, RecDate,  'correct', 'congruent', Frequency_Band);
                    save(fullfile(savedir, save_file_name), 'PSI_congruent');
                    save_file_name = sprintf('%s_%s_%s_%s_%s_PSI.mat', Animal, RecDate, 'correct', 'incongruent', Frequency_Band);
                    save(fullfile(savedir, save_file_name), 'PSI_incongruent')
             
                    save_file_name = sprintf('%s_%s_%s_%s_%s_Coh.mat', Animal, RecDate, 'correct', 'congruent', Frequency_Band);
                    save(fullfile(savedir, save_file_name), 'Coh_congruent');
                    save_file_name = sprintf('%s_%s_%s_%s_%s_Coh.mat', Animal, RecDate, 'correct', 'incongruent', Frequency_Band);
                    save(fullfile(savedir, save_file_name), 'Coh_incongruent');

                    clear_list = {'freq_congruent', 'freq_incongruent', 'cfg','PSI_congruent', 'PSI_incongruent', 'Coh_congruent', 'Coh_incongruent'};
                    clear(clear_list{:});
                    clear clear_list;


     end

