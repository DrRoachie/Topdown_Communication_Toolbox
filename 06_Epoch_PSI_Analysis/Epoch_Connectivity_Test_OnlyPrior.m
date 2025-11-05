
%% Epoch_Connectivity_Test

% This code tests how connectivity between PFC and AC changes as a function
% of epoch (LED illumination versus testTone presentation). For each
% session, PSI during preCueOnset and testTone onset is calculated. Then
% the session-wise PSI values are visualized and statistically compared.
% Though code for OnlyPretone Trials (neural LED followed by pretones) is
% present, we chose to do this analysis on OnlyPrior trials (informatative
% LED followed by silence). 
 
% Step 1: Get all the channels that are significantly modulated by the
% presentation of the informative LED ('OnlyPrior Trials, preCueOnset')
% regardless of Animal ('MrM' and 'MrCassius', SNR ('AllSNR') and Behavior
% ('Correct' and 'Wrong' Trials). These Shared_Channels will be used in the
% calculation of both 

% Step 2: Calculate the PSI for each significant channel pair for each Epoch. 

% Note: Because we are testing the effect of epoch, the number of trials
% should be the same between the two conditions ('preCueOnset' and
% 'TestToneOnset'), and thus the PSI values are directly comparable.
% However, Lalitta's preprocessed data has different trials as a function
% of epoch. I am not sure why this is...here, we just truncate the larger
% data frame to match the smaller one. 

% Step 3: Organize PSI values by channel pair. 

% Step 4: Average the PSI values 

% Step 5: Plot the mean PSI values per channel pair (heat maps)

% Step 6a: Generate histograms of all PSI values (regardless of spatial location) reflecting directionality 

% Step 6b: Between conditions, subtract the PSI values of a given
% direction. Make a histogram of the residuals. 

% Step 7: Test if either direction deviates from zero, and also if the two
% directions deviate from each other. 

%% Set parameters for the preCue Condition

Condition         = 'OnlyPrior';                                                  % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];   % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 

%% Calculate PSI Spectra during the preCueOnset epoch; PriorOnly Trials, Collapsed Across Behavior (Correct and Wrong) and SNR (High and Low)
% ahead of the development of this code, we found that there was not a
% strong influence of Behavior/SNR on PSI, so we collapse across them here.
% 

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
    end


     
    if strcmp(Frequency_Band, 'beta')

        session_info = {'MrCassius', '190404';
                        'MrCassius', '190413';
                        'MrCassius', '190416';
                        'MrCassius', '190418';
                        'MrCassius', '190419';
                        'MrM', '190422';
                        'MrCassius', '190517';
                        'MrM', '190525';
                        'MrM', '190527';
                        'MrM', '190601';
                        'MrCassius', '190703';
                        'MrCassius', '190713';
                        'MrCassius', '190718';
                        'MrM', '190719';
                        'MrCassius', '190723';};

     end

     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 

        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir = fullfile('D:\04_Epoc_Cut_bipolar', Animal, 'preCue', RecDate);                                 % epoch data, with the time chunk that you want isolated
        chandir = fullfile('D:\05_Significant_Channels_bipolar', Animal, 'preCue', RecDate);     % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')  % make folders if they don't exist already
           mkdir(savedir);
        end   

        preCue_correct_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'preCueOnset', 'OnlyPrior', 'Correct', Frequency_Band);
        preCue_correct_files = dir(fullfile(chandir, preCue_correct_SigChans_fn));

        if  isempty(preCue_correct_files)

            % generate null output files if there are no sig. channels for Condition_1
            % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, 'OnlyPrior', RecDate, Frequency_Band, 'Correct'));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');

        else % get the labels for the PFC channels and AC channels for correct trials 

        preCue_correct_SigChans = load(fullfile(chandir, preCue_correct_SigChans_fn), 'significant_channels');
        preCue_correct_modifiedCellArray   = regexprep(preCue_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        preCue_correct_PFC_chans  =  preCue_correct_modifiedCellArray(startsWith(preCue_correct_modifiedCellArray, '*PFC'));
        preCue_correct_AC_chans   =  preCue_correct_modifiedCellArray(startsWith(preCue_correct_modifiedCellArray, '*AC'));

        end

        % get the labels for the PFC channels and AC channels for wrong trials 

        preCue_wrong_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'preCueOnset', 'OnlyPrior', 'Wrong', Frequency_Band);
        preCue_wrong_files = dir(fullfile(chandir, preCue_wrong_SigChans_fn));


        if  isempty(preCue_wrong_files)

          % generate null output files if there are no sig. channels for Condition_2
          % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

          save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, 'preCueOnset', Frequency_Band,  'Wrong'));
          fid4 = fopen(save_file_name, 'w');
          fclose('all');

        else

        preCue_wrong_SigChans = load(fullfile(chandir, preCue_wrong_SigChans_fn), 'significant_channels');
        preCue_wrong_modifiedCellArray   = regexprep(preCue_wrong_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        preCue_wrong_PFC_chans  =  preCue_wrong_modifiedCellArray(startsWith(preCue_wrong_modifiedCellArray, '*PFC'));
        preCue_wrong_AC_chans   =  preCue_wrong_modifiedCellArray(startsWith(preCue_wrong_modifiedCellArray, '*AC'));

        % now we build a matrix of all of the shared channel pairs between
        % correct and wrong trials

        % Find common elements
        combined_PFC_chans =  unique([preCue_correct_PFC_chans; preCue_wrong_PFC_chans], 'stable');
        combined_AC_chans  =  unique([preCue_correct_AC_chans;  preCue_wrong_AC_chans],  'stable');

        % Get all pairwise combinations
        [PFC_comb, AC_comb] = ndgrid(combined_PFC_chans, combined_AC_chans);

         % Reshape into Nx2 cell array
         Shared_ChannelPairs = [PFC_comb(:), AC_comb(:)]; % these are the channel pairs that we are going to the analysis on, these will change as we loop through sessions.

         % clear variables used to generate shared channel pairs
         clear_list = {'AC_comb', 'combined_AC_chans', 'combined_PFC_chans', 'preCue_correct_files', 'preCue_wrong_files','PFC_comb', 'preCue_correct_PFC_chans', 'preCue_wrong_PFC_chans', ...
                       'preCue_correct_AC_chans', 'preCue_wrong_AC_chans', 'preCue_correct_modifiedCellArray', 'preCue_wrong_modifiedCellArray', ...
                       'preCue_correct_SigChans_fn', 'preCue_wrong_SigChans_fn', 'preCue_correct_SigChans', 'preCue_wrong_SigChans'};

         clear(clear_list{:});
         clear clear_list;

         end        

         % STEP 2: Load preCue data and calculate PSI

         % extract session information
         baseFileName = sessions.name;
         fullFileName = fullfile(sessions.folder, baseFileName);
         fprintf(1, 'Now reading %s/n', fullFileName);
         V = load(fullfile(fullFileName));

         % set parameters required for fieldtrip functions and data selection
         params.choice = V.choice;
         params.err = V.err;
         params.pretone = V.pretone;
         params.pretoneLength = V.pretoneLength;
         params.prior = cell2char(V.prior); 
         params.SNR = V.SNR;

         % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

         V = ConvertStim(V, RecDate);

         % Adds new field in data structure for OnlyPrior Trials that
         % indicates congruence between target and LED

            if strcmp(Condition,'OnlyPrior') == 1

             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.prior));  % Preallocate as a cell array

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

             params.congruency = V.congruency;

            end

            % Adds new field in data structure for OnlyPretone Trials that
            % indicates congruence between target and pretone

            if strcmp(Condition,'OnlyPretone') == 1

             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
             V.pretone = cellstr(V.pretone);

             % Loop through each element in V.pretone and V.target

                 for i = 1:length(V.pretone)
                    if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                        V.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                        V.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(V.pretone{i}, 'N')
                        V.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                        V.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

             params.congruency = V.congruency;

            end

            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = [];                                    % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = []';                                   % Select all trials regarless of congruency
            preCueOnset_data     = selectData(V.data,params,iSelect);     % 

             % Step 2: Load testToneOnset data just for the purpose of the making
             % sure that PSI is calculated using the same number of trials.

            % Set directory and get a list of all files in the folder with the desired file name pattern.
            datadir = fullfile('D:\04_Epoc_Cut_bipolar', Animal, 'testTone', RecDate);                             % epoch data, with the time chunk that you want isolated
            chandir = fullfile('D:\05_Significant_Channels_bipolar', Animal, 'testTone', RecDate);           % the sig channels that are output LFP_Spectral_Analysis
            savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

            sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
            addpath(genpath(datadir));
            addpath(genpath(chandir));

            % make save folders and directory
            savedir = fullfile(savedir, RecDate, Animal);

                if ~exist(savedir, 'dir')  % make folders if they don't exist already
                   mkdir(savedir);
                end   

             % extract session information
             baseFileName = sessions.name;
             fullFileName = fullfile(sessions.folder, baseFileName);
             fprintf(1, 'Now reading %s/n', fullFileName);
             V = load(fullfile(fullFileName));

             % set parameters required for fieldtrip functions and data selection
             params.choice = V.choice;
             params.err = V.err;
             params.pretone = V.pretone;
             params.pretoneLength = V.pretoneLength;
             params.prior = cell2char(V.prior); 
             params.SNR = V.SNR;

             % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

             V = ConvertStim(V, RecDate);

             % Adds new field in data structure for OnlyPrior Trials that
             % indicates congruence between target and LED

                if strcmp(Condition,'OnlyPrior') == 1

                     % Initialize V.congruency as a cell array as a new selection factor for select data function  

                     V.congruency = cell(size(V.prior));  % Preallocate as a cell array

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

                     params.congruency = V.congruency;

                end

            % Adds new field in data structure for OnlyPretone Trials that
            % indicates congruence between target and pretone

                if strcmp(Condition,'OnlyPretone') == 1

                 % Initialize V.congruency as a cell array as a new selection factor for select data function  

                 V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
                 V.pretone = cellstr(V.pretone);

                 % Loop through each element in V.pretone and V.target

                     for i = 1:length(V.pretone)
                        if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                            V.congruency{i} = 'congruent';  % If both are 'H'
                        elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                            V.congruency{i} = 'congruent';  % If both are 'L'
                        elseif strcmp(V.pretone{i}, 'N')
                            V.congruency{i} = 'neutral';  % If prior is 'N'
                        elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                            V.congruency{i} = 'incongruent';  % If prior and target are different
                        end
                     end

                 params.congruency = V.congruency;

                 end

            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = [];                                    % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = []';                                   % Select all trials regarless of congruency
            testToneOnset_data   = selectData(V.data,params,iSelect);   

            % Truncate the larger data frame to match smaller values. This
            %  ensures that the testTone epoch and the preCue
            % epoch have the same number of trials.

            numTrials_preCue     = size(preCueOnset_data.trial, 2);
            numTrials_testTone   = size(testToneOnset_data.trial, 2);
            minTrials            = min(numTrials_preCue, numTrials_testTone);

            preCueOnset_data.trial      = preCueOnset_data.trial(:, 1:minTrials);
            preCueOnset_data.time       = preCueOnset_data.time(:, 1:minTrials);
            preCueOnset_data.sampleinfo = preCueOnset_data.sampleinfo(1:minTrials, :);

            testToneOnset_data.trial      = testToneOnset_data.trial(:, 1:minTrials);
            testToneOnset_data.time       = testToneOnset_data.time(:, 1:minTrials);
            testToneOnset_data.sampleinfo = testToneOnset_data.sampleinfo(1:minTrials, :);

            % calculate PSI for the preCue epoch

            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.taper      = 'dpss';
            cfg.output     = 'fourier';
            cfg.tapsmofrq  = 4; 
            cfg.pad        = 1;
            cfg.foilim     = [0 100];

            freq_preCueOnset       = ft_freqanalysis(cfg, preCueOnset_data);

            cfg             = [];
            cfg.method      = 'psi';
            cfg.bandwidth   = 2;
            cfg.channelcmb  = Shared_ChannelPairs; % representing the all shared channel pairs for the 

            PSI_preCueOnset         = ft_connectivityanalysis(cfg, freq_preCueOnset);
            PSI_preCueOnset.dof     = length(preCueOnset_data.trial);

            save_file_name = sprintf('%s_%s_%s_%s_PSI_data.mat', Animal, RecDate, 'preCueOnset', Frequency_Band);
            save(fullfile(savedir, save_file_name), 'PSI_preCueOnset');

            clear_list = {' freq_preCueOnset',  'cfg','PSI_preCueOnset'};
            clear(clear_list{:});
            clear clear_list;

    end

%% clear everything 

clearvars -except Condition Frequency_Band SNR 


%% Calculate PSI Spectra during the testToneOnset epoch; PriorOnly Trials, Collapsed Across Behavior (Correct and Wrong) and SNR (High and Low)

    % Step 1: Get list of modulated channels 

    if strcmp(Frequency_Band, 'theta') || strcmp(Frequency_Band, 'alpha')

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

    end

    

    if strcmp(Frequency_Band, 'beta')

    session_info = {'MrCassius', '190404';
                    'MrCassius', '190413';
                    'MrCassius', '190416';
                    'MrCassius', '190418';
                    'MrCassius', '190419';
                    'MrM', '190422';
                    'MrCassius', '190517';
                    'MrM', '190525';
                    'MrM', '190527';
                    'MrM', '190601';
                    %'MrCassius', '190703';
                    'MrCassius', '190713';
                    'MrCassius', '190718';
                    'MrM', '190719';
                    'MrCassius', '190723';};

     end

    for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 

        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir = fullfile('D:\04_Epoc_Cut_bipolar', Animal, 'testTone', RecDate);                                 % epoch data, with the time chunk that you want isolated
        chandir = fullfile('D:\05_Significant_Channels_bipolar', Animal, 'testTone', RecDate);     % the sig channels that are output LFP_Spectral_Analysis
        savedir = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);
        if ~exist(savedir, 'dir')  % make folders if they don't exist already
           mkdir(savedir);
        end   

        testTone_correct_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'testToneOnset', 'OnlyPrior', 'Correct', Frequency_Band);
        testTone_correct_files = dir(fullfile(chandir, testTone_correct_SigChans_fn));

        if  isempty(testTone_correct_files)

            % generate null output files if there are no sig. channels for Condition_1
            % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, 'OnlyPrior', RecDate, Frequency_Band, 'Correct'));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');

        else % get the labels for the PFC channels and AC channels 


        testTone_correct_SigChans = load(fullfile(chandir, testTone_correct_SigChans_fn), 'significant_channels');
        testTone_correct_modifiedCellArray   = regexprep(testTone_correct_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        testTone_correct_PFC_chans  =  testTone_correct_modifiedCellArray(startsWith(testTone_correct_modifiedCellArray, '*PFC'));
        testTone_correct_AC_chans   =  testTone_correct_modifiedCellArray(startsWith(testTone_correct_modifiedCellArray, '*AC'));

        end

        testTone_wrong_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, 'testToneOnset', 'OnlyPrior', 'Wrong', Frequency_Band);
        testTone_wrong_files = dir(fullfile(chandir, testTone_wrong_SigChans_fn));


        if  isempty(testTone_wrong_files)

          % generate null output files if there are no sig. channels for Condition_2
          % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

          save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, 'testToneOnset', Frequency_Band,  'Wrong'));
          fid4 = fopen(save_file_name, 'w');
          fclose('all');

        else

        testTone_wrong_SigChans = load(fullfile(chandir, testTone_wrong_SigChans_fn), 'significant_channels');
        testTone_wrong_modifiedCellArray   = regexprep(testTone_wrong_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

        % Separate elements that start with *PFC and *AC
        testTone_wrong_PFC_chans  =  testTone_wrong_modifiedCellArray(startsWith(testTone_wrong_modifiedCellArray, '*PFC'));
        testTone_wrong_AC_chans   =  testTone_wrong_modifiedCellArray(startsWith(testTone_wrong_modifiedCellArray, '*AC'));

        % now we build a matrix of all of the shared channel pairs 

        % Find common elements
        combined_PFC_chans =  unique([testTone_correct_PFC_chans; testTone_wrong_PFC_chans], 'stable');
        combined_AC_chans  =  unique([testTone_correct_AC_chans;  testTone_wrong_AC_chans],  'stable');

        % Get all pairwise combinations
        [PFC_comb, AC_comb] = ndgrid(combined_PFC_chans, combined_AC_chans);

         % Reshape into Nx2 cell array
         Shared_ChannelPairs = [PFC_comb(:), AC_comb(:)];

         % clear variables used to generate shared channel pairs
         clear_list = {'AC_comb', 'combined_AC_chans', 'combined_PFC_chans', 'testTone_correct_files', 'testTone_wrong_files','PFC_comb', 'testTone_correct_PFC_chans', 'testTone_wrong_PFC_chans', ...
                       'testTone_correct_AC_chans', 'testTone_wrong_AC_chans', 'testTone_correct_modifiedCellArray', 'testTone_wrong_modifiedCellArray', ...
                       'testTone_correct_SigChans_fn', 'testTone_wrong_SigChans_fn', 'testTone_correct_SigChans', 'testTone_wrong_SigChans'};

         clear(clear_list{:});
         clear clear_list;

         end        

         % Load preCue data 

        % Need the preCue data to ensure that the Epochs are the same
        % number of trials 
   
        datadir = fullfile('D:\04_Epoc_Cut_bipolar', Animal, 'preCue', RecDate);                                 % epoch data, with the time chunk that you want isolated
        chandir = fullfile('D:\05_Significant_Channels_bipolar', Animal, 'preCue', RecDate);     % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

        sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
        addpath(genpath(datadir));
        addpath(genpath(chandir));

        % make save folders and directory
        savedir = fullfile(savedir, RecDate, Animal);

                if ~exist(savedir, 'dir')  % make folders if they don't exist already
                   mkdir(savedir);
                end   


         % extract session information
         baseFileName = sessions.name;
         fullFileName = fullfile(sessions.folder, baseFileName);
         fprintf(1, 'Now reading %s/n', fullFileName);
         V = load(fullfile(fullFileName));

         % set parameters required for fieldtrip functions and data selection
         params.choice = V.choice;
         params.err = V.err;
         params.pretone = V.pretone;
         params.pretoneLength = V.pretoneLength;
         params.prior = cell2char(V.prior); 
         params.SNR = V.SNR;

         % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

         V = ConvertStim(V, RecDate);

         % Adds new field in data structure for OnlyPrior Trials that
         % indicates congruence between target and LED

            if strcmp(Condition,'OnlyPrior') == 1

             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.prior));  % Preallocate as a cell array

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

             params.congruency = V.congruency;

            end

            % Adds new field in data structure for OnlyPretone Trials that
            % indicates congruence between target and pretone

            if strcmp(Condition,'OnlyPretone') == 1

             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
             V.pretone = cellstr(V.pretone);

             % Loop through each element in V.pretone and V.target

                 for i = 1:length(V.pretone)
                    if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                        V.congruency{i} = 'congruent';  % If both are 'H'
                    elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                        V.congruency{i} = 'congruent';  % If both are 'L'
                    elseif strcmp(V.pretone{i}, 'N')
                        V.congruency{i} = 'neutral';  % If prior is 'N'
                    elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                        V.congruency{i} = 'incongruent';  % If prior and target are different
                    end
                 end

             params.congruency = V.congruency;

            end

            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = [];                                    % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = []';                                   % Select all trials regarless of congruency
            preCueOnset_data     = selectData(V.data,params,iSelect);     % 

             % Load testToneOnset data 

            % Set directory and get a list of all files in the folder with the desired file name pattern.
            datadir = fullfile('D:\04_Epoc_Cut_bipolar', Animal, 'testTone', RecDate);                             % epoch data, with the time chunk that you want isolated
            chandir = fullfile('D:\05_Significant_Channels_bipolar', Animal, 'testTone', RecDate);           % the sig channels that are output LFP_Spectral_Analysis
            savedir = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

            sessions = dir(fullfile(datadir,'*.mat')); % bad naming convention but only ever one session at a time
            addpath(genpath(datadir));
            addpath(genpath(chandir));

            % make save folders and directory
            savedir = fullfile(savedir, RecDate, Animal);

                if ~exist(savedir, 'dir')  % make folders if they don't exist already
                   mkdir(savedir);
                end   

             % extract session information
             baseFileName = sessions.name;
             fullFileName = fullfile(sessions.folder, baseFileName);
             fprintf(1, 'Now reading %s/n', fullFileName);
             V = load(fullfile(fullFileName));

             % set parameters required for fieldtrip functions and data selection
             params.choice = V.choice;
             params.err = V.err;
             params.pretone = V.pretone;
             params.pretoneLength = V.pretoneLength;
             params.prior = cell2char(V.prior); 
             params.SNR = V.SNR;

             % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 

             V = ConvertStim(V, RecDate);

             % Adds new field in data structure for OnlyPrior Trials that
             % indicates congruence between target and LED

                if strcmp(Condition,'OnlyPrior') == 1

                     % Initialize V.congruency as a cell array as a new selection factor for select data function  

                     V.congruency = cell(size(V.prior));  % Preallocate as a cell array

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

                     params.congruency = V.congruency;

                end

            % Adds new field in data structure for OnlyPretone Trials that
            % indicates congruence between target and pretone

                if strcmp(Condition,'OnlyPretone') == 1

                 % Initialize V.congruency as a cell array as a new selection factor for select data function  

                 V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
                 V.pretone = cellstr(V.pretone);

                 % Loop through each element in V.pretone and V.target

                     for i = 1:length(V.pretone)
                        if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                            V.congruency{i} = 'congruent';  % If both are 'H'
                        elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                            V.congruency{i} = 'congruent';  % If both are 'L'
                        elseif strcmp(V.pretone{i}, 'N')
                            V.congruency{i} = 'neutral';  % If prior is 'N'
                        elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                            V.congruency{i} = 'incongruent';  % If prior and target are different
                        end
                     end

                 params.congruency = V.congruency;

                 end

            iSelect              = setStimulusCondition('OnlyPrior');
            iSelect.err          = [];                                    % Select both correct and wrong trial
            iSelect.SNR          = SNR;                                   % Select all trials regardless of SNR value
            iSelect.congruency   = []';                                   % Select all trials regarless of congruency
            testToneOnset_data   = selectData(V.data,params,iSelect);   

            % Truncate the larger data frame to match smaller values 

            numTrials_preCue     = size(preCueOnset_data.trial, 2);
            numTrials_testTone   = size(testToneOnset_data.trial, 2);
            minTrials            = min(numTrials_preCue, numTrials_testTone);

            preCueOnset_data.trial      = preCueOnset_data.trial(:, 1:minTrials);
            preCueOnset_data.time       = preCueOnset_data.time(:, 1:minTrials);
            preCueOnset_data.sampleinfo = preCueOnset_data.sampleinfo(1:minTrials, :);

            testToneOnset_data.trial      = testToneOnset_data.trial(:, 1:minTrials);
            testToneOnset_data.time       = testToneOnset_data.time(:, 1:minTrials);
            testToneOnset_data.sampleinfo = testToneOnset_data.sampleinfo(1:minTrials, :);

            % calculate PSI for both epochs 

            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.taper      = 'dpss';
            cfg.output     = 'fourier';
            cfg.tapsmofrq  = 4; 
            cfg.pad        = 1;
            cfg.foilim     = [0 100];

            freq_testToneOnset     = ft_freqanalysis(cfg, testToneOnset_data);

            cfg             = [];
            cfg.method      = 'psi';
            cfg.bandwidth   = 2;
            cfg.channelcmb  = Shared_ChannelPairs; 

            PSI_testToneOnset       = ft_connectivityanalysis(cfg, freq_testToneOnset);
            PSI_testToneOnset.dof   = length(testToneOnset_data.trial);

            save_file_name = sprintf('%s_%s_%s_%s_PSI_data.mat', Animal, RecDate, 'testToneOnset', Frequency_Band);
            save(fullfile(savedir, save_file_name), 'PSI_testToneOnset');

            clear_list = {'freq_testToneOnset',  'cfg','PSI_testToneOnset'};
            clear(clear_list{:});
            clear clear_list;       

    end

