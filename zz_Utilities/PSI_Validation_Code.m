%% Expectation_Connectivity_Test

% This code aims to validate and test how to inturpret PSI values. 
% In short, it is difficult to tell based off the field trip documentation
% how to know which direction corresponds to negative values and positive
% values. 

% Here we compare real data to simulated data to ascertain how field trip
% assigns directionality to our data structure. 

% Step 1 is to run our own data. We compare the field trip output structure
% between using cfg.channelcmb and running the code with all channels. With
% cfg.channelcmb it will populate a visualization with a PSI values, but
% due to the lay out it does not tell you which direction is which. When
% you run all channels you find that value in the upper right quadrant
% where negative values denote AC to PFC communication (because PFC was
% entered first) 

% Step 2 we validate Step 1 using simulated data, where two AC channels
% lead two PFC channels by 100 ms. In this case, we know the directionality
% of the interaction and we can correctly assign direction to the sign of
% the PSI values. 


%% Step 1

Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];   % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 

    % Step 1: Get list of modulated channels 

    % Because theta/alpha & beta have the chance to use different channels,
    % you have to make the code Frequency_Band dependent. 

    if strcmp(Frequency_Band, 'theta') || strcmp(Frequency_Band, 'alpha') % theta and alpha have the same shared pairs according to LFP_spec_eval, so they can use the same set/

        session_info = {'MrCassius', '190330';};
  
    end


     for rd = 1:length(session_info(:,1))

        % Establish the subset of data that you want to process. 
        Animal           = session_info{rd,1};                  % Options: 'MrCassius', 'MrM'; 
        RecDate          = session_info{rd,2};                  % Options: 'YYMMDD'; 

        % Set directory and get a list of all files in the folder with the desired file name pattern.
        datadir           = fullfile('D:\04_Epoc_Cut', Animal, 'testTone', RecDate);                                 % epoch data, with the time chunk that you want isolated
        chandir           = fullfile('D:\05_Significant_Channels_epoch', Animal, 'testTone', RecDate);     % the sig channels that are output LFP_Spectral_Analysis
        savedir           = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Connectivity_Whole';  % Path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)

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


               % Calculate PSI for both Conditions 
        
                    cfg            = [];
                    cfg.method     = 'mtmfft';
                    cfg.taper      = 'dpss';
                    cfg.output     = 'fourier';
                    cfg.tapsmofrq  = 4; 
                    cfg.pad        = 1;
                    cfg.foilim     = [0 100];
        
                    freq_OnlyPrior      = ft_freqanalysis(cfg, OnlyPrior_trials);
                    freq_OnlyPretone    = ft_freqanalysis(cfg, OnlyPretone_trials);
        
                    cfg             = [];
                    cfg.method      = 'psi';
                    cfg.bandwidth   = 2;
                    PSI_OnlyPrior         = ft_connectivityanalysis(cfg, freq_OnlyPrior);

                    cfg             = [];
                    cfg.method      = 'psi';
                    cfg.bandwidth   = 2;
                    %cfg.channelcmb  = Shared_ChannelPairs; % representing the all shared channel pairs for the 
                    PSI_OnlyPrior_channelcmb         = ft_connectivityanalysis(cfg, freq_OnlyPrior);


     end

%% visualization

cfg = [];
cfg.parameter = 'psispctrm';
ft_connectivityplot(cfg, PSI_OnlyPrior );

cfg = [];
cfg.parameter = 'psispctrm';
ft_connectivityplot(cfg, PSI_OnlyPrior_channelcmb  );



%%

% Visualization indicates that the correct channel pairs were processed,
% according to channelcmb parameter. And that actually the visualization
% strategy just duplicates the requested directionality. 

%%

pfc5_ac3 = PSI_OnlyPrior_channelcmb.psispctrm  (3, 21, 1);  % PFC-5 is at index 5, AC-3 is at index 21 ~ this one matches my data
ac3_pfc5 = PSI_OnlyPrior_channelcmb.psispctrm  (21, 3, 1);  % PFC-5 is at index 5, AC-3 is at index 21 


%% Step 2: simulate data

fsample     = 1000;          % sampling rate (Hz)
nsample     = fsample * 30;  % total samples for 30 s of data

% simulate some data
fsample = 1000;
nsample = fsample*30;
 
dat = randn(1,nsample+100);
 
% data.trial{1}(1,:) = dat(1,1:nsample) + 0.1.*randn(1,nsample);
% data.trial{1}(2,:) = dat(1,1:nsample) + 0.1.*randn(1,nsample);
% data.trial{1}(3,:) = dat(1,100+(1:nsample))+ 0.1.*randn(1,nsample);
% data.trial{1}(4,:) = dat(1,100+(1:nsample))+ 0.1.*randn(1,nsample);
data.trial{1}(1,:) = dat(1,1:nsample);
data.trial{1}(2,:) = dat(1,1:nsample);
data.trial{1}(3,:) = dat(1,100+(1:nsample));
data.trial{1}(4,:) = dat(1,100+(1:nsample));
data.time{1} = (1:nsample)./fsample;
data.label   = {'D1_PFC_ch05';'D1_PFC_ch06';'D2_AC_ch03';'D2_AC_ch04'}; % channel 2 is leading channel 1

data.fsample = fsample;

% cut into 2 second snippets
cfg = [];
cfg.length = 2;
data = ft_redefinetrial(cfg, data);

% spectral decomposition
  cfg            = [];
                    cfg.method     = 'mtmfft';
                    cfg.taper      = 'dpss';
                    cfg.output     = 'fourier';
                    cfg.tapsmofrq  = 4; 
                
                    cfg.foilim     = [0 100];

freq = ft_freqanalysis(cfg, data);
 
%% connectivity estimation
cfg = [];
cfg.method = 'psi';
cfg.bandwidth = 5;
cfg.channelcmb = {'D1_PFC_ch05','D2_AC_ch03'; 'D1_PFC_ch06','D2_AC_ch03';'D1_PFC_ch05', 'D2_AC_ch04';'D1_PFC_ch06','D2_AC_ch04'}; 
psi = ft_connectivityanalysis(cfg, freq);
 
%% visualization
cfg = [];
cfg.parameter = 'psispctrm';
ft_connectivityplot(cfg, psi);

%%


pfc1_ac1 = psi.psispctrm  (1, 3, :);  % PFC-5 is at index 5, AC-3 is at index 21 ~ this one matches my data

