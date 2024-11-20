%% Script Summary

% Takes the pre-processed Taku data, and for each experimental condition
% (PriorOnly and PretoneOnly) calculates the coherence and granger spectra
% and tests if they are significantly different from one another. 

% Version 2 is built for speed. Basically we use the wrapper to only get
% the Xcorr and Maris results, no plotting, no saving big chunks of data. 


% Establish what sessions and conditions you want to process, we piecewise
% run the 5-6 sessions at a time to prevent critical time loss during
% processing. 

session_info = {'MrCassius', '190605' 'preCueOnset', 'Correct', {'alpha'};
                'MrCassius', '190703' 'preCueOnset', 'Correct', {'alpha'}};

% loop 1: Run for all sessions
for rd = 1:length(session_info(:,1))

    % Establish the subset of data that you want to process. 
    Animal           = session_info{rd,1};                % Options: 'MrCassius', 'MrM'; 
    RecDate          = session_info{rd,2};                % Options: 'YYMMDD'; 
    Epoch            = session_info{rd,3};                % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
    Behavior         = session_info{rd,4};                % Options:  'Correct', Wrong, due to memory constraints                          
    RandomIteration  = 300;                               % Options: the number of times you want to allocate the random partition, and generate the test statistic
    freq_band_list   = session_info{rd,5};
    threshold_value  = 0.5;       % the threshold used for the monte carlo test between PriorOnly and PretoneOnly Conditions
    
    % Set directory and get a list of all files in the folder with the desired file name pattern.
    datadir = fullfile('D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2', Animal, extractBefore(Epoch, 'Onset'), RecDate); % epoch data, with the time chunk that you want isolated
    chandir = fullfile('D:\00_Significant_chans', Animal, extractBefore(Epoch, 'Onset'), RecDate);                                             % the sig channels that are output LFP_Spectral_Analysis
    specdir = fullfile('D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed', Animal, Epoch);                             % this is where the uncut preprocessed data is stored; you need the full time series for the spectogram
    savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\';                                        % path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
    
    sessions = dir(fullfile(datadir,'*.mat')); % bad naming conventionm but only ever one session at a time
    addpath(genpath(datadir));
    addpath(genpath(chandir));
    
    % make save folders and directory
    savedir = fullfile(savedir, RecDate, Animal, extractBefore(Epoch, 'Onset'), Behavior);
    if ~exist(savedir, 'dir')  % make folders if they don't exist already
       mkdir(savedir);
    end
      
    % loop 2: for every frequency band run the following 
    for fb = 1:length(freq_band_list)
        
        Frequency_Band = freq_band_list{fb};
        
        % Fetch the OnlyPrior Channels (output of ChanWise_SpecEval)
        
        PriorOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, 'OnlyPrior', Behavior, Frequency_Band);
        PriorOnly_files = dir(fullfile(chandir, PriorOnly_SigChans_fn));
                         
        if  isempty(PriorOnly_files)

            % generate null output files if there are no sig. channels for PriorOnly Condition
            % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

            save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band,Behavior));
            fid4 = fopen(save_file_name, 'w');
            fclose('all');

        else % get the labels for the PFC channels and AC channels 
       
             
               PriorOnly_SigChans = load(fullfile(chandir, PriorOnly_SigChans_fn), 'significant_channels');
               PriorOnly_modifiedCellArray   = regexprep(PriorOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
            
               % Separate elements that start with *PFC and *AC
               PriorOnly_PFC_chans  =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*PFC'));
               PriorOnly_AC_chans   =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*AC'));
     
               % Fetch the OnlyPretone Channels 
              
               PretoneOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, 'OnlyPretone', Behavior, Frequency_Band);
               PretoneOnly_files = dir(fullfile(chandir, PretoneOnly_SigChans_fn));

         
           if isempty(PretoneOnly_files)


              % generate null output files if there are no sig. channels for PretoneOnly Condition
              % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty
                                                
              save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
              fid4 = fopen(save_file_name, 'w');
              fclose('all');
                                                    
           else

             PretoneOnly_SigChans = load(fullfile(chandir, PretoneOnly_SigChans_fn), 'significant_channels');
             PretoneOnly_modifiedCellArray   = regexprep(PretoneOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
                        
             % Separate elements that start with *PFC and *AC
             PretoneOnly_PFC_chans  =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*PFC'));
             PretoneOnly_AC_chans   =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*AC'));
                    
             % Get the Channel Pairs that are shared by both the OnlyPrior and OnlyPretone Condition 
             % now we build a matrix of all of the shared channel pairs 

             % Find common elements
             common_PFC_chans = intersect(PriorOnly_PFC_chans, PretoneOnly_PFC_chans);
             common_AC_chans  = intersect(PriorOnly_AC_chans, PretoneOnly_AC_chans);

             % Get all pairwise combinations
             [PFC_comb, AC_comb] = ndgrid(common_PFC_chans, common_AC_chans);
            
             % Reshape into Nx2 cell array
             Shared_ChannelPairs = [PFC_comb(:), AC_comb(:)];
            
             % generate null output files if there are no sig. channels for PretoneOnly Condition
             % the code could still carry on even if both conditions had sig channels, but should be stopped if there are zero shared channels pairs

                if isempty(Shared_ChannelPairs)

                    save_file_name = fullfile(savedir, sprintf('NOSHAREDPAIRS_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band,Behavior));
                    fid4 = fopen(save_file_name, 'w');
                    fclose('all');

                else

                    % clear variables used to generate shared channel pairs
                    clear_list = {'AC_comb', 'common_AC_chans', 'common_PFC_chans', 'files', 'PFC_comb', 'PriorOnly_PFC_chans', 'PretoneOnly_PFC_chans', ...
                                  'PriorOnly_AC_chans', 'PretoneOnly_AC_chans', 'PriorOnly_modifiedCellArray', 'PretoneOnly_modifiedCellArray', ...
                                  'PriorOnly_SigChans_fn', 'PretoneOnly_SigChans_fn', 'PriorOnly_SigChans', 'PretoneOnly_SigChans'};
                    clear(clear_list{:});
                    clear clear_list;


                    %loop 3: for each shared channel pair 

                    for cp = 1:length(Shared_ChannelPairs)      % this loop skips if Shared_ChannelPairs is empty because length = 0

                        tic;
    
                        fprintf('\nNow running channel pair %s out of %s\n', num2str(cp), num2str(length(Shared_ChannelPairs)));
    
                        Current_ChanPair = Shared_ChannelPairs(cp,:);
                        Current_ChanPair_lb = strcat(Current_ChanPair{1}, '_', Current_ChanPair{2});
                        Current_ChanPair_lb = strrep(Current_ChanPair_lb, '*', '');

                        % calculate xcorr reults 
    
                            if strcmp(Behavior,'Correct') == 1
        
                            [OnlyPrior_shuffled_xcorr_result, OnlyPrior_shuffled_lagsResults,OnlyPrior_xcorr_result, OnlyPrior_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, 'OnlyPrior', Current_ChanPair);
                            [OnlyPretone_shuffled_xcorr_result, OnlyPretone_shuffled_lagsResults, OnlyPretone_xcorr_result, OnlyPretone_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, 'OnlyPretone', Current_ChanPair);
        
                            save_file_name = sprintf('%s_%s_%s_%s_%s_%s_XCorr_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            save(fullfile(savedir, save_file_name), 'OnlyPrior_shuffled_xcorr_result', 'OnlyPrior_shuffled_lagsResults',        ...
                                 'OnlyPrior_xcorr_result','OnlyPrior_lagsResults','OnlyPretone_shuffled_xcorr_result','OnlyPretone_shuffled_lagsResults','OnlyPretone_xcorr_result', ...
                                  'OnlyPretone_lagsResults');
        
                            % drop xcorr
        
                            clear_list = {'OnlyPrior_shuffled_xcorr_result', 'OnlyPrior_shuffled_lagsResults', 'OnlyPrior_xcorr_result', 'OnlyPrior_lagsResults', ...
                                          'OnlyPretone_shuffled_xcorr_result','OnlyPretone_shuffled_lagsResults','OnlyPretone_xcorr_result', 'OnlyPretone_lagsResults', ...
                                          'Fig1_name'};
                            clear(clear_list{:});
                            clear clear_list;
                            
                            end
            
                     end 

                            % Run Maris statistical test on Granger and Coherence
                            RunMaris_v3(Animal, RecDate, Epoch, Shared_ChannelPairs, Frequency_Band, Behavior, datadir, savedir, sessions, RandomIteration, threshold_value);
                 end          
           end 

        end
    
    end
                 
end
