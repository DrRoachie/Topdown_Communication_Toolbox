%% Script Summary

% Takes the pre-processed Taku data, and for each experimental condition
% (PriorOnly and PretoneOnly) calculates the coherence and granger spectra
% and tests if they are significantly different from one another. 



% Establish what sessions and conditions you want to process, we piecewise
% run the 5-6 sessions at a time to prevent critical time loss during
% processing. 

session_info = {'MrCassius', '190330' 'testToneOnset', 'Wrong', {'theta'};
                'MrCassius', '190404' 'testToneOnset', 'Wrong', {'theta'};
                'MrCassius', '190413' 'testToneOnset', 'Wrong', {'theta'};
                'MrCassius', '190416' 'testToneOnset', 'Wrong', {'theta'};
                'MrCassius', '190418' 'testToneOnset', 'Wrong', {'theta'}};


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
      
         tic;
        
         % loop 2: for every frequency band run the following 
         for fb = 1:length(freq_band_list)
        
               Frequency_Band = freq_band_list{fb};
        
                  % Fetch the OnlyPrior Channels (output of ChanWise_SpecEval)
            
                    Condition = 'OnlyPrior';

                    PriorOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, 'OnlyPrior', Behavior, Frequency_Band);
                                files = dir(fullfile(chandir, PriorOnly_SigChans_fn));
                                    
                                if ~isempty(files)
        
                                PriorOnly_SigChans = load(fullfile(chandir, PriorOnly_SigChans_fn), 'significant_channels');
                                PriorOnly_modifiedCellArray   = regexprep(PriorOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
        
                                % Separate elements that start with *PFC and *AC
                                PriorOnly_PFC_chans  =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*PFC'));
                                PriorOnly_AC_chans   =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*AC'));
                                
                                end
        
                                % generate null output files if there are no sig. channels for PriorOnly Condition
                                % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty

                                if isempty(files)
                                                Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_coh.txt', Animal, RecDate, Epoch, Frequency_Band, Behavior));
                                                Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_zthresh.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band, Behavior));
                                                Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_montecarlo.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band, Behavior));
                                                save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band,Behavior));
                                                fid1 = fopen(Fig1_name, 'w');
                                                fid2 = fopen(Fig2_name, 'w');
                                                fid3 = fopen(Fig3_name, 'w');
                                                fid4 = fopen(save_file_name, 'w');
                                                fclose('all');
                                else
                   
        
                  % Fetch the OnlyPretone Channels 
        
                    Condition = 'OnlyPretone';
                    
                    PretoneOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, 'OnlyPretone', Behavior, Frequency_Band);
                                files = dir(fullfile(chandir, PretoneOnly_SigChans_fn));
                                    
                                if ~isempty(files)
                                PretoneOnly_SigChans = load(fullfile(chandir, PretoneOnly_SigChans_fn), 'significant_channels');
                                PretoneOnly_modifiedCellArray   = regexprep(PretoneOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
        
                                % Separate elements that start with *PFC and *AC
                                PretoneOnly_PFC_chans  =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*PFC'));
                                PretoneOnly_AC_chans   =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*AC'));
                                end
        
                                % generate null output files if there are no sig. channels for PretoneOnly Condition
                                % this is supposed to stop the code because there is no point continuing if there either of the conditions are empty
                                if isempty(files)
                                                Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_coh.txt', Animal, RecDate, Epoch, Frequency_Band, Behavior));
                                                Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_zthresh.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
                                                Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_montecarlo.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
                                                save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
                                                fid1 = fopen(Fig1_name, 'w');
                                                fid2 = fopen(Fig2_name, 'w');
                                                fid3 = fopen(Fig3_name, 'w');
                                                fid4 = fopen(save_file_name, 'w');
                                                fclose('all');
                                else
        
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
                                                Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_coh.txt', Animal, RecDate, Epoch, Frequency_Band, Behavior));
                                                Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_zthresh.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band, Behavior));
                                                Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_montecarlo.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band, Behavior));
                                                save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band,Behavior));
                                                fid1 = fopen(Fig1_name, 'w');
                                                fid2 = fopen(Fig2_name, 'w');
                                                fid3 = fopen(Fig3_name, 'w');
                                                fid4 = fopen(save_file_name, 'w');
                                                fclose('all');
                     else
        
                    % clear variables used to generate shared channel pairs
                    clear_list = {'AC_comb', 'common_AC_chans', 'common_PFC_chans', 'files', 'PFC_comb', 'PriorOnly_PFC_chans', 'PretoneOnly_PFC_chans', ...
                                  'PriorOnly_AC_chans', 'PretoneOnly_AC_chans', 'PriorOnly_modifiedCellArray', 'PretoneOnly_modifiedCellArray', ...
                                  'PriorOnly_SigChans_fn', 'PretoneOnly_SigChans_fn', 'PriorOnly_SigChans', 'PretoneOnly_SigChans'};
                    clear(clear_list{:});
                    clear clear_list;
        
                    % load and filter data for PriorOnly and PretoneOnly
                    fName = sprintf('%s-%s_bdLFP_%s_ft.mat', Animal, RecDate, Epoch);
                    V = load(fullfile(datadir,fName));
        
                    params.choice           = V.choice;
                    params.err              = V.err;
                    params.pretone          = V.pretone;
                    params.pretoneLength    = V.pretoneLength;
                    params.prior            = cell2char(V.prior); 
                    params.SNR              = V.SNR;
                    
                    % separate data into right and wrong trials 
                    
                    iSelect = setStimulusCondition('OnlyPrior');    
        
                        if strcmp(Behavior,'Correct') == 1
                            iSelect.err = 'c'; % choose correct trials
                        elseif strcmp(Behavior,'Wrong')
                            iSelect.err = 'w';
                        end
        
                    data_prior = selectData(V.data,params,iSelect);
        
                    iSelect = setStimulusCondition('OnlyPretone');              
                    
                        if strcmp(Behavior,'Correct') == 1
                            iSelect.err = 'c'; % choose correct trials
                        elseif strcmp(Behavior,'Wrong')
                            iSelect.err = 'w';
                        end
        
                    data_pretone = selectData(V.data,params,iSelect);
        
                    % compute cross spectral density matrix for PriorOnly and PretoneOnly conditions
        
                    cfg            = [];
                    cfg.method     = 'mtmfft';
                    cfg.taper      = 'dpss';
                    cfg.output     = 'fourier';
                    cfg.tapsmofrq  = 4; 
                    cfg.pad        = 1;
                    cfg.foilim     = [0 100];
                    
                    freq_prior     = ft_freqanalysis(cfg, data_prior);
                    freq_pretone   = ft_freqanalysis(cfg, data_pretone);
        
                    clear_list = {'fName', 'V', 'params', 'iSelect', 'cfg'};
                    clear(clear_list{:});
                    clear clear_list;
        
                          % Spectrogram to show the LFP spectra of the current channel pair 
        
                          if strcmp(Behavior,'Correct') == 1
                            [tfreq_OnlyPrior_c]   = GetSpectrogram(specdir, Animal, RecDate, Epoch, 'OnlyPrior', Behavior);
                            [tfreq_OnlyPretone_c] = GetSpectrogram(specdir, Animal, RecDate, Epoch, 'OnlyPretone', Behavior);
                          end
                    
                          if strcmp(Behavior,'Wrong') == 1
                            [tfreq_OnlyPrior_w]   = GetSpectrogram(specdir, Animal, RecDate, Epoch, 'OnlyPrior', Behavior);
                            [tfreq_OnlyPretone_w] = GetSpectrogram(specdir, Animal, RecDate, Epoch, 'OnlyPretone', Behavior);
                          end

                            % for every channel pair calculate and save the meta data you need to describe the communication between PFC and AC
                            for cp = 1:length(Shared_ChannelPairs)      % this loop skips if Shared_ChannelPairs is empty because length = 0

                            tic;

                            fprintf('\nNow running channel pair %s out of %s\n', num2str(cp), num2str(length(Shared_ChannelPairs)));

                            Current_ChanPair = Shared_ChannelPairs(cp,:);
                            Current_ChanPair_lb = strcat(Current_ChanPair{1}, '_', Current_ChanPair{2});
                            Current_ChanPair_lb = strrep(Current_ChanPair_lb, '*', '');

                            % plot the spectrograms and save
                            if strcmp(Behavior,'Correct') == 1
                            plot_spectrogram_ChanPair(tfreq_OnlyPrior_c, Current_ChanPair, Frequency_Band)
                            plot_spectrogram_ChanPair(tfreq_OnlyPretone_c, Current_ChanPair, Frequency_Band)
                            end 

                            if strcmp(Behavior,'Wrong') == 1
                            plot_spectrogram_ChanPair(tfreq_OnlyPrior_w, Current_ChanPair, Frequency_Band)
                            plot_spectrogram_ChanPair(tfreq_OnlyPretone_w, Current_ChanPair, Frequency_Band)
                            end 

                            Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Coh_Prior_Spectrogram.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            Fig2_name      = sprintf('%s_%s_%s_%s_%s_%s_Coh_Pretone_Spectrogram.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            saveas(1, fullfile(savedir, Fig1_name));
                            saveas(2, fullfile(savedir, Fig2_name));
                            close all


                            % Calculate the coherence spectra of PriorOnly and PretoneOnly Trials 

                            if strcmp(Behavior,'Correct') == 1
                           
                            [PFC_AC_Coh_Spectrum_PriorOnly_c, ~, ~]                        = getDirectoryCoherence(sessions, freq_prior, data_prior, Current_ChanPair, Behavior);    % only take coh spectrum, ignore other returned values (unused in this script)
                            [PFC_AC_Coh_Spectrum_PretoneOnly_c, ~, ~]                      = getDirectoryCoherence(sessions, freq_pretone, data_pretone, Current_ChanPair, Behavior);
                            [descriptives_onlyprior_coh_c, descriptives_onlypretone_coh_c] = GetCohDescriptives(PFC_AC_Coh_Spectrum_PriorOnly_c, PFC_AC_Coh_Spectrum_PretoneOnly_c); % get the coherence descriptions ready
                           
                            Plot_CoherenceSpec(Frequency_Band, descriptives_onlyprior_coh_c, descriptives_onlypretone_coh_c)                               % plot the coherence spectra of PriorOnly and PretoneOnly

                            Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Coh_Spec.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

                            saveas(1, fullfile(savedir, Fig1_name));
                            close all

                            save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Coh_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            save(fullfile(savedir, save_file_name), 'tfreq_OnlyPrior_c','tfreq_OnlyPretone_c');

                            end

                            if strcmp(Behavior,'Wrong') == 1
                           
                            [PFC_AC_Coh_Spectrum_PriorOnly_w, ~, ~]                        = getDirectoryCoherence(sessions, freq_prior, data_prior, Current_ChanPair, Behavior);    % only take coh spectrum, ignore other returned values (unused in this script)
                            [PFC_AC_Coh_Spectrum_PretoneOnly_w, ~, ~]                      = getDirectoryCoherence(sessions, freq_pretone, data_pretone, Current_ChanPair, Behavior);
                            [descriptives_onlyprior_coh_w, descriptives_onlypretone_coh_w] = GetCohDescriptives(PFC_AC_Coh_Spectrum_PriorOnly_w, PFC_AC_Coh_Spectrum_PretoneOnly_w); % get the coherence descriptions ready
                           
                            Plot_CoherenceSpec(Frequency_Band, descriptives_onlyprior_coh_w, descriptives_onlypretone_coh_w)                               % plot the coherence spectra of PriorOnly and PretoneOnly

                            Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Coh_Spec.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

                            saveas(1, fullfile(savedir, Fig1_name));
                            close all

                            save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Coh_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            save(fullfile(savedir, save_file_name), 'tfreq_OnlyPrior_w','tfreq_OnlyPretone_w');
                            
                            end

                            % clear tfreq, the coh spectra, the results of the coherence Maris test

                            clear_list = {'tfreq_OnlyPrior_c', 'tfreq_OnlyPretone_c', 'Fig1_name', 'Fig2_name', 'Fig3_name', 'PFC_AC_Coh_Spectrum_PriorOnly_c', 'PFC_AC_Coh_Spectrum_PretoneOnly_c', ...
                                          'descriptives_onlyprior_coh_c', 'descriptives_onlypretone_coh_c', 'tfreq_OnlyPrior_w', 'tfreq_OnlyPretone_w', 'PFC_AC_Coh_Spectrum_PriorOnly_w', 'PFC_AC_Coh_Spectrum_PretoneOnly_w', ... 
                                          'descriptives_onlyprior_coh_w', 'descriptives_onlypretone_coh_w'};


                            clear(clear_list{:});
                            clear clear_list;
                                                      

                           % Calculate the granger spectra of PriorOnly and PretoneOnly Trials 
 

                            if strcmp(Behavior,'Correct') == 1

                            Send_Chan_lb        = Current_ChanPair{1};
                            Rec_Chan_lb         = Current_ChanPair{2};

                            [Granger_PriorOnly_c, STD_Granger_PriorOnly_c]      = getGrangerSpec(freq_prior, data_prior, Send_Chan_lb, Rec_Chan_lb, Frequency_Band);        % grangerspctrm(1,2,:) is sender --> receiver
                            [Granger_PretoneOnly_c, STD_Granger_PretoneOnly_c]  = getGrangerSpec(freq_pretone, data_pretone, Send_Chan_lb, Rec_Chan_lb, Frequency_Band);    % grangerspctrm(2,1,:) is receiver --> sender

                            % Compute Null Threshold
                            Bstrp_Iteration         = 100;

                            Null_Thresh_PriorOnly_c   = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, 'PriorOnly', Frequency_Band, Animal, RecDate, Epoch, datadir);
                            Null_Thresh_PretoneOnly_c = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, 'PretoneOnly', Frequency_Band, Animal, RecDate, Epoch, datadir);
                
                            % Plot
                            Plot_Granger(Send_Chan_lb, Rec_Chan_lb, Granger_PriorOnly_c, Granger_PretoneOnly_c, STD_Granger_PriorOnly_c, STD_Granger_PretoneOnly_c, ...
                            Null_Thresh_PriorOnly_c, Null_Thresh_PretoneOnly_c);

                            % save data/figures
                            Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_Spec_PFC_to_AC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            Fig2_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_Spec_AC_to_PFC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

                            saveas(1, fullfile(savedir, Fig1_name));
                            saveas(2, fullfile(savedir, Fig2_name));
                            save(fullfile(savedir, save_file_name), 'Granger_PriorOnly_c', 'Granger_PretoneOnly_c', ...
                                 'STD_Granger_PriorOnly_c','STD_Granger_PretoneOnly_c','Null_Thresh_PriorOnly_c','Null_Thresh_PretoneOnly_c');

                            % drop granger 

                            clear_list = {'Send_Chan_lb', 'Rec_Chan_lb', 'Granger_PriorOnly_c', 'Granger_PretoneOnly_c', ...
                                 'STD_Granger_PriorOnly_c','STD_Granger_PretoneOnly_c','Null_Thresh_PriorOnly_c','Null_Thresh_PretoneOnly_c', 'Bstrp_Iteration', ...
                                 'Fig1_name', 'Fig2_name'};
                            clear(clear_list{:});
                            clear clear_list;

                            close all;

                            end

                            if strcmp(Behavior,'Wrong') == 1

                            Send_Chan_lb        = Current_ChanPair{1};
                            Rec_Chan_lb         = Current_ChanPair{2};

                            [Granger_PriorOnly_w, STD_Granger_PriorOnly_w]      = getGrangerSpec(freq_prior, data_prior, Send_Chan_lb, Rec_Chan_lb, Frequency_Band);        % grangerspctrm(1,2,:) is sender --> receiver
                            [Granger_PretoneOnly_w, STD_Granger_PretoneOnly_w]  = getGrangerSpec(freq_pretone, data_pretone, Send_Chan_lb, Rec_Chan_lb, Frequency_Band);    % grangerspctrm(2,1,:) is receiver --> sender

                            % Compute Null Threshold
                            Bstrp_Iteration         = 100;

                            Null_Thresh_PriorOnly_w   = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, 'PriorOnly', Frequency_Band, Animal, RecDate, Epoch, datadir);
                            Null_Thresh_PretoneOnly_w = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, 'PretoneOnly', Frequency_Band, Animal, RecDate, Epoch, datadir);
                
                            % Plot
                            Plot_Granger(Send_Chan_lb, Rec_Chan_lb, Granger_PriorOnly_w, Granger_PretoneOnly_w, STD_Granger_PriorOnly_w, STD_Granger_PretoneOnly_w, ...
                            Null_Thresh_PriorOnly_w, Null_Thresh_PretoneOnly_w);

                            % save data/figures
                            Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_Spec_PFC_to_AC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            Fig2_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_Spec_AC_to_PFC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

                            saveas(1, fullfile(savedir, Fig1_name));
                            saveas(2, fullfile(savedir, Fig2_name));
                            save(fullfile(savedir, save_file_name), 'Granger_PriorOnly_w', 'Granger_PretoneOnly_w', ...
                                 'STD_Granger_PriorOnly_w','STD_Granger_PretoneOnly_w','Null_Thresh_PriorOnly_w','Null_Thresh_PretoneOnly_w');

                            % drop granger 

                            clear_list = {'Send_Chan_lb', 'Rec_Chan_lb', 'Granger_PriorOnly_w', 'Granger_PretoneOnly_w', ...
                                 'STD_Granger_PriorOnly_w','STD_Granger_PretoneOnly_w','Null_Thresh_PriorOnly_w','Null_Thresh_PretoneOnly_w', 'Bstrp_Iteration', ...
                                 'Fig1_name', 'Fig2_name'};
                            clear(clear_list{:});
                            clear clear_list;

                            close all;

                            end

                            % calculate xcorr reults 

                            if strcmp(Behavior,'Correct') == 1

                            [OnlyPrior_shuffled_xcorr_result, OnlyPrior_shuffled_lagsResults,OnlyPrior_xcorr_result, OnlyPrior_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, 'OnlyPrior', Current_ChanPair);
                            [OnlyPretone_shuffled_xcorr_result, OnlyPretone_shuffled_lagsResults, OnlyPretone_xcorr_result, OnlyPretone_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, 'OnlyPretone', Current_ChanPair);

                            %  plot xcorr results  
                            figure;
                            hold on;
                            plot(OnlyPrior_lagsResults, OnlyPrior_xcorr_result.avg, 'g', 'LineWidth', 1.5);
                            fill([OnlyPrior_lagsResults, fliplr(OnlyPrior_lagsResults)], [OnlyPrior_xcorr_result.upper, fliplr(OnlyPrior_xcorr_result.lower)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                            plot(OnlyPrior_shuffled_lagsResults, OnlyPrior_shuffled_xcorr_result.avg, 'k', 'LineWidth', 1.5);
                            fill([OnlyPrior_shuffled_lagsResults, fliplr(OnlyPrior_lagsResults)], [OnlyPrior_shuffled_xcorr_result.upper, fliplr(OnlyPrior_shuffled_xcorr_result.lower)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                            plot(OnlyPretone_lagsResults, OnlyPretone_xcorr_result.avg, 'b', 'LineWidth', 1.5);
                            fill([OnlyPretone_lagsResults, fliplr(OnlyPretone_lagsResults)], [OnlyPretone_xcorr_result.upper, fliplr(OnlyPretone_xcorr_result.lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                            plot(OnlyPretone_shuffled_lagsResults, OnlyPretone_shuffled_xcorr_result.avg, 'r', 'LineWidth', 1.5);
                            fill([OnlyPretone_shuffled_lagsResults, fliplr(OnlyPretone_shuffled_lagsResults)], [OnlyPretone_shuffled_xcorr_result.upper, fliplr(OnlyPretone_shuffled_xcorr_result.lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                            xlabel('Lag');
                            ylabel('Average Trial-by-Trial Cross-Correlation');
                            title('PFC-to-AC Cross-Correlation with 95% CI');
                            %legend('LED Only', '' ,'NULL','', '','','Pretone Only','')
                            hold off;

                            % save figure

                            Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_XCorr.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                            saveas(1, fullfile(savedir, Fig1_name));

                            close all

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

                            elapsedTime = toc;
                            disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);

                            end
                            end

                            
                            % run Maris statistical test on Granger and Coherence

                            if length(Shared_ChannelPairs) ~= 0    % only run Maris statistical test if Shared_ChannelPairs is not empty
                                RunMaris_v3(Animal, RecDate, Epoch, Shared_ChannelPairs, Frequency_Band, Behavior, datadir, savedir, sessions, RandomIteration, threshold_value);
                            end
                        
                        end      
                        end
                        end
                 
end

end
elapsedTime = toc;
disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);