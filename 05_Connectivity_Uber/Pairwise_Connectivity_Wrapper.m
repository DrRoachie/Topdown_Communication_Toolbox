%% Script Summary

% Takes the pre-processed Taku data and calculates coherence. This
% particular file takes the coherence and averages across 4 groups of
% layers.  


%% Set directory and get a list of all files in the folder with the desired file name pattern.

datadir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2\MrCassius\testTone\190418'; % The cut epoc data 
chandir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_05_25_Analysis\00_Significant_chans\MrCassius\testTone\190418'; % the sig channels that are output
specdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed\MrCassius\testToneonset'; % this is where the uncut preprocessed data is stored; you need the full time series for the spectogram
savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_06_10_Analysis';
sessions = dir(fullfile(datadir,'*.mat'));
addpath(genpath(datadir));
addpath(genpath(chandir));
isSave = 0;  

%% Establish the subset of data that you want to process. 
Animal           = 'MrCassius';                  % Options: 'MrCassius', 'MrM'; 
RecDate          = '190418';        
Epoch            = 'testToneOnset';                % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
Behavior         = 'Correct';                    % Options:  'Correct', Wrong, due to memory constraints                          
RandomIteration  = 300;                        % options: the number of times you want to allocate the random partition, and generate the test statistic
%freq_band_list   = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};
freq_band_list   = {'beta'};
threshold_value  = 0.5;       % the threshold used for the monte carlo test between PriorOnly and PretoneOnly Conditions

%% for every frequency band run 

 for fb = 1:length(freq_band_list)

       Frequency_Band = freq_band_list{fb};

          % Fetch the OnlyPrior Channels (output of ChanWise_SpecEval)
    
            Condition = 'OnlyPrior';
            
            PriorOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
                        files = dir(fullfile(chandir, PriorOnly_SigChans_fn));
                            
                        if ~isempty(files)

                        PriorOnly_SigChans = load(fullfile(chandir, PriorOnly_SigChans_fn), 'significant_channels');
                        PriorOnly_modifiedCellArray   = regexprep(PriorOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

                        % Separate elements that start with *PFC and *AC
                        PriorOnly_PFC_chans  =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*PFC'));
                        PriorOnly_AC_chans   =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*AC'));
                        
                        end

             
                        if isempty(files)
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_zthresh.txt', RecDate, Frequency_Band, Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_montecarlo.txt', RecDate, Frequency_Band, Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_data.txt', RecDate, Frequency_Band,Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose(fid1);
                                        fclose(fid2);
                                        fclose(fid3);
                                        fclose(fid4);
                        else
           

          % Fetch the OnlyPretone Channels 

            Condition = 'OnlyPretone';
            
            PretoneOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
                        files = dir(fullfile(chandir, PretoneOnly_SigChans_fn));
                            
                        if ~isempty(files)
                        PretoneOnly_SigChans = load(fullfile(chandir, PretoneOnly_SigChans_fn), 'significant_channels');
                        PretoneOnly_modifiedCellArray   = regexprep(PretoneOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

                        % Separate elements that start with *PFC and *AC
                        PretoneOnly_PFC_chans  =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*PFC'));
                        PretoneOnly_AC_chans   =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*AC'));
                        end

                        if isempty(files)
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_zthresh.txt', RecDate, Frequency_Band,  Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_montecarlo.txt', RecDate, Frequency_Band,  Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_data.txt', RecDate, Frequency_Band,  Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose(fid1);
                                        fclose(fid2);
                                        fclose(fid3);
                                        fclose(fid4);
                        else

            % Get the Channel Pairs that are shared by both the OnlyPrior and OnlyPretone Condition 

            % Find common elements
            common_PFC_chans = intersect(PriorOnly_PFC_chans, PretoneOnly_PFC_chans);
            common_AC_chans = intersect(PriorOnly_AC_chans, PretoneOnly_AC_chans);

            % Get all pairwise combinations
            [PFC_comb, AC_comb] = ndgrid(common_PFC_chans, common_AC_chans);
            
            % Reshape into Nx2 cell array
            Shared_ChannelPairs = [PFC_comb(:), AC_comb(:)];

                            % for every channel pair calculate and save the meta data you need to describe the communication between PFC and AC
                            
                            for cp = 1:length(Shared_ChannelPairs)

                            Current_ChanPair = Shared_ChannelPairs(cp,:);
                            Current_ChanPair_lb = strcat(Current_ChanPair{1}, '_', Current_ChanPair{2});
                            Current_ChanPair_lb = strrep(Current_ChanPair_lb, '*', '');

                            % calculate xcorr reults 

                             Condition = 'OnlyPrior';
                            [OnlyPrior_shuffled_xcorr_result, OnlyPrior_shuffled_lagsResults,OnlyPrior_xcorr_result, OnlyPrior_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, Condition, Current_ChanPair);
                            
                             Condition = 'OnlyPretone';
                            [OnlyPretone_shuffled_xcorr_result, OnlyPretone_shuffled_lagsResults, OnlyPretone_xcorr_result, OnlyPretone_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, Condition, Current_ChanPair);
                            
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

                            Fig6_name      = sprintf('xcorr_%s_%s_%s_%s_%s_xcorr.fig', RecDate, Frequency_Band, Epoch, Current_ChanPair_lb, Behavior);
                            saveas(1, fullfile(savedir, Fig6_name));

                            close

                            end
                        end                        
                        end
 end