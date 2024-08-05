%% Script Summary

% Takes the pre-processed Taku data and calculates coherence. This
% particular file takes the coherence and averages across 4 groups of
% layers. 


%% Set directory and get a list of all files in the folder with the desired file name pattern.

isSave = 0;
plot   = 1;
save_file_name = 'Statistic_Test_Results.mat';
mainDirectory = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\04_Epoc_Cut\MrCassius\LED';
savedir = 'C:\Users\Corey Roach\Desktop\2024_01_18_Analysis\MrCassius\PFCdeep_ACupmid'; % Options: end of directory should be one of the four chan clusters
sessions = dir(mainDirectory);
sessions = sessions([sessions.isdir] & ~ismember({sessions.name}, {'.', '..'}));

for s = 1:length(sessions)
    
    currentsession = fullfile(mainDirectory, sessions(s).name);
    matFile        = dir(fullfile(currentsession, '*.mat'));
    addpath(genpath(currentsession));
    

            % Establish the subset of data that you want to process. 
            
            Animal          = 'MrCassius';                  
            Epoch           = 'preCueOnset';              % options: 'preCueOnset', 'moveOnset', & 'testToneOnset';
            Channels        = 'PFCdeep_ACupmid';           % options: 'allChan', 'PFClowmid_ACupmid','PFCdeep_ACupmid', 'PFCupper_ACupper', or 'PFCdeep_ACdeep';
            RandomIteration = 300;                        % options: the number of times you want to allocate the random partition, and generate the test statistic;
            Behavior        = 'Both';                     % options: 'Both', 'Correct', Wrong;
            
            % Get PriorOnly Coherence 
            
            Condition = 'OnlyPrior';                          % options: 'OnlyPrior', 'OnlyPretone', & 'Both';
            [PFC_AC_Coh_Spectrum_PriorOnly_c,PFC_AC_Coh_Spectrum_PriorOnly_w,PFC_AC_Coh_fd_PriorOnly_c, PFC_AC_Coh_fd_PriorOnly_w, TrialInfo_PriorOnly_c,TrialInfo_PriorOnly_w] = getDirectoryCoherence(matFile, Condition, Channels);
            
            
            % Get PretoneOnly Coherence 
            
            Condition = 'OnlyPretone';                          % options: 'OnlyPrior', 'OnlyPretone', & 'Both';
            [PFC_AC_Coh_Spectrum_PretoneOnly_c,PFC_AC_Coh_Spectrum_PretoneOnly_w,PFC_AC_Coh_fd_PretoneOnly_c, PFC_AC_Coh_fd_PretoneOnly_w, TrialInfo_PretoneOnly_c,TrialInfo_PretoneOnly_w] = getDirectoryCoherence(matFile, Condition, Channels);
            
            % Use the coherence spectra of the prior and pretone conditions to calculate the Maris z-statistic 
            
            data_Maris_zstatistic_C = []; 
            data_Maris_zstatistic_W = []; 
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1 
                
                for k = 1:length(matFile)
                Coh_test_statistic_a_C = atanh(abs(PFC_AC_Coh_Spectrum_PriorOnly_c(k).cohspctrm)) - (1/(PFC_AC_Coh_Spectrum_PriorOnly_c(k).dof-2)); % derived from Maris et al 2007 Eq 1 
                Coh_test_statistic_b_C = atanh(abs(PFC_AC_Coh_Spectrum_PretoneOnly_c(k).cohspctrm)) - (1/(PFC_AC_Coh_Spectrum_PretoneOnly_c(k).dof-2));
                Coh_test_statistic_c_C = sqrt((1/(PFC_AC_Coh_Spectrum_PriorOnly_c(k).dof-2))+(1/(PFC_AC_Coh_Spectrum_PretoneOnly_c(k).dof-2)));
                Coh_test_statistic_Z_C = Coh_test_statistic_a_C-Coh_test_statistic_b_C/Coh_test_statistic_c_C;
                Coh_test_statistic_Zavg_C = mean( Coh_test_statistic_Z_C , 1);                                                                          % collapse across channel pairs 
                data_Maris_zstatistic_C = [data_Maris_zstatistic_C; Coh_test_statistic_Zavg_C];
                end
            end
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1 
                for k = 1:length(matFile)
                Coh_test_statistic_a_W = atanh(abs(PFC_AC_Coh_Spectrum_PriorOnly_w(k).cohspctrm)) - (1/(PFC_AC_Coh_Spectrum_PriorOnly_w(k).dof-2)); % derived from Maris et al 2007 Eq 1 
                Coh_test_statistic_b_W = atanh(abs(PFC_AC_Coh_Spectrum_PretoneOnly_w(k).cohspctrm)) - (1/(PFC_AC_Coh_Spectrum_PretoneOnly_w(k).dof-2));
                Coh_test_statistic_c_W = sqrt((1/(PFC_AC_Coh_Spectrum_PriorOnly_w(k).dof-2))+(1/(PFC_AC_Coh_Spectrum_PretoneOnly_w(k).dof-2)));
                Coh_test_statistic_Z_W = Coh_test_statistic_a_W-Coh_test_statistic_b_W/Coh_test_statistic_c_W;
                Coh_test_statistic_Zavg_W = mean( Coh_test_statistic_Z_W, 1);                                                                          % collapse across channel pairs 
                data_Maris_zstatistic_W = [data_Maris_zstatistic_W; Coh_test_statistic_Zavg_W];
                end
            end
                                                                                                    
            % Pool Prior and Pretone Trials 
            
            [RandomPartition_c_1_Size, RandomPartition_c_2_Size, RandomPartition_w_1_Size, RandomPartition_w_2_Size, All_Trials_c, All_Trials_w] = PoolPriorPretone(currentsession, matFile, Animal, Epoch);
            
            % Using the pooled data, for each permutation divide the data into two subsets and calculate the z-statistic between the two subsets. 
            
            Monte_Carlo_Array_Maris_C   = [];
            Monte_Carlo_Array_Maris_W   = [];
            
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
            
                for i = 1:RandomIteration
                
                % generate a random partition equal to the number of trials in the condition with more trials 
                
                    totalTrials_C = size(All_Trials_c.trial, 2);
                    randomCols_1_C = randperm(totalTrials_C, RandomPartition_c_1_Size); 
                    randomRows_1_C = randomCols_1_C';
                
                    for r = 1:RandomPartition_c_1_Size
                    Random_Partition_c_1.time(:, r)        = All_Trials_c.time(:, randomCols_1_C(r));
                    Random_Partition_c_1.trial(:, r)       = All_Trials_c.trial(:, randomCols_1_C(r));
                    Random_Partition_c_1.sampleinfo(r, :)  = All_Trials_c.sampleinfo(randomRows_1_C(r), :);
                    end
                    
                    Random_Partition_c_1.label      = All_Trials_c(1).label;
                    Random_Partition_c_1.fsample    = All_Trials_c(1).fsample;
                
                    remainingCols_C = setdiff(1:totalTrials_C, randomCols_1_C);
                    randomCols_2_C  = remainingCols_C(randperm(length(remainingCols_C)));
                    randomRows_2_C  = randomCols_2_C'; 
                
                    for rr = 1:RandomPartition_c_2_Size
                    Random_Partition_c_2.time(:, rr)        = All_Trials_c.time(:, randomCols_2_C(rr));
                    Random_Partition_c_2.trial(:, rr)       = All_Trials_c.trial(:, randomCols_2_C(rr));
                    Random_Partition_c_2.sampleinfo(rr, :)  = All_Trials_c.sampleinfo(randomRows_2_C(rr), :);
                    end
                
                    Random_Partition_c_2.label      = All_Trials_c(1).label;
                    Random_Partition_c_2.fsample    = All_Trials_c(1).fsample;
                
                        % non-parametric computation of the cross-spectral density matrix (slow)
                            
                            cfg            = [];
                            cfg.method     = 'mtmfft';
                            cfg.taper      = 'hanning';
                            cfg.output     = 'fourier';
                            %cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
                            cfg.pad        = 1;
                            cfg.foilim     = [1 100];
                            
                            freq_1_C        = ft_freqanalysis(cfg, Random_Partition_c_1);
                            freq_2_C        = ft_freqanalysis(cfg, Random_Partition_c_2);
                            
                                % calculate the coherence spectra in partition 1 and partition 2 
                    
                                if strcmp(Channels,'allChan') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb = {'*PFC*','*AC*'};
                            
                                coh_PFC_AC_1_C = ft_connectivityanalysis(cfg, freq_1_C);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_C.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_C);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_C
                        
                                coh_PFC_AC_2_C = ft_connectivityanalysis(cfg, freq_2_C);       % calculate the coherence spectra of partition 2 
                                clear freq_2_C
                               
                                C_1_C = reshape_coherence(coh_PFC_AC_1_C, eID_a, eID_b);            % reshape partition 1
                                C_2_C = reshape_coherence(coh_PFC_AC_2_C, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_C
                                clear coh_PFC_AC_2_C
                                end 
                                
                        
                                if strcmp(Channels,'PFClowmid_ACupmid') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch06' '*PFC_ch11'
                                                   '*AC_ch06' '*PFC_ch12'
                                                   '*AC_ch06' '*PFC_ch13'
                                                   '*AC_ch06' '*PFC_ch14'
                                                   '*AC_ch06' '*PFC_ch15'
                                                   '*AC_ch07' '*PFC_ch11'
                                                   '*AC_ch07' '*PFC_ch12'
                                                   '*AC_ch07' '*PFC_ch13'
                                                   '*AC_ch07' '*PFC_ch14'
                                                   '*AC_ch07' '*PFC_ch15'
                                                   '*AC_ch08' '*PFC_ch11'
                                                   '*AC_ch08' '*PFC_ch12'
                                                   '*AC_ch08' '*PFC_ch13'
                                                   '*AC_ch08' '*PFC_ch14'
                                                   '*AC_ch08' '*PFC_ch15'
                                                   '*AC_ch09' '*PFC_ch11'
                                                   '*AC_ch09' '*PFC_ch12'
                                                   '*AC_ch09' '*PFC_ch13'
                                                   '*AC_ch09' '*PFC_ch14'
                                                   '*AC_ch09' '*PFC_ch15'
                                                   '*AC_ch10' '*PFC_ch11'
                                                   '*AC_ch10' '*PFC_ch12'
                                                   '*AC_ch10' '*PFC_ch13'
                                                   '*AC_ch10' '*PFC_ch14'
                                                   '*AC_ch10' '*PFC_ch15'};   
                                               
                                coh_PFC_AC_1_C = ft_connectivityanalysis(cfg, freq_1_C);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_C.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_C);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_C
                        
                                coh_PFC_AC_2_C = ft_connectivityanalysis(cfg, freq_2_C);       % calculate the coherence spectra of partition 2 
                                clear freq_2_C
                               
                                C_1_C = reshape_coherence(coh_PFC_AC_1_C, eID_a, eID_b);            % reshape partition 1
                                C_2_C = reshape_coherence(coh_PFC_AC_2_C, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_C
                                clear coh_PFC_AC_2_C
                                end 
                            
                                if strcmp(Channels,'PFCdeep_ACupmid') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch06' '*PFC_ch16'
                                                   '*AC_ch06' '*PFC_ch17'
                                                   '*AC_ch06' '*PFC_ch18'
                                                   '*AC_ch06' '*PFC_ch19'
                                                   '*AC_ch06' '*PFC_ch20'
                                                   '*AC_ch07' '*PFC_ch16'
                                                   '*AC_ch07' '*PFC_ch17'
                                                   '*AC_ch07' '*PFC_ch18'
                                                   '*AC_ch07' '*PFC_ch19'
                                                   '*AC_ch07' '*PFC_ch20'
                                                   '*AC_ch08' '*PFC_ch16'
                                                   '*AC_ch08' '*PFC_ch17'
                                                   '*AC_ch08' '*PFC_ch18'
                                                   '*AC_ch08' '*PFC_ch19'
                                                   '*AC_ch08' '*PFC_ch20'
                                                   '*AC_ch09' '*PFC_ch16'
                                                   '*AC_ch09' '*PFC_ch17'
                                                   '*AC_ch09' '*PFC_ch18'
                                                   '*AC_ch09' '*PFC_ch19'
                                                   '*AC_ch09' '*PFC_ch20'
                                                   '*AC_ch10' '*PFC_ch16'
                                                   '*AC_ch10' '*PFC_ch17'
                                                   '*AC_ch10' '*PFC_ch18'
                                                   '*AC_ch10' '*PFC_ch19'
                                                   '*AC_ch10' '*PFC_ch20'};   
                                               
                                coh_PFC_AC_1_C = ft_connectivityanalysis(cfg, freq_1_C);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_C.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_C);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_C
                        
                                coh_PFC_AC_2_C = ft_connectivityanalysis(cfg, freq_2_C);       % calculate the coherence spectra of partition 2 
                                clear freq_2_C
                               
                                C_1_C = reshape_coherence(coh_PFC_AC_1_C, eID_a, eID_b);            % reshape partition 1
                                C_2_C = reshape_coherence(coh_PFC_AC_2_C, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_C
                                clear coh_PFC_AC_2_C
                                end 
                                
                            
                                if strcmp(Channels,'PFCdeep_ACdeep') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch16' '*PFC_ch16'
                                                   '*AC_ch16' '*PFC_ch17'
                                                   '*AC_ch16' '*PFC_ch18'
                                                   '*AC_ch16' '*PFC_ch19'
                                                   '*AC_ch16' '*PFC_ch20'
                                                   '*AC_ch17' '*PFC_ch16'
                                                   '*AC_ch17' '*PFC_ch17'
                                                   '*AC_ch17' '*PFC_ch18'
                                                   '*AC_ch17' '*PFC_ch19'
                                                   '*AC_ch17' '*PFC_ch20'
                                                   '*AC_ch18' '*PFC_ch16'
                                                   '*AC_ch18' '*PFC_ch17'
                                                   '*AC_ch18' '*PFC_ch18'
                                                   '*AC_ch18' '*PFC_ch19'
                                                   '*AC_ch18' '*PFC_ch20'
                                                   '*AC_ch19' '*PFC_ch16'
                                                   '*AC_ch19' '*PFC_ch17'
                                                   '*AC_ch19' '*PFC_ch18'
                                                   '*AC_ch19' '*PFC_ch19'
                                                   '*AC_ch19' '*PFC_ch20'
                                                   '*AC_ch20' '*PFC_ch16'
                                                   '*AC_ch20' '*PFC_ch17'
                                                   '*AC_ch20' '*PFC_ch18'
                                                   '*AC_ch20' '*PFC_ch19'
                                                   '*AC_ch20' '*PFC_ch20'};   
                                               
                                coh_PFC_AC_1_C = ft_connectivityanalysis(cfg, freq_1_C);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_C.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_C);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_C
                        
                                coh_PFC_AC_2_C = ft_connectivityanalysis(cfg, freq_2_C);       % calculate the coherence spectra of partition 2 
                                clear freq_2_C
                               
                                C_1_C = reshape_coherence(coh_PFC_AC_1_C, eID_a, eID_b);            % reshape partition 1
                                C_2_C = reshape_coherence(coh_PFC_AC_2_C, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_C
                                clear coh_PFC_AC_2_C
                                end 
                        
                                if strcmp(Channels,'PFCupper_ACupper') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch03' '*PFC_ch03'
                                                   '*AC_ch03' '*PFC_ch04'
                                                   '*AC_ch03' '*PFC_ch05'
                                                   '*AC_ch03' '*PFC_ch06'
                                                   '*AC_ch03' '*PFC_ch07'
                                                   '*AC_ch04' '*PFC_ch03'
                                                   '*AC_ch04' '*PFC_ch04'
                                                   '*AC_ch04' '*PFC_ch05'
                                                   '*AC_ch04' '*PFC_ch06'
                                                   '*AC_ch04' '*PFC_ch07'
                                                   '*AC_ch05' '*PFC_ch03'
                                                   '*AC_ch05' '*PFC_ch04'
                                                   '*AC_ch05' '*PFC_ch05'
                                                   '*AC_ch05' '*PFC_ch06'
                                                   '*AC_ch05' '*PFC_ch07'
                                                   '*AC_ch06' '*PFC_ch03'
                                                   '*AC_ch06' '*PFC_ch04'
                                                   '*AC_ch06' '*PFC_ch05'
                                                   '*AC_ch06' '*PFC_ch06'
                                                   '*AC_ch06' '*PFC_ch07'
                                                   '*AC_ch07' '*PFC_ch03'
                                                   '*AC_ch07' '*PFC_ch04'
                                                   '*AC_ch07' '*PFC_ch05'
                                                   '*AC_ch07' '*PFC_ch06'
                                                   '*AC_ch07' '*PFC_ch07'};   
                                               
                                coh_PFC_AC_1_C = ft_connectivityanalysis(cfg, freq_1_C);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_C.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_C);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_C
                        
                                coh_PFC_AC_2_C = ft_connectivityanalysis(cfg, freq_2_C);       % calculate the coherence spectra of partition 2 
                                clear freq_2_C
                               
                                C_1_C = reshape_coherence(coh_PFC_AC_1_C, eID_a, eID_b);            % reshape partition 1
                                C_2_C = reshape_coherence(coh_PFC_AC_2_C, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_C
                                clear coh_PFC_AC_2_C
                                end 
                    
                              
                                % Maris et al 2007 z-statistic 
                               
                                Coh_test_statistic_a_C = atanh(abs(C_1_C.cohspctrm)) - (1/(C_1_C.dof-2));                        % taken from Maris et al 2007 Eq 1 
                                Coh_test_statistic_b_C = atanh(abs(C_2_C.cohspctrm)) - (1/(C_2_C.dof-2));
                                Coh_test_statistic_c_C = sqrt((1/(C_1_C.dof-2))+(1/(C_2_C.dof-2)));
                                Coh_test_statistic_Z_C = Coh_test_statistic_a_C-Coh_test_statistic_b_C/Coh_test_statistic_c_C;
                                Coh_test_statistic_Zavg_C = mean( Coh_test_statistic_Z_C, 1);                                  % this line collaspes across channel-pair comparisons
                                
                                Monte_Carlo_Array_Maris_C = [Monte_Carlo_Array_Maris_C; Coh_test_statistic_Zavg_C];
                    
                end
            
            end
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1
            
                for i = 1:RandomIteration
                
                % generate a random partition equal to the number of trials in the condition with more trials 
                
                    totalTrials_W = size(All_Trials_w.trial, 2);
                    randomCols_1_W = randperm(totalTrials_W, RandomPartition_w_1_Size); 
                    randomRows_1_W = randomCols_1_W';
                
                    for r = 1:RandomPartition_w_1_Size
                    Random_Partition_w_1.time(:, r)        = All_Trials_w.time(:, randomCols_1_W(r));
                    Random_Partition_w_1.trial(:, r)       = All_Trials_w.trial(:, randomCols_1_W(r));
                    Random_Partition_w_1.sampleinfo(r, :)  = All_Trials_w.sampleinfo(randomRows_1_W(r), :);
                    end
                    
                    Random_Partition_w_1.label      = All_Trials_w(1).label;
                    Random_Partition_w_1.fsample    = All_Trials_w(1).fsample;
                
                    remainingCols_W = setdiff(1:totalTrials_W, randomCols_1_W);
                    randomCols_2_W  = remainingCols_W(randperm(length(remainingCols_W)));
                    randomRows_2_W  = randomCols_2_W'; 
                
                    for rr = 1:RandomPartition_w_2_Size
                    Random_Partition_w_2.time(:, rr)        = All_Trials_w.time(:, randomCols_2_W(rr));
                    Random_Partition_w_2.trial(:, rr)       = All_Trials_w.trial(:, randomCols_2_W(rr));
                    Random_Partition_w_2.sampleinfo(rr, :)  = All_Trials_w.sampleinfo(randomRows_2_W(rr), :);
                    end
                
                    Random_Partition_w_2.label      = All_Trials_w(1).label;
                    Random_Partition_w_2.fsample    = All_Trials_w(1).fsample;
                
                        % non-parametric computation of the cross-spectral density matrix (slow)
                            
                            cfg            = [];
                            cfg.method     = 'mtmfft';
                            cfg.taper      = 'dpss';
                            cfg.output     = 'fourier';
                            cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
                            cfg.pad        = 1;
                            cfg.foilim     = [1 100];
                            
                            freq_1_W       = ft_freqanalysis(cfg, Random_Partition_w_1);
                            freq_2_W       = ft_freqanalysis(cfg, Random_Partition_w_2);
                            
                                % calculate the coherence spectra in partition 1 and partition 2 
                    
                                if strcmp(Channels,'allChan') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb = {'*PFC*','*AC*'};
                            
                                coh_PFC_AC_1_W = ft_connectivityanalysis(cfg, freq_1_W);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_W.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_W);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_W
                        
                                coh_PFC_AC_2_W = ft_connectivityanalysis(cfg, freq_2_W);       % calculate the coherence spectra of partition 2 
                                clear freq_2_W
                               
                                C_1_W = reshape_coherence(coh_PFC_AC_1_W, eID_a, eID_b);            % reshape partition 1
                                C_2_W = reshape_coherence(coh_PFC_AC_2_W, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_W
                                clear coh_PFC_AC_2_W
                                end 
                                
                        
                                if strcmp(Channels,'PFClowmid_ACupmid') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch06' '*PFC_ch11'
                                                   '*AC_ch06' '*PFC_ch12'
                                                   '*AC_ch06' '*PFC_ch13'
                                                   '*AC_ch06' '*PFC_ch14'
                                                   '*AC_ch06' '*PFC_ch15'
                                                   '*AC_ch07' '*PFC_ch11'
                                                   '*AC_ch07' '*PFC_ch12'
                                                   '*AC_ch07' '*PFC_ch13'
                                                   '*AC_ch07' '*PFC_ch14'
                                                   '*AC_ch07' '*PFC_ch15'
                                                   '*AC_ch08' '*PFC_ch11'
                                                   '*AC_ch08' '*PFC_ch12'
                                                   '*AC_ch08' '*PFC_ch13'
                                                   '*AC_ch08' '*PFC_ch14'
                                                   '*AC_ch08' '*PFC_ch15'
                                                   '*AC_ch09' '*PFC_ch11'
                                                   '*AC_ch09' '*PFC_ch12'
                                                   '*AC_ch09' '*PFC_ch13'
                                                   '*AC_ch09' '*PFC_ch14'
                                                   '*AC_ch09' '*PFC_ch15'
                                                   '*AC_ch10' '*PFC_ch11'
                                                   '*AC_ch10' '*PFC_ch12'
                                                   '*AC_ch10' '*PFC_ch13'
                                                   '*AC_ch10' '*PFC_ch14'
                                                   '*AC_ch10' '*PFC_ch15'};   
                                               
                                coh_PFC_AC_1_W = ft_connectivityanalysis(cfg, freq_1_W);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_W.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_W);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_W
                        
                                coh_PFC_AC_2_W = ft_connectivityanalysis(cfg, freq_2_W);       % calculate the coherence spectra of partition 2 
                                clear freq_2_W
                               
                                C_1_W = reshape_coherence(coh_PFC_AC_1_W, eID_a, eID_b);            % reshape partition 1
                                C_2_W = reshape_coherence(coh_PFC_AC_2_W, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_W
                                clear coh_PFC_AC_2_W
                                end 
                            
                                if strcmp(Channels,'PFCdeep_ACupmid') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch06' '*PFC_ch16'
                                                   '*AC_ch06' '*PFC_ch17'
                                                   '*AC_ch06' '*PFC_ch18'
                                                   '*AC_ch06' '*PFC_ch19'
                                                   '*AC_ch06' '*PFC_ch20'
                                                   '*AC_ch07' '*PFC_ch16'
                                                   '*AC_ch07' '*PFC_ch17'
                                                   '*AC_ch07' '*PFC_ch18'
                                                   '*AC_ch07' '*PFC_ch19'
                                                   '*AC_ch07' '*PFC_ch20'
                                                   '*AC_ch08' '*PFC_ch16'
                                                   '*AC_ch08' '*PFC_ch17'
                                                   '*AC_ch08' '*PFC_ch18'
                                                   '*AC_ch08' '*PFC_ch19'
                                                   '*AC_ch08' '*PFC_ch20'
                                                   '*AC_ch09' '*PFC_ch16'
                                                   '*AC_ch09' '*PFC_ch17'
                                                   '*AC_ch09' '*PFC_ch18'
                                                   '*AC_ch09' '*PFC_ch19'
                                                   '*AC_ch09' '*PFC_ch20'
                                                   '*AC_ch10' '*PFC_ch16'
                                                   '*AC_ch10' '*PFC_ch17'
                                                   '*AC_ch10' '*PFC_ch18'
                                                   '*AC_ch10' '*PFC_ch19'
                                                   '*AC_ch10' '*PFC_ch20'};   
                                               
                                coh_PFC_AC_1_W = ft_connectivityanalysis(cfg, freq_1_W);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_W.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_W);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_W
                        
                                coh_PFC_AC_2_W = ft_connectivityanalysis(cfg, freq_2_W);       % calculate the coherence spectra of partition 2 
                                clear freq_2_W
                               
                                C_1_W = reshape_coherence(coh_PFC_AC_1_W, eID_a, eID_b);            % reshape partition 1
                                C_2_W = reshape_coherence(coh_PFC_AC_2_W, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_W
                                clear coh_PFC_AC_2_W
                                end 
                                
                            
                                if strcmp(Channels,'PFCdeep_ACdeep') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch16' '*PFC_ch16'
                                                   '*AC_ch16' '*PFC_ch17'
                                                   '*AC_ch16' '*PFC_ch18'
                                                   '*AC_ch16' '*PFC_ch19'
                                                   '*AC_ch16' '*PFC_ch20'
                                                   '*AC_ch17' '*PFC_ch16'
                                                   '*AC_ch17' '*PFC_ch17'
                                                   '*AC_ch17' '*PFC_ch18'
                                                   '*AC_ch17' '*PFC_ch19'
                                                   '*AC_ch17' '*PFC_ch20'
                                                   '*AC_ch18' '*PFC_ch16'
                                                   '*AC_ch18' '*PFC_ch17'
                                                   '*AC_ch18' '*PFC_ch18'
                                                   '*AC_ch18' '*PFC_ch19'
                                                   '*AC_ch18' '*PFC_ch20'
                                                   '*AC_ch19' '*PFC_ch16'
                                                   '*AC_ch19' '*PFC_ch17'
                                                   '*AC_ch19' '*PFC_ch18'
                                                   '*AC_ch19' '*PFC_ch19'
                                                   '*AC_ch19' '*PFC_ch20'
                                                   '*AC_ch20' '*PFC_ch16'
                                                   '*AC_ch20' '*PFC_ch17'
                                                   '*AC_ch20' '*PFC_ch18'
                                                   '*AC_ch20' '*PFC_ch19'
                                                   '*AC_ch20' '*PFC_ch20'};   
                                               
                                coh_PFC_AC_1_W = ft_connectivityanalysis(cfg, freq_1_W);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_W.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_W);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_W
                        
                                coh_PFC_AC_2_W = ft_connectivityanalysis(cfg, freq_2_W);       % calculate the coherence spectra of partition 2 
                                clear freq_2_W
                               
                                C_1_W = reshape_coherence(coh_PFC_AC_1_W, eID_a, eID_b);            % reshape partition 1
                                C_2_W = reshape_coherence(coh_PFC_AC_2_W, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_W
                                clear coh_PFC_AC_2_W
                                end 
                        
                                if strcmp(Channels,'PFCupper_ACupper') == 1
                                cfg            = [];
                                cfg.method     = 'coh';
                                cfg.channelcmb =  {'*AC_ch03' '*PFC_ch03'
                                                   '*AC_ch03' '*PFC_ch04'
                                                   '*AC_ch03' '*PFC_ch05'
                                                   '*AC_ch03' '*PFC_ch06'
                                                   '*AC_ch03' '*PFC_ch07'
                                                   '*AC_ch04' '*PFC_ch03'
                                                   '*AC_ch04' '*PFC_ch04'
                                                   '*AC_ch04' '*PFC_ch05'
                                                   '*AC_ch04' '*PFC_ch06'
                                                   '*AC_ch04' '*PFC_ch07'
                                                   '*AC_ch05' '*PFC_ch03'
                                                   '*AC_ch05' '*PFC_ch04'
                                                   '*AC_ch05' '*PFC_ch05'
                                                   '*AC_ch05' '*PFC_ch06'
                                                   '*AC_ch05' '*PFC_ch07'
                                                   '*AC_ch06' '*PFC_ch03'
                                                   '*AC_ch06' '*PFC_ch04'
                                                   '*AC_ch06' '*PFC_ch05'
                                                   '*AC_ch06' '*PFC_ch06'
                                                   '*AC_ch06' '*PFC_ch07'
                                                   '*AC_ch07' '*PFC_ch03'
                                                   '*AC_ch07' '*PFC_ch04'
                                                   '*AC_ch07' '*PFC_ch05'
                                                   '*AC_ch07' '*PFC_ch06'
                                                   '*AC_ch07' '*PFC_ch07'};   
                                               
                                coh_PFC_AC_1_W = ft_connectivityanalysis(cfg, freq_1_W);       % calculate the coherence spectra of partition 1 
                                f   = coh_PFC_AC_1_W.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1_W);                % get the electrode ID of the 1st and 2nd columns  
                                clear freq_1_W
                        
                                coh_PFC_AC_2_W = ft_connectivityanalysis(cfg, freq_2_W);       % calculate the coherence spectra of partition 2 
                                clear freq_2_W
                               
                                C_1_W = reshape_coherence(coh_PFC_AC_1_W, eID_a, eID_b);            % reshape partition 1
                                C_2_W = reshape_coherence(coh_PFC_AC_2_W, eID_a, eID_b);            % reshape partition 2
                                clear coh_PFC_AC_1_W
                                clear coh_PFC_AC_2_W
                                end 
                    
                              
                                % Maris et al 2007 z-statistic 
                               
                                Coh_test_statistic_a_W = atanh(abs(C_1_W.cohspctrm)) - (1/(C_1_W.dof-2));                        % taken from Maris et al 2007 Eq 1 
                                Coh_test_statistic_b_W = atanh(abs(C_2_W.cohspctrm)) - (1/(C_2_W.dof-2));
                                Coh_test_statistic_c_W = sqrt((1/(C_1_W.dof-2))+(1/(C_2_W.dof-2)));
                                Coh_test_statistic_Z_W = Coh_test_statistic_a_W-Coh_test_statistic_b_W/Coh_test_statistic_c_W;
                                Coh_test_statistic_Zavg_W = mean( Coh_test_statistic_Z_W, 1);                                  % this line collaspes across channel-pair comparisons
                                
                                Monte_Carlo_Array_Maris_W = [Monte_Carlo_Array_Maris_W; Coh_test_statistic_Zavg_W];
                    
                end
            
            end
            % take the data, threshold, and get the max cluster value (the test statistic)
            
            threshold_value = .5;
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1 
            [data_coh_Maris_Z_C, data_z_statistic_Maris_threshold_C, all_cluster_sums_Maris_C, all_clusters_Maris_C, data_max_cluster_sum_Maris_C] = getDataClusterMax(data_Maris_zstatistic_C,threshold_value);
            end
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1
            [data_coh_Maris_Z_W, data_z_statistic_Maris_threshold_W, all_cluster_sums_Maris_W, all_clusters_Maris_W, data_max_cluster_sum_Maris_W] = getDataClusterMax(data_Maris_zstatistic_W,threshold_value);
            end
            
            % take each permutation, threshold, and get the max cluster values (the test statistic). Now you will have the monte carlo distribution comprised of a test-statistic from each permutation
            
            threshold_value = .5;
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
            [Monte_coh_Z_Maris_C,permutation_zthresholds_Maris_C,thresholded_Permutation_Maris_C, monte_all_cluster_sums_Maris_C, monte_all_clusters_Maris_C, monte_max_cluster_sum_Maris_C] = getPermClusterMax(Monte_Carlo_Array_Maris_C,threshold_value);
            end
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1
            [Monte_coh_Z_Maris_W,permutation_zthresholds_Maris_W,thresholded_Permutation_Maris_W, monte_all_cluster_sums_Maris_W, monte_all_clusters_Maris_W, monte_max_cluster_sum_Maris_W] = getPermClusterMax(Monte_Carlo_Array_Maris_W,threshold_value);
            end
            
            % plot the coherence spectra of the pretone and prior conditions 
            
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
            % calculate the descriptive of the OnlyPrior Condition for correct trials 
             descriptives_onlyprior_coh_avg_c           = squeeze(mean(PFC_AC_Coh_Spectrum_PriorOnly_c.cohspctrm,1));
             descriptives_onlyprior_coh_std_c           = squeeze(std(PFC_AC_Coh_Spectrum_PriorOnly_c.cohspctrm,1));
             descriptives_onlyprior_coh_c_sqrtlength    = sqrt(length(PFC_AC_Coh_Spectrum_PriorOnly_c.cohspctrm));
             descriptives_onlyprior_coh_c_stdErr        = descriptives_onlyprior_coh_std_c/descriptives_onlyprior_coh_c_sqrtlength;
             descriptives_onlyprior_coh_c_CIupper       = descriptives_onlyprior_coh_avg_c  + descriptives_onlyprior_coh_c_stdErr;      
             descriptives_onlyprior_coh_c_CIlower       = descriptives_onlyprior_coh_avg_c  - descriptives_onlyprior_coh_c_stdErr;
            
             % calculate the descriptives of the OnlyPretone Condition for correct trials 
             descriptives_onlypretone_coh_avg_c           = squeeze(mean(PFC_AC_Coh_Spectrum_PretoneOnly_c.cohspctrm,1));
             descriptives_onlypretone_coh_std_c           = squeeze(std(PFC_AC_Coh_Spectrum_PretoneOnly_c.cohspctrm,1));
             descriptives_onlypretone_coh_c_sqrtlength    = sqrt(length(PFC_AC_Coh_Spectrum_PretoneOnly_c.cohspctrm));
             descriptives_onlypretone_coh_c_stdErr        = descriptives_onlypretone_coh_std_c/descriptives_onlypretone_coh_c_sqrtlength;
             descriptives_onlypretone_coh_c_CIupper       = descriptives_onlypretone_coh_avg_c  + descriptives_onlypretone_coh_c_stdErr;      
             descriptives_onlypretone_coh_c_CIlower       = descriptives_onlypretone_coh_avg_c  - descriptives_onlypretone_coh_c_stdErr;
            end
            
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1
             % calculate the descriptive of the OnlyPrior Condition for wrong trials 
             descriptives_onlyprior_coh_avg_w          = squeeze(mean(PFC_AC_Coh_Spectrum_PriorOnly_w.cohspctrm,1));
             descriptives_onlyprior_coh_std_w           = squeeze(std(PFC_AC_Coh_Spectrum_PriorOnly_w.cohspctrm,1));
             descriptives_onlyprior_coh_w_sqrtlength    = sqrt(length(PFC_AC_Coh_Spectrum_PriorOnly_w.cohspctrm));
             descriptives_onlyprior_coh_w_stdErr        = descriptives_onlyprior_coh_std_w/descriptives_onlyprior_coh_w_sqrtlength;
             descriptives_onlyprior_coh_w_CIupper       = descriptives_onlyprior_coh_avg_w  + descriptives_onlyprior_coh_w_stdErr;      
             descriptives_onlyprior_coh_w_CIlower       = descriptives_onlyprior_coh_avg_w  - descriptives_onlyprior_coh_w_stdErr;
            
             % calculate the descriptives of the OnlyPretone Condition for wrong trials 
             descriptives_onlypretone_coh_avg_w           = squeeze(mean(PFC_AC_Coh_Spectrum_PretoneOnly_w.cohspctrm,1));
             descriptives_onlypretone_coh_std_w           = squeeze(std(PFC_AC_Coh_Spectrum_PretoneOnly_w.cohspctrm,1));
             descriptives_onlypretone_coh_w_sqrtlength    = sqrt(length(PFC_AC_Coh_Spectrum_PretoneOnly_w.cohspctrm));
             descriptives_onlypretone_coh_w_stdErr        = descriptives_onlypretone_coh_std_w/descriptives_onlypretone_coh_w_sqrtlength;
             descriptives_onlypretone_coh_w_CIupper       = descriptives_onlypretone_coh_avg_w  + descriptives_onlypretone_coh_w_stdErr;      
             descriptives_onlypretone_coh_w_CIlower       = descriptives_onlypretone_coh_avg_w  - descriptives_onlypretone_coh_w_stdErr;
            end
            
            %
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
                if plot == 1
            
            figure
            
            subplot(2,3,1)
            hold on
            ciplot(descriptives_onlyprior_coh_c_CIupper(1, 5:8) ,descriptives_onlyprior_coh_c_CIlower(1, 5:8), 'blue');
            plot(descriptives_onlyprior_coh_avg_c(1, 5:8), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_c_CIupper(1, 5:8) ,descriptives_onlypretone_coh_c_CIlower(1, 5:8), 'green');
            plot(descriptives_onlypretone_coh_avg_c(1, 5:8), 'color', [0 0 0]);
            title ('theta')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 4])
            hold off
            
            subplot(2,3,2)
            hold on
            ciplot(descriptives_onlyprior_coh_c_CIupper(1, 8:14) ,descriptives_onlyprior_coh_c_CIlower(1, 8:14), 'blue');
            plot(descriptives_onlyprior_coh_avg_c(1, 8:14), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_c_CIupper(1, 8:14) ,descriptives_onlypretone_coh_c_CIlower(1, 8:14), 'green');
            plot(descriptives_onlypretone_coh_avg_c(1, 8:14), 'color', [0 0 0]);
            title ('alpha')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 7])
            hold off
            
            subplot(2,3,3)
            hold on
            ciplot(descriptives_onlyprior_coh_c_CIupper(1, 15:30) ,descriptives_onlyprior_coh_c_CIlower(1, 15:30), 'blue');
            plot(descriptives_onlyprior_coh_avg_c(1, 15:30), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_c_CIupper(1, 15:30) ,descriptives_onlypretone_coh_c_CIlower(1, 15:30), 'green');
            plot(descriptives_onlypretone_coh_avg_c(1, 15:30), 'color', [0 0 0]);
            title ('beta')
            legend('LED Only', '' ,'Pretone Only','')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 16])
            hold off
            
            subplot(2,3,4)
            hold on
            ciplot(descriptives_onlyprior_coh_c_CIupper(1, 31:54) ,descriptives_onlyprior_coh_c_CIlower(1, 31:54), 'blue');
            plot(descriptives_onlyprior_coh_avg_c(1, 31:54), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_c_CIupper(1, 31:54) ,descriptives_onlypretone_coh_c_CIlower(1, 31:54), 'green');
            plot(descriptives_onlypretone_coh_avg_c(1, 31:54), 'color', [0 0 0]);
            title ('gamma')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 24])
            hold off
            
            subplot(2,3,5)
            hold on
            ciplot(descriptives_onlyprior_coh_c_CIupper(1, 66:100) ,descriptives_onlyprior_coh_c_CIlower(1, 66:100), 'blue');
            plot(descriptives_onlyprior_coh_avg_c(1, 66:100), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_c_CIupper(1, 66:100) ,descriptives_onlypretone_coh_c_CIlower(1, 66:100), 'green');
            plot(descriptives_onlypretone_coh_avg_c(1, 66:100), 'color', [0 0 0]);
            title ('high gamma')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 35])
            hold off
            
            subplot(2,3,6)
            hold on
            ciplot(descriptives_onlyprior_coh_c_CIupper,descriptives_onlyprior_coh_c_CIlower, 'blue');
            plot(descriptives_onlyprior_coh_avg_c, 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_c_CIupper,descriptives_onlypretone_coh_c_CIlower, 'green');
            plot(descriptives_onlypretone_coh_avg_c, 'color', [0 0 0]);
            title ('full spectra')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 101])
            hold off

                end
            end
            
            %
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1
                if plot == 1
            figure
            
            subplot(2,3,1)
            hold on
            ciplot(descriptives_onlyprior_coh_w_CIupper(1, 5:8) ,descriptives_onlyprior_coh_w_CIlower(1, 5:8), 'blue');
            plot(descriptives_onlyprior_coh_avg_w(1, 5:8), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_w_CIupper(1, 5:8) ,descriptives_onlypretone_coh_w_CIlower(1, 5:8), 'green');
            plot(descriptives_onlypretone_coh_avg_w(1, 5:8), 'color', [0 0 0]);
            title ('theta')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 4])
            hold off
            
            subplot(2,3,2)
            hold on
            ciplot(descriptives_onlyprior_coh_w_CIupper(1, 8:14) ,descriptives_onlyprior_coh_w_CIlower(1, 8:14), 'blue');
            plot(descriptives_onlyprior_coh_avg_w(1, 8:14), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_w_CIupper(1, 8:14) ,descriptives_onlypretone_coh_w_CIlower(1, 8:14), 'green');
            plot(descriptives_onlypretone_coh_avg_w(1, 8:14), 'color', [0 0 0]);
            title ('alpha')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 7])
            hold off
            
            subplot(2,3,3)
            hold on
            ciplot(descriptives_onlyprior_coh_w_CIupper(1, 15:30) ,descriptives_onlyprior_coh_w_CIlower(1, 15:30), 'blue');
            plot(descriptives_onlyprior_coh_avg_w(1, 15:30), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_w_CIupper(1, 15:30) ,descriptives_onlypretone_coh_w_CIlower(1, 15:30), 'green');
            plot(descriptives_onlypretone_coh_avg_w(1, 15:30), 'color', [0 0 0]);
            title ('beta')
            legend('LED Only', '' ,'Pretone Only','')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 16])
            hold off
            
            subplot(2,3,4)
            hold on
            ciplot(descriptives_onlyprior_coh_w_CIupper(1, 31:54) ,descriptives_onlyprior_coh_w_CIlower(1, 31:54), 'blue');
            plot(descriptives_onlyprior_coh_avg_w(1, 31:54), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_w_CIupper(1, 31:54) ,descriptives_onlypretone_coh_w_CIlower(1, 31:54), 'green');
            plot(descriptives_onlypretone_coh_avg_w(1, 31:54), 'color', [0 0 0]);
            title ('gamma')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 24])
            hold off
            
            subplot(2,3,5)
            hold on
            ciplot(descriptives_onlyprior_coh_w_CIupper(1, 66:100) ,descriptives_onlyprior_coh_w_CIlower(1, 66:100), 'blue');
            plot(descriptives_onlyprior_coh_avg_w(1, 66:100), 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_w_CIupper(1, 66:100) ,descriptives_onlypretone_coh_w_CIlower(1, 66:100), 'green');
            plot(descriptives_onlypretone_coh_avg_w(1, 66:100), 'color', [0 0 0]);
            title ('high gamma')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 35])
            hold off
            
            subplot(2,3,6)
            hold on
            ciplot(descriptives_onlyprior_coh_w_CIupper,descriptives_onlyprior_coh_w_CIlower, 'blue');
            plot(descriptives_onlyprior_coh_avg_w, 'color', [0 0 0]);
            ciplot(descriptives_onlypretone_coh_w_CIupper,descriptives_onlypretone_coh_w_CIlower, 'green');
            plot(descriptives_onlypretone_coh_avg_w, 'color', [0 0 0]);
            title ('full spectra')
            ylabel ('coherence')
            xlabel ('frequency (Hz)')
            xlim([1 101])
            hold off
                end 
            end
            
            % plot the monte carlo distribution of the Maris test-statistics (cluster max) with the data for correct trials 
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
                if plot == 1 
            figure 
            subplot(2,3,1)
            histogram(monte_max_cluster_sum_Maris_C.theta, 100);
            hold on 
            line([data_max_cluster_sum_Maris_C.theta ,data_max_cluster_sum_Maris_C.theta ], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('theta')
            
            subplot(2,3,2)
            histogram(monte_max_cluster_sum_Maris_C.alpha, 100);
            hold on 
            line([data_max_cluster_sum_Maris_C.alpha,data_max_cluster_sum_Maris_C.alpha ], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('alpha')
            
            subplot(2,3,3)
            histogram(monte_max_cluster_sum_Maris_C.beta, 100);
            hold on 
            line([data_max_cluster_sum_Maris_C.beta ,data_max_cluster_sum_Maris_C.beta], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('beta')
            legend('permutation', 'data')
            
            subplot(2,3,4)
            histogram(monte_max_cluster_sum_Maris_C.gamma, 100);
            hold on 
            line([data_max_cluster_sum_Maris_C.gamma ,data_max_cluster_sum_Maris_C.gamma], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('gamma')
            
            subplot(2,3,5)
            histogram(monte_max_cluster_sum_Maris_C.highGamma, 100);
            hold on 
            line([data_max_cluster_sum_Maris_C.highGamma ,data_max_cluster_sum_Maris_C.highGamma], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('high gamma')
            
            subplot(2,3,6)
            histogram(monte_max_cluster_sum_Maris_C.all, 100);
            hold on 
            line([data_max_cluster_sum_Maris_C.all ,data_max_cluster_sum_Maris_C.all], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('all')
                end
            end
            
            % plot the monte carlo distribution of the Maris test-statistics (cluster max) with the data for wrong trials 
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')== 1
                if plot == 1 
            
            figure 
            subplot(2,3,1)
            histogram(monte_max_cluster_sum_Maris_W.theta, 100);
            hold on 
            line([data_max_cluster_sum_Maris_W.theta ,data_max_cluster_sum_Maris_W.theta ], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('theta')
            
            subplot(2,3,2)
            histogram(monte_max_cluster_sum_Maris_W.alpha, 100);
            hold on 
            line([data_max_cluster_sum_Maris_W.alpha,data_max_cluster_sum_Maris_W.alpha ], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('alpha')
            
            subplot(2,3,3)
            histogram(monte_max_cluster_sum_Maris_W.beta, 100);
            hold on 
            line([data_max_cluster_sum_Maris_W.beta ,data_max_cluster_sum_Maris_W.beta], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('beta')
            legend('permutation', 'data')
            
            subplot(2,3,4)
            histogram(monte_max_cluster_sum_Maris_W.gamma, 100);
            hold on 
            line([data_max_cluster_sum_Maris_W.gamma ,data_max_cluster_sum_Maris_W.gamma], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('gamma')
            
            subplot(2,3,5)
            histogram(monte_max_cluster_sum_Maris_W.highGamma, 100);
            hold on 
            line([data_max_cluster_sum_Maris_W.highGamma ,data_max_cluster_sum_Maris_W.highGamma], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('high gamma')
            
            subplot(2,3,6)
            histogram(monte_max_cluster_sum_Maris_W.all, 100);
            hold on 
            line([data_max_cluster_sum_Maris_W.all ,data_max_cluster_sum_Maris_W.all], ylim, 'Color', 'r', 'LineWidth', 2)
            hold off
            title('all')
                end
            end
            
            % calculate the p-value for correct trials 
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
               
                pvalue_C.theta = sum(monte_max_cluster_sum_Maris_C.theta <= data_max_cluster_sum_Maris_C.theta) / numel(monte_max_cluster_sum_Maris_C.theta);
                pvalue_C.alpha = sum(monte_max_cluster_sum_Maris_C.alpha <= data_max_cluster_sum_Maris_C.alpha) / numel(monte_max_cluster_sum_Maris_C.alpha);
                pvalue_C.beta = sum(monte_max_cluster_sum_Maris_C.beta <= data_max_cluster_sum_Maris_C.beta) / numel(monte_max_cluster_sum_Maris_C.beta);
                pvalue_C.gamma = sum(monte_max_cluster_sum_Maris_C.gamma <= data_max_cluster_sum_Maris_C.gamma) / numel(monte_max_cluster_sum_Maris_C.gamma);
                pvalue_C.highGamma = sum(monte_max_cluster_sum_Maris_C.highGamma <= data_max_cluster_sum_Maris_C.highGamma) / numel(monte_max_cluster_sum_Maris_C.highGamma);
                pvalue_C.all = sum(monte_max_cluster_sum_Maris_C.all <= data_max_cluster_sum_Maris_C.all) / numel(monte_max_cluster_sum_Maris_C.all);
                
                if pvalue_C.theta > .5 
                pvalue_C.theta = 1 - pvalue_C.theta;
                end 
                
                if pvalue_C.alpha > .5 
                pvalue_C.alpha = 1 - pvalue_C.alpha;
                end 
                
                if pvalue_C.beta > .5 
                pvalue_C.beta = 1 - pvalue_C.beta;
                end 
                
                if pvalue_C.gamma > .5 
                pvalue_C.gamma = 1 - pvalue_C.gamma;
                end 
                
                if pvalue_C.highGamma > .5 
                pvalue_C.highGamma = 1 - pvalue_C.highGamma;
                end
                
                if pvalue_C.all > .5 
                pvalue_C.all = 1 - pvalue_C.all;
                end
            
        
            end
            % calculate the p-value for wrong trials 
            
            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1
            
                pvalue_W.theta = sum(monte_max_cluster_sum_Maris_W.theta <= data_max_cluster_sum_Maris_W.theta) / numel(monte_max_cluster_sum_Maris_W.theta);
                pvalue_W.alpha = sum(monte_max_cluster_sum_Maris_W.alpha <= data_max_cluster_sum_Maris_W.alpha) / numel(monte_max_cluster_sum_Maris_W.alpha);
                pvalue_W.beta = sum(monte_max_cluster_sum_Maris_W.beta <= data_max_cluster_sum_Maris_W.beta) / numel(monte_max_cluster_sum_Maris_W.beta);
                pvalue_W.gamma = sum(monte_max_cluster_sum_Maris_W.gamma <= data_max_cluster_sum_Maris_W.gamma) / numel(monte_max_cluster_sum_Maris_W.gamma);
                pvalue_W.highGamma = sum(monte_max_cluster_sum_Maris_W.highGamma <= data_max_cluster_sum_Maris_W.highGamma) / numel(monte_max_cluster_sum_Maris_W.highGamma);
                pvalue_W.all = sum(monte_max_cluster_sum_Maris_W.all <= data_max_cluster_sum_Maris_W.all) / numel(monte_max_cluster_sum_Maris_W.all);
                
                if pvalue_W.theta > .5 
                pvalue_W.theta = 1 - pvalue_W.theta;
                end 
                
                if pvalue_W.alpha > .5 
                pvalue_W.alpha = 1 - pvalue_W.alpha;
                end 
                
                if pvalue_W.beta > .5 
                pvalue_W.beta = 1 - pvalue_W.beta;
                end 
                
                if pvalue_W.gamma > .5 
                pvalue_W.gamma = 1 - pvalue_W.gamma;
                end 
                
                if pvalue_W.highGamma > .5 
                pvalue_W.highGamma = 1 - pvalue_W.highGamma;
                end
                
                if pvalue_W.all > .5 
                pvalue_W.all = 1 - pvalue_W.all;
                end
            end
            
            %

            outputDirectory = fullfile(savedir, sessions(s).name);
                    if ~exist(outputDirectory, 'dir')
                     mkdir(outputDirectory);
                    end

            if isSave == 1 
            
                save(fullfile(outputDirectory, save_file_name), 'data_coh_Maris_Z_C', 'data_z_statistic_Maris_threshold_C', 'all_cluster_sums_Maris_C', 'all_clusters_Maris_C', 'data_max_cluster_sum_Maris_C',...
                                                    'data_coh_Maris_Z_W', 'data_z_statistic_Maris_threshold_W', 'all_cluster_sums_Maris_W', 'all_clusters_Maris_W', 'data_max_cluster_sum_Maris_W',...
                                                    'data_Maris_zstatistic_C', 'data_Maris_zstatistic_W', 'threshold_value', 'Monte_coh_Z_Maris_C','permutation_zthresholds_Maris_C','thresholded_Permutation_Maris_C',...
                                                    'monte_all_cluster_sums_Maris_C', 'monte_all_clusters_Maris_C', 'monte_max_cluster_sum_Maris_C','Monte_coh_Z_Maris_W','permutation_zthresholds_Maris_W',...
                                                    'thresholded_Permutation_Maris_W', 'monte_all_cluster_sums_Maris_W', 'monte_all_clusters_Maris_W', 'monte_max_cluster_sum_Maris_W', 'pvalue_C','pvalue_W');
            
            
  end
        end
          
