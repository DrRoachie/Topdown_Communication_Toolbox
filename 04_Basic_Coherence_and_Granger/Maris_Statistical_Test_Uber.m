%% Script Summary

% Takes the pre-processed Taku data and calculates coherence. This
% particular file takes the coherence and averages across 4 groups of
% layers. 


%% Set directory and get a list of all files in the folder with the desired file name pattern.

datadir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2\MrCassius\testTone\190418';
chandir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_05_25_Analysis\01_Layer_Groupings\MrCassius\testTone\190418';
savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_05_20_Analysis\02_Statistical_Test_Results\MrCassius\testTone\190418';
sessions = dir(fullfile(datadir,'*.mat'));
addpath(genpath(datadir));
addpath(genpath(chandir));
isSave = 0;  

%% Establish the subset of data that you want to process. 
Animal         = 'MrCassius';                  % Options: 'MrCassius', 'MrM'; 
RecDate        = '190418';        
Epoch          = 'testToneOnset';                % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
Behavior       = 'Correct';                    % Options:  'Correct', Wrong, due to memory constraints                          
RandomIteration  = 300;                        % options: the number of times you want to allocate the random partition, and generate the test statistic

%%

freq_band_list          = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};
LayerGrouping_List     = {'PFClowmid_ACupmid';'PFCupmid_ACupmid'; 'PFCdeep_ACupmid'; 'PFCdeep_ACdeep'; 'PFCupper_ACupper'};
%LayerGrouping_List      = {'PFClowmid_ACupmid'};
threshold_value         = 0.5;       % the threshold used for the monte carlo test between PriorOnly and PretoneOnly Conditions

if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1

    for fb = 1:length(freq_band_list)

       Frequency_Band = freq_band_list{fb};
  
       % Get the channel groupings with only spectrally significant channels for both PretoneOnly and PriorOnly trial
            PriorOnly_NO_SigGroups_fn   = sprintf('PriorOnly_NO_SigGroups_%s_%s_%s_%s.txt', RecDate, Epoch, Behavior, Frequency_Band);
            PretoneOnly_NO_SigGroups_fn = sprintf('PretoneOnly_NO_SigGroups_%s_%s_%s_%s.txt', RecDate, Epoch, Behavior, Frequency_Band);

            nullcheck_1 = dir(fullfile(chandir, PriorOnly_NO_SigGroups_fn));
            nullcheck_2 = dir(fullfile(chandir, PretoneOnly_NO_SigGroups_fn));

            if ~isempty(nullcheck_1)||~isempty(nullcheck_2)

                        for lc = 1:length(LayerGrouping_List)
            
                            LayerGroup = LayerGrouping_List{lc};
                         
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_zthresh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_montecarlo.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_data.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose(fid1);
                                        fclose(fid2);
                                        fclose(fid3);
                                        fclose(fid4);

                                        clear nullcheck_1
                                        clear nullcheck_2
                                  
                        end

             else

                PriorOnly_SigChans_c_fn   = sprintf('PriorOnly_GroupedSigChans_%s_%s_%s_%s.mat', RecDate, Epoch, Behavior, Frequency_Band);
                files = dir(fullfile(chandir, PriorOnly_SigChans_c_fn));
                
                if ~isempty(files)
                   PriorOnly_SigChans_c = load(fullfile(chandir, PriorOnly_SigChans_c_fn), 'PriorOnly_Chans_Filtered');         % for each frequency band, pull the spectrally relevant channels in the PriorOnly Condition 
                end
            
                PretoneOnly_SigChans_c_fn = sprintf('PretoneOnly_GroupedSigChans_%s_%s_%s_%s.mat', RecDate, Epoch, Behavior, Frequency_Band);
                files = dir(fullfile(chandir, PretoneOnly_SigChans_c_fn));
                
                if ~isempty(files)
                    PretoneOnly_SigChans_c = load(fullfile(chandir, PretoneOnly_SigChans_c_fn), 'PretoneOnly_Chans_Filtered');   % for each frequency band, pull the spectrally relevant channels in the PretoneOnly Condition 
                end

                        % Carry out the Coherence Test for every layer grouping, graph, and save data
                      
                        for lc = 1:length(LayerGrouping_List)
            
                            LayerGroup = LayerGrouping_List{lc};
                            Condition = 'OnlyPrior';                          % options: 'OnlyPrior', 'OnlyPretone', & 'Both';
                            Channels  = PriorOnly_SigChans_c.PriorOnly_Chans_Filtered.(LayerGroup);
            
                                    if isempty(Channels)
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_zthresh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_montecarlo.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_data.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose(fid1);
                                        fclose(fid2);
                                        fclose(fid3);
                                        fclose(fid4);
                                    else
            
                            [PFC_AC_Coh_Spectrum_PriorOnly_c,PFC_AC_Coh_fd_PriorOnly_c,TrialInfo_PriorOnly_c] = getDirectoryCoherence(sessions, Condition, Channels, Behavior);              % Get the PriorOnly Coherence 
                       
                            Condition = 'OnlyPretone';                          % options: 'OnlyPrior', 'OnlyPretone', & 'Both';
                            Channels  = PretoneOnly_SigChans_c.PretoneOnly_Chans_Filtered.(LayerGroup);
            
                                    if isempty(Channels)
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_zthresh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_montecarlo.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_data.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose(fid1);
                                        fclose(fid2);
                                        fclose(fid3);
                                        fclose(fid4);
                                    else
            
                            [PFC_AC_Coh_Spectrum_PretoneOnly_c,PFC_AC_Coh_fd_PretoneOnly_c,TrialInfo_PretoneOnly_c] = getDirectoryCoherence(sessions, Condition, Channels, Behavior);         % Get PretoneOnly Coherence 
                                                                                   
                            [data_coh_Maris_Z_c, data_z_statistic_Maris_threshold_c, all_cluster_sums_Maris_c, all_clusters_Maris_c, data_max_cluster_sum_Maris_c,... 
                            Monte_coh_Z_Maris_c, permutation_zthresholds_Maris_c, thresholded_Permutation_Maris_c,...
                            monte_all_cluster_sums_Maris_c, monte_all_clusters_Maris_c, monte_max_cluster_sum_Maris_c, pvalue_c] = RunMaris(PFC_AC_Coh_Spectrum_PriorOnly_c, PFC_AC_Coh_Spectrum_PretoneOnly_c,...  % Run the Statistical Test between the OnlyPrior and OnlyPretone Conditions
                                                                                                                                   datadir, sessions, Animal, Behavior, Epoch, Channels, threshold_value, RandomIteration);
                        
                            [descriptives_onlyprior_coh_c, descriptives_onlypretone_coh_c] = GetCohDesciptives(PFC_AC_Coh_Spectrum_PriorOnly_c, PFC_AC_Coh_Spectrum_PretoneOnly_c); % get the coherence descriptions ready
                            Plot_CoherenceSpec(Frequency_Band, descriptives_onlyprior_coh_c, descriptives_onlypretone_coh_c)                                                        % plot the coherence spectra of PriorOnly and PretoneOnly
                            Plot_ThreshedZ(Frequency_Band, data_coh_Maris_Z_c, data_z_statistic_Maris_threshold_c)                                                                  % plot the z-spectra of the two conditions with threshold
                            Plot_MonteCarloDist(Frequency_Band,monte_max_cluster_sum_Maris_c,data_max_cluster_sum_Maris_c)                                                          % plot monte-carlo distibution with data line 
            
                            Fig1_name      = sprintf('%s_%s_%s_%s_coh.fig', RecDate, Frequency_Band, LayerGroup, Behavior);
                            Fig2_name      = sprintf('%s_%s_%s_%s_zthresh.fig', RecDate, Frequency_Band, LayerGroup, Behavior);
                            Fig3_name      = sprintf('%s_%s_%s_%s_montecarlo.fig', RecDate, Frequency_Band, LayerGroup, Behavior);
                            save_file_name = sprintf('%s_%s_%s_%s_data.mat', RecDate, Frequency_Band, LayerGroup, Behavior);
                            
                            saveas(1, fullfile(savedir, Fig1_name));
                            saveas(2, fullfile(savedir, Fig2_name));
                            saveas(3, fullfile(savedir, Fig3_name));
            
                            save(fullfile(savedir, save_file_name), 'data_coh_Maris_Z_c', 'data_z_statistic_Maris_threshold_c', 'all_cluster_sums_Maris_c', 'all_clusters_Maris_c', 'data_max_cluster_sum_Maris_c',...
                                                                    'threshold_value', 'Monte_coh_Z_Maris_c','permutation_zthresholds_Maris_c','thresholded_Permutation_Maris_c',...
                                                                    'monte_all_cluster_sums_Maris_c', 'monte_all_clusters_Maris_c', 'monte_max_cluster_sum_Maris_c', 'pvalue_c');
            
                            close all

                                  end
                                end
                              end
            
                        end
            
                 end
            
             end



      
            %% plot the coherence spectra, the z spectra with threshold, and the monte carlo distribution with data line 
            % % 
            % if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1
            % [descriptives_onlyprior_coh_c, descriptives_onlypretone_coh_c] = GetCohDesciptives(PFC_AC_Coh_Spectrum_PriorOnly_c, PFC_AC_Coh_Spectrum_PretoneOnly_c);
            % Plot_CoherenceSpec(Frequency_Band, descriptives_onlyprior_coh_c, descriptives_onlypretone_coh_c)
            % Plot_ThreshedZ(Frequency_Band, data_coh_Maris_Z_c, data_z_statistic_Maris_threshold_c)
            % Plot_MonteCarloDist(Frequency_Band,monte_max_cluster_sum_Maris_c,data_max_cluster_sum_Maris_c)
            % end



