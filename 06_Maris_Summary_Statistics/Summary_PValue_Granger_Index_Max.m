
% This code takes the granger spectra from channel pairs that are
% significantly modulated by context, and tests the correlation between
% directionalities (PFC-to-AC) or (AC-to-PFC) within context condition (OnlyPrior or OnlyPretone) or
% between context conditions within a single direction. 

%Note: there is also code that calculates a simple granger index and uses
%that as a metric to compare across conditions. 

%This code also carries out a paired t-test between the two conditions

%% Define data

Frequency_Band  = 'theta';               % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'Granger';             % 'Granger' or 'Coherence'
Animal          = 'MrM';                 % 'MrCassius' and/or 'MrM'
Direction       = 'P2A';                 % 'P2A' or 'A2P'      
freq_band_list  = {'theta'};
animals         = {'MrM'};
Epoch           = 'testToneOnset';

%% load histogram data arrays and max channel pairs
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

maxdir          = 'D:\2024_09_27_Analysis\XCorr_Histogram_Data\COHERENCE_testTone_xcorr_histogram_data';
path_to_data    = 'D:\2024_09_27_Analysis\XCorr_Histogram_Data\COHERENCE_testTone_xcorr_histogram_data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;

%%

rootdir = 'D:\2024_09_27_Analysis';
sessions = dir(fullfile(rootdir, '19*'));

ChanPair_Array_PFC_AC = zeros(20,20);   % used to store passing p-value count for granger (bidirectional)
ChanPair_Array_AC_PFC = zeros(20,20);

%% Initialize the arrays for comparing the Granger Index values across
% condition

PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs      = [];
PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs    = [];
PriorOnly_Granger_Index_P2A_Sig_ChannelPairs         = [];
PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs       = [];

PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs      = [];
PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs    = [];
PriorOnly_Granger_Index_A2P_Sig_ChannelPairs         = [];
PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs       = [];

% Initialize the arrays for comparing the non-normalized granger values
% across directionality  

PriorOnly_Granger_Value_P2A_Sig_ChannelPairs         = [];
PriorOnly_Granger_Value_A2P_Sig_ChannelPairs         = [];

PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs       = [];
PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs       = [];

PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs      = [];
PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs      = [];

PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs    = [];
PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs    = [];


%% pull the channel labels of the strongest connection according to the spatial distribution of channel pairs whose coherence is signigicantly modulated by context

for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    fName   = sprintf('pvalue_max_channels_%s_%s_Coh.mat', Animal, Frequency_Band);
    load(fullfile(maxdir, fName));
    max_results_PFC_AC = max_results;
    max_results_PFC_AC = max_results_PFC_AC(max_results_PFC_AC(:, 1)==1, :); % only get the strongest connection
    clear max_results;

    fName   = sprintf('pvalue_max_channels_%s_%s_Coh.mat', Animal, Frequency_Band);
    load(fullfile(maxdir, fName));
    max_results_AC_PFC = max_results;
    max_results_AC_PFC = max_results_AC_PFC(max_results_AC_PFC(:, 1)==1, :);
    clear max_results;

    max_chans_PFC_AC      = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3)];  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
    max_chans_AC_PFC      = [max_results_AC_PFC(:,3) max_results_AC_PFC(:,2)];  % IMPORTANT NOTE: whereas in the xcorr histogram code we need to flip the columns to adjust direction and extract the unique rows, granger needs to preserve the directionality

end


%% Loop through all the sessions and pull all granger files in the PFC to AC direction that are in the max channel list regardless of whether it is significant or not.

PFC_AC_matching_files = {};

for i = 1:length(sessions) 

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};

        if exist(fullfile(rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
            
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), ['*' Frequency_Band '*Granger_data.mat']));
             
                 if strcmp(Direction, 'P2A')
                    
                     matching_files = {};

                    for k = 1:length(files)
                        
                            % Extract the file name from the 'name' field for the current row
                            file_name = files(k).name;
    
                            % Use regular expressions to extract PFC and AC channel numbers
                            pfc_match = regexp(file_name, 'PFC_ch(\d+)', 'tokens');
                            ac_match  = regexp(file_name, 'AC_ch(\d+)', 'tokens');
                    
                            if ~isempty(pfc_match) && ~isempty(ac_match)
                                
                                % Convert the extracted channel numbers to numeric format
                                pfc_ch = str2double(pfc_match{1}{1});
                                ac_ch  = str2double(ac_match{1}{1});
                                
                                % Check if this pair exists in max_chans
                                is_in_max_chans_PFC_AC = any(ismember(max_chans_PFC_AC, [pfc_ch, ac_ch], 'rows'));
                                
                                if is_in_max_chans_PFC_AC
                                    % Append the matching file name to the cell array
                                    matching_files{end+1} = file_name;
                                    
                                end
                            
                            end
                       end

                 end

                if ~isempty(matching_files)
                PFC_AC_matching_files = [PFC_AC_matching_files; matching_files']; %#ok<AGROW> % Transpose to ensure it's a column array
                end

        end
   
    end
    
end


%% Go through all the channels in the strongest connection and pull any channel pair that were contextual modulated

        %Loop through the matching files
        
if strcmp(Direction, 'P2A')

        for i = 1:length(PFC_AC_matching_files)

            % Get the current file name
            matching_file = PFC_AC_matching_files{i};

            % Extract the session ID from the file name (the second component, e.g., '190722')
            split_file_name = split(matching_file, '_'); 

            % Components of the file name (e.g., MrM_190722_testToneOnset_PFC_ch19_AC_ch18_theta_Correct_Granger_data.mat)
            animal_id = split_file_name{1};    % e.g., 'MrM'
            session_id = split_file_name{2};   % e.g., '190722'
            epoch_full = split_file_name{3};        % e.g., 'testToneOnset'
            epoch = extractBefore(epoch_full, 'Onset');   % Extract just 'testTone'

            % Construct the full path to the data file
            session_dir = fullfile(rootdir, session_id, animal_id, epoch);

            % Construct the full path to the data file within the session directory
            file_path = fullfile(session_dir, matching_file);

             % Optional: Print the file path to check if it's correctly built
            disp(file_path); % Ensure this shows the full, correct path

            % Check if the file actually exists
            if isfile(file_path)  % This checks if the file actually exists
                fprintf('Found file: %s\n', file_path);
            end

            Granger_Data = load(file_path); % load the current file from the list of files found in the maximum connection

            if isfield( Granger_Data, 'pvalue_PFC_AC') && isfield( Granger_Data, 'pvalue_AC_PFC')
              
                    fprintf('Variables pvalue_PFC_AC and pvalue_AC_PFC loaded successfully. Proceeding with analysis...\n');
                    
                    pvalue_PFC_AC = Granger_Data.pvalue_PFC_AC.(Frequency_Band);
               
                        % Sig Channel Pairs Section

                            if pvalue_PFC_AC <= 0.05


                                % calculate PriorOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              

                                element_1 = [];
                                    for z = 1:length(Granger_Data.Granger_PriorOnly.grangerspctrm)
                                    element_1 = [element_1; Granger_Data.Granger_PriorOnly.grangerspctrm(:,:, z)];
                                    end
                                element_1 = element_1(element_1(:,1) ~= 0, 1);
                                element_1 = mean(element_1); % collapse across frequencies within band

                                PriorOnly_Granger_Value_P2A_Sig_ChannelPairs = [PriorOnly_Granger_Value_P2A_Sig_ChannelPairs, element_1];

                                % calculate PriorOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    

                                element_2 = [];
                                    for z = 1:length(Granger_Data.Granger_PriorOnly.grangerspctrm)
                                    element_2 = [element_2; Granger_Data.Granger_PriorOnly.grangerspctrm(:,:, z)];
                                    end
                                element_2 = element_2(element_2(:,2) ~= 0, 2);
                                element_2 = mean(element_2);

                                PriorOnly_Granger_Value_A2P_Sig_ChannelPairs = [PriorOnly_Granger_Value_A2P_Sig_ChannelPairs, element_2];

                                % calculate PriorOnly granger index 

                                granger_index = (element_1 - element_2)/(element_1 + element_2);
                                PriorOnly_Granger_Index_P2A_Sig_ChannelPairs = [PriorOnly_Granger_Index_P2A_Sig_ChannelPairs; granger_index];

                                % calculate PretoneOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              

                                element_1 = [];
                                    for z = 1:length(Granger_Data.Granger_PretoneOnly.grangerspctrm)
                                    element_1 = [element_1; Granger_Data.Granger_PretoneOnly.grangerspctrm(:,:, z)];
                                    end
                                element_1 = element_1(element_1(:,1) ~= 0, 1);
                                element_1 = mean(element_1);

                                PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs = [PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs, element_1];

                                % calculate PretoneOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    

                                element_2 = [];
                                    for z = 1:length(Granger_Data.Granger_PretoneOnly.grangerspctrm)
                                    element_2 = [element_2; Granger_Data.Granger_PretoneOnly.grangerspctrm(:,:, z)];
                                    end
                                element_2 = element_2(element_2(:,2) ~= 0, 2);
                                element_2 = mean(element_2);

                                PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs = [PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs, element_2];

                                % calculate PretoneOnly granger index 

                                granger_index = (element_1 - element_2)/(element_1 + element_2);
                                PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs = [PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs; granger_index];
                          
                            end

                        % NonSig Channel Pair Section 
                    
                        if pvalue_PFC_AC > 0.05
                        
                            % calculate PriorOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
        
                            element_1 = [];
                                for z = 1:length(Granger_Data.Granger_PriorOnly.grangerspctrm)
                                element_1 = [element_1; Granger_Data.Granger_PriorOnly.grangerspctrm(:,:, z)];
                                end
                            element_1 = element_1(element_1(:,1) ~= 0, 1);
                            element_1 = mean(element_1); % collapse across frequencies within band
        
                            PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs = [PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs, element_1];
        
                            % calculate PriorOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
        
                            element_2 = [];
                                for z = 1:length(Granger_Data.Granger_PriorOnly.grangerspctrm)
                                element_2 = [element_2; Granger_Data.Granger_PriorOnly.grangerspctrm(:,:, z)];
                                end
                            element_2 = element_2(element_2(:,2) ~= 0, 2);
                            element_2 = mean(element_2);
        
                            PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs = [PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs, element_2];
        
                            % calculate PriorOnly granger index 
        
                            granger_index = (element_1 - element_2)/(element_1 + element_2);
                            PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs = [PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs; granger_index];
        
                            % calculate PretoneOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
        
                            element_1 = [];
                                for z = 1:length(Granger_Data.Granger_PretoneOnly.grangerspctrm)
                                element_1 = [element_1; Granger_Data.Granger_PretoneOnly.grangerspctrm(:,:, z)];
                                end
                            element_1 = element_1(element_1(:,1) ~= 0, 1);
                            element_1 = mean(element_1);
        
                            PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs = [PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs, element_1];
        
                            % calculate PretoneOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
        
                            element_2 = [];
                                for z = 1:length(Granger_Data.Granger_PretoneOnly.grangerspctrm)
                                element_2 = [element_2; Granger_Data.Granger_PretoneOnly.grangerspctrm(:,:, z)];
                                end
                            element_2 = element_2(element_2(:,2) ~= 0, 2);
                            element_2 = mean(element_2);
        
                            PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs = [PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs, element_2];
        
                            % calculate PretoneOnly granger index 
        
                            granger_index = (element_1 - element_2)/(element_1 + element_2);
                            PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs = [PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs; granger_index];

                end

                clear Granger_Data
            end
        end
end
      
% %% plot the Granger Index for the P2A direction 
% 
%                     % plot the significant channels 
% 
%                     figure(1)
% 
%                     % scatter plot of prior Condition granger values 
%                     hold on; % Retain current plot when adding new plots
%                     % Scatter plot
%                     scatter(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs);
% 
%                     % Add labels
%                     xlabel('PriorOnly PFC-2-AC Granger Index');
%                     ylabel('PretoneOnly PFC-2-AC Granger Index');
%                     % Create the title string using sprintf
%                     titleString = sprintf('%s - %s: Granger Index of Significantly Modulated Channel Pairs', Animal, Frequency_Band);
%                     % Set the title of the figure
%                     title(titleString);
% 
%                     % Set axes to be equal
%                     axis equal;
% 
%                      % Set the origin to (0,0)
%                     % Determine the range for the axes
%                     xLimits = [min(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs), max(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs)];
%                     yLimits = [min(PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs), max(PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs)];
% 
%                     % Apply the determined limits
%                     xlim(xLimits);
%                     ylim(yLimits);
% 
%                     % Add a (1,1) linear line starting at the origin (0,0)
% 
%                     plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
% 
%                     % Perform a paired t-test
%                     [h, p] = ttest(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs);
% 
%                     % Display p-value on the plot
%                     text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
% 
%                     hold off; % Release hold on the current plot
% 
%                     figure(2)
% 
%                     % plot the non-significant channels 
% 
%                     % scatter plot of prior Condition granger values 
%                     hold on; % Retain current plot when adding new plots
%                     % Scatter plot
%                     scatter(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs, PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs);
% 
%                     % Add labels
%                     xlabel('PriorOnly PFC-2-AC Granger Index');
%                     ylabel('PretoneOnly PFC-2-AC Granger Index');
%                     % Create the title string using sprintf
%                     titleString = sprintf('%s - %s: Granger Index of Non-Significantly Modulated Channel Pairs', Animal, Frequency_Band);
%                     % Set the title of the figure
%                     title(titleString);
% 
%                     % Set axes to be equal
%                     axis equal;
% 
%                      % Set the origin to (0,0)
%                     % Determine the range for the axes
%                     xLimits = [min(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs), max(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs)];
%                     yLimits = [min(PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs), max(PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs)];
% 
%                     % Apply the determined limits
%                     xlim(xLimits);
%                     ylim(yLimits);
% 
%                     % Add a (1,1) linear line starting at the origin (0,0)
% 
%                     plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
% 
%                     % Perform a paired t-test
%                     [h, p] = ttest(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs, PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs);
% 
%                     % Display p-value on the plot
%                     text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
%                     hold off; % Release hold on the current plot
% 
%  %% Plot the Significant Granger Values for the P2A direction 
% 
%                     figure(3)
% 
%                     % within condition, plot the granger values
% 
%                     % scatter plot of prior Condition granger values 
%                     hold on; % Retain current plot when adding new plots
%                     % Scatter plot
%                     scatter(PriorOnly_Granger_Value_P2A_Sig_ChannelPairs, PriorOnly_Granger_Value_A2P_Sig_ChannelPairs);
% 
%                     ax = gca;
%                     ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
%                     ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
% 
%                     % Add labels
%                     xlabel('PriorOnly A2P Granger Value');
%                     ylabel('PriorOnly P2A Granger Value');
%                     % Create the title string using sprintf
%                     titleString = sprintf('%s - %s: Directionality of Significant Granger Values PriorOnly', Animal, Frequency_Band);
%                     % Set the title of the figure
%                     title(titleString);
% 
%                     % % Set axes to be equal
%                     % axis equal;
% 
%                      % Set the origin to (0,0)
%                     % Determine the range for the axes
%                     xLimits = [min(PriorOnly_Granger_Value_P2A_Sig_ChannelPairs), max(PriorOnly_Granger_Value_P2A_Sig_ChannelPairs)];
%                     yLimits = [min(PriorOnly_Granger_Value_A2P_Sig_ChannelPairs), max(PriorOnly_Granger_Value_A2P_Sig_ChannelPairs)];
% 
%                     % Apply the determined limits
%                     xlim(xLimits);
%                     ylim(yLimits);
% 
%                     plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
% 
%                     % Perform a paired t-test
%                     [h, p] = ttest(PriorOnly_Granger_Value_P2A_Sig_ChannelPairs, PriorOnly_Granger_Value_A2P_Sig_ChannelPairs);
% 
%                     % Display p-value on the plot
%                     text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
%                     hold off
% 
% 
%                     figure(4)
% 
%                     % within condition, plot the granger values
% 
%                     % scatter plot of prior Condition granger values 
%                     hold on; % Retain current plot when adding new plots
%                     % Scatter plot
%                     scatter(PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs);
% 
%                     ax = gca;
%                     ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
%                     ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
% 
%                     % Add labels
%                     xlabel('PretoneOnly A2P Granger Value');
%                     ylabel('PretoneOnly P2A Granger Value');
%                     % Create the title string using sprintf
%                     titleString = sprintf('%s - %s: Directionality of Significant Granger Values PretoneOnly', Animal, Frequency_Band);
%                     % Set the title of the figure
%                     title(titleString);
%                     % 
%                     % Set axes to be equal
%                     % axis equal;
% 
%                      % Set the origin to (0,0)
%                     % Determine the range for the axes
%                     xLimits = [min(PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs), max(PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs)];
%                     yLimits = [min(PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs), max(PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs)];
% 
%                     % Apply the determined limits
%                     xlim(xLimits);
%                     ylim(yLimits);
% 
%                     plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
% 
%                     % Perform a paired t-test
%                     [h, p] = ttest(PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs', PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs');
% 
%                     % Display p-value on the plot
%                     text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
%                     hold off
% 
% %% Plot the NonSignificant Granger Values for the P2A direction 
% 
%                     figure(5)
% 
%                     % within condition, plot the granger values
% 
%                     % scatter plot of prior Condition granger values 
%                     hold on; % Retain current plot when adding new plots
%                     % Scatter plot
%                     scatter(PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs, PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs);
% 
%                     ax = gca;
%                     ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
%                     ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
% 
%                     % Add labels
%                     xlabel('PriorOnly A2P Granger Value');
%                     ylabel('PriorOnly P2A Granger Value');
%                     % Create the title string using sprintf
%                     titleString = sprintf('%s - %s: Directionality of NonSignificant Granger Values PriorOnly', Animal, Frequency_Band);
%                     % Set the title of the figure
%                     title(titleString);
% 
%                     % % Set axes to be equal
%                     % axis equal;
% 
%                      % Set the origin to (0,0)
%                     % Determine the range for the axes
%                     xLimits = [min(PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs), max(PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs)];
%                     yLimits = [min(PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs), max(PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs)];
% 
%                     % Apply the determined limits
%                     xlim(xLimits);
%                     ylim(yLimits);
% 
%                     plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
% 
%                     % Perform a paired t-test
%                     [h, p] = ttest(PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs, PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs);
% 
%                     % Display p-value on the plot
%                     text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
%                     hold off
% 
%                     figure(6)
% 
%                     % within condition, plot the granger values
% 
%                     % scatter plot of prior Condition granger values 
%                     hold on; % Retain current plot when adding new plots
%                     % Scatter plot
%                     scatter(PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs, PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs);
% 
%                     ax = gca;
%                     ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
%                     ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
% 
%                     % Add labels
%                     xlabel('PretoneOnly A2P Granger Value');
%                     ylabel('PretoneOnly P2A Granger Value');
%                     % Create the title string using sprintf
%                     titleString = sprintf('%s - %s: Directionality of NonSignificant Granger Values PretoneOnly', Animal, Frequency_Band);
%                     % Set the title of the figure
%                     title(titleString);
%                     % 
%                     % Set axes to be equal
%                     %axis equal;
% 
%                      % Set the origin to (0,0)
%                     % Determine the range for the axes
%                     xLimits = [min(PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs), max(PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs)];
%                     yLimits = [min(PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs), max(PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs)];
% 
%                     % Apply the determined limits
%                     xlim(xLimits);
%                     ylim(yLimits);
% 
%                     plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
% 
%                     % Perform a paired t-test
%                     [h, p] = ttest(PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs', PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs');
% 
%                     % Display p-value on the plot
%                     text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
%                     hold off

%%

PriorOnly_Granger_Index_AllStongPairs         = vertcat(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs, PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs);
PretoneOnly_Granger_Index_AllStongPairs       = vertcat(PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs);

PriorOnly_Granger_P2A_Values_AllStongPairs    = horzcat(PriorOnly_Granger_Value_P2A_Sig_ChannelPairs, PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs);
PriorOnly_Granger_A2P_Values_AllStongPairs    = horzcat(PriorOnly_Granger_Value_A2P_Sig_ChannelPairs, PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs);
PretoneOnly_Granger_P2A_Values_AllStongPairs  = horzcat(PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs);
PretoneOnly_Granger_A2P_Values_AllStongPairs  = horzcat(PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs, PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs);

%% 



                    % figure(1)
                    % 
                    % % within condition, plot the granger values
                    % 
                    % % scatter plot of prior Condition granger values 
                    % hold on; % Retain current plot when adding new plots
                    % % Scatter plot
                    % scatter(PriorOnly_Granger_Index_AllStongPairs, PretoneOnly_Granger_Index_AllStongPairs);
                    % 
                    % ax = gca;
                    % ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
                    % ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis
                    % 
                    % % Add labels
                    % xlabel('PriorOnly A2P Granger Index');
                    % ylabel('PriorOnly P2A Granger Index');
                    % % Create the title string using sprintf
                    % titleString = sprintf('%s - %s: Directionality of A.S.P. Granger Index PriorOnly', Animal, Frequency_Band);
                    % % Set the title of the figure
                    % title(titleString);
                    % 
                    % % % Set axes to be equal
                    % % axis equal;
                    % 
                    %  % Set the origin to (0,0)
                    % % Determine the range for the axes
                    % xLimits = [min(PriorOnly_Granger_Index_AllStongPairs), max(PriorOnly_Granger_Index_AllStongPairs)];
                    % yLimits = [min(PretoneOnly_Granger_Index_AllStongPairs), max(PretoneOnly_Granger_Index_AllStongPairs)];
                    % 
                    % % Apply the determined limits
                    % xlim(xLimits);
                    % ylim(yLimits);
                    % 
                    % plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1
                    % 
                    % % Perform a paired t-test
                    % [h, p] = ttest(PriorOnly_Granger_Index_AllStongPairs, PretoneOnly_Granger_Index_AllStongPairs);
                    % 
                    % % Display p-value on the plot
                    % text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');
                    % hold off

%%

figure (1)

                    set(gcf, 'PaperUnits', 'inches'); % Change to 'centimeters' if desired
                    set(gcf, 'PaperPosition', [0 0 6 6]); % [left, bottom, width, height] in inches

                    % within condition, plot the granger values

                    % scatter plot of prior Condition granger values 
                    hold on; % Retain current plot when adding new plots
                    
                    % Scatter plot
                    scatter(PriorOnly_Granger_P2A_Values_AllStongPairs , PriorOnly_Granger_A2P_Values_AllStongPairs);

                    ax = gca;
                    ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
                    ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis

                    % Add labels
                    xlabel('LED Trials A2P Granger Value', 'FontSize', 14);
                    ylabel('LED Trials P2A Granger Value', 'FontSize', 14);
                    % Create the title string using sprintf
                    titleString = sprintf('%s - %s: Directionality of A.S.P. Granger Values PriorOnly', Animal, Frequency_Band);
                    % Set the title of the figure
                    %title(titleString);

                    % % Set axes to be equal
                    % axis equal;

                     % Set the origin to (0,0)
                    % Determine the range for the axes
                    xLimits = [min(PriorOnly_Granger_P2A_Values_AllStongPairs), max(PriorOnly_Granger_P2A_Values_AllStongPairs)];
                    yLimits = [min(PriorOnly_Granger_A2P_Values_AllStongPairs), max(PriorOnly_Granger_A2P_Values_AllStongPairs)];

                    % Apply the determined limits
                    xlim(xLimits);
                    ylim(yLimits);

                    plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Perform a paired t-test
                    [p, h, stats] = ranksum(PriorOnly_Granger_P2A_Values_AllStongPairs, PriorOnly_Granger_A2P_Values_AllStongPairs);

                 % Display p-value on the plot in the desired location
                    text(xLimits(2) * 0.95, yLimits(2) * 0.95, sprintf('p-value = %.4f', p), ...
                         'FontSize', 14, 'Color', 'b', 'HorizontalAlignment', 'right');
                                        
                    hold off

                    %%

figure(2)

                    set(gcf, 'PaperUnits', 'inches'); % Change to 'centimeters' if desired
                    set(gcf, 'PaperPosition', [0 0 6 6]); % [left, bottom, width, height] in inches

                    % within condition, plot the granger values

                    % scatter plot of prior Condition granger values 
                    hold on; % Retain current plot when adding new plots
                    % Scatter plot
                    scatter(PretoneOnly_Granger_P2A_Values_AllStongPairs , PretoneOnly_Granger_A2P_Values_AllStongPairs);

                    ax = gca;
                    ax.XAxis.Exponent = 0;  % Disable scientific notation for X-axis
                    ax.YAxis.Exponent = 0;  % Disable scientific notation for Y-axis

                    % Add labels
                    xlabel('Pretone Trials A2P Granger Value', 'FontSize', 14);
                    ylabel('Pretone Trials P2A Granger Value', 'FontSize', 14);
                    % Create the title string using sprintf
                    titleString = sprintf('%s - %s: Directionality of A.S.P. Granger Values PretoneOnly', Animal, Frequency_Band);
                    % Set the title of the figure
                    % title(titleString);

                    % % Set axes to be equal
                    % axis equal;

                     % Set the origin to (0,0)
                    % Determine the range for the axes
                    xLimits = [min(PretoneOnly_Granger_P2A_Values_AllStongPairs), max(PretoneOnly_Granger_P2A_Values_AllStongPairs)];
                    yLimits = [min(PretoneOnly_Granger_A2P_Values_AllStongPairs), max(PretoneOnly_Granger_A2P_Values_AllStongPairs)];

                    % Apply the determined limits
                    xlim(xLimits);
                    ylim(yLimits);

                    plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Perform a paired t-test
                    [p, h, stats] = ranksum(PretoneOnly_Granger_P2A_Values_AllStongPairs, PretoneOnly_Granger_A2P_Values_AllStongPairs);

                    % Display p-value on the plot in the desired location
                    text(xLimits(2) * 0.95, yLimits(2) * 0.95, sprintf('p-value = %.4f', p), ...
                         'FontSize', 14, 'Color', 'b', 'HorizontalAlignment', 'right');
                                        
                    hold off
