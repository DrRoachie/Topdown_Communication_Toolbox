
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
Epoch  = 'testToneOnset';

%% load histogram data arrays and max channel pairs
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

maxdir          = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_24_Analysis\XCorr_Histogram_Data\GRANGER_testTone_xcorr_histogram_data';
path_to_data    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_24_Analysis\XCorr_Histogram_Data\GRANGER_testTone_xcorr_histogram_data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;

%%

%rootdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis';
rootdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_24_Analysis';
sessions = dir(fullfile(rootdir, '19*'));

ChanPair_Array_PFC_AC = zeros(20,20);   % used to store passing p-value count for granger (bidirectional)
ChanPair_Array_AC_PFC = zeros(20,20);

% Initialize the arrays for comparing the Granger Index values across
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


%% pull the channel labels of the strongest connection

for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    fName = sprintf('pvalue_max_channels_%s_%s_PFC_AC.mat', Animal, freq_band);
    load(fullfile(maxdir, fName));
    max_results_PFC_AC = max_results;
    max_results_PFC_AC = max_results_PFC_AC(max_results_PFC_AC(:, 1)==1, :);
    clear max_results;

    fName = sprintf('pvalue_max_channels_%s_%s_AC_PFC.mat', Animal, freq_band);
    load(fullfile(maxdir, fName));
    max_results_AC_PFC = max_results;
    max_results_AC_PFC = max_results_AC_PFC(max_results_AC_PFC(:, 1)==1, :);
    clear max_results;

    max_chans_PFC_AC      = [max_results_PFC_AC(:,2) max_results_PFC_AC(:,3)];  % send chan is column 2 (send chan is PFC in PFC to AC direction, send chan is AC in AC to PFC direction)
    max_chans_AC_PFC      = [max_results_AC_PFC(:,2) max_results_AC_PFC(:,3)]; % IMPORTANT NOTE: whereas in the xcorr histogram code we need to flip the columns to adjust direction and extract the unique rows, granger needs to preserve the directionality

end

%% row-wise loops through the max-chan list and  pulls all files that in the top cut

% Assuming 'files' is a structure array, and each row contains a field 'name'
% Also assuming 'max_chans' is a 94x2 array with PFC and AC channel pairs

% Loop through the rows of the structure 'files'

if strcmp(Direction, 'P2A')

    matching_files = {};

    for i = 1:length(files)
        % Extract the file name from the 'name' field for the current row
        file_name = files(i).name;
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
        
            PFC_AC_matching_files = matching_files';

        end
    end
end


if strcmp(Direction, 'A2P')

    matching_files = {};

    for i = 1:length(files)
        % Extract the file name from the 'name' field for the current row
        file_name = files(i).name;
        
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
             
            AC_PFC_matching_files = matching_files';

        end
    end
end

%% pulls all the significant and non-significant channel pairs 

for i = 1:length(sessions) 

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};

        if exist(fullfile(rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
            
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), ['*' Frequency_Band '*Granger_data.mat']));
                
              if strcmp(Direction, 'P2A')
    
                    keep_idx = false(length(files), 1);  % Initially mark all entries as false (don't keep them)
    
                    for f = 1:length(files)                    % Loop through the 'files' structure and compare each file name with 'matching_files'

                    current_file_name = files(i).name;  % Extract the file name from the current 'files' structure
        
                    % Check if the current file name is in 'matching_files'
                    if any(strcmp(current_file_name, PFC_AC_matching_files))
                       keep_idx(i) = true; % If it is in 'matching_files', mark it to keep (set to true)
                    end

                    % Filter the 'files' structure to keep only those in 'matching_files'
                    filtered_files = files(keep_idx);

                    end 
              end           
 
        else
            continue;
        end

        if ~isempty(files)
           fprintf('Processing session %s: \n', num2str(RecDate));
        end
        
        p_count = 0;
        
        for k = 1:length(files) % runs if there are any processed channel pairs in given frequency band
        
            file_name = files(k).name;

            % pull channel pair info from file name
            tokens = regexp(file_name, '(.*?)_', 'tokens');
            if isempty(tokens)
                continue;
            end

            Send_Cort   = tokens{4}{1}; 
            Send_Num    = tokens{5}{1};
            Rec_Cort    = tokens{6}{1};
            Rec_Num     = tokens{7}{1};

            Send_Chan   = [Send_Cort '_' Send_Num];
            Rec_Chan    = [Rec_Cort '_' Rec_Num];

            % pull all present data files for channel pair
            if strcmp(Statistic, 'Granger')
                channel_files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), [extractBefore(file_name, 'Granger') '*']));
            end

            % check if Maris statistical test has been run on the channel pair yet (pair can have Coh_data/Granger_data file without
            % having statistical test results stored in that file since RunMaris_v3.m appends data from statistical test)
            
            if length(channel_files) ~= 15      % number of files for a channel pair is 15 after Maris runs
                continue;
            end
            
            % load p value of frequency band
            datadir = fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), file_name);

            if strcmp(Statistic, 'Granger')
    
                % load in p-value
                load(datadir,'pvalue_PFC_AC','pvalue_AC_PFC');
                pvalue_PFC_AC = pvalue_PFC_AC.(Frequency_Band);
                pvalue_AC_PFC = pvalue_AC_PFC.(Frequency_Band);

        
               if strcmp(Direction, 'P2A')

                    % ***************************************SIG P2A SECTION START******************************* %

                    if pvalue_PFC_AC <= 0.05

                        % basically you need 
                    load(datadir, 'Granger_PriorOnly');
                    load(datadir, 'Granger_PretoneOnly')

                    % calculate PriorOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              

                    element_1 = [];
                    for z = 1:length(Granger_PriorOnly.grangerspctrm)
                    element_1 = [element_1; Granger_PriorOnly.grangerspctrm(:,:, z)];
                    end

                    element_1 = element_1(element_1(:,1) ~= 0, 1);
                    element_1 = mean(element_1); % collapse across frequencies within band

                    PriorOnly_Granger_Value_P2A_Sig_ChannelPairs = [PriorOnly_Granger_Value_P2A_Sig_ChannelPairs, element_1];

                    % calculate PriorOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    

                    element_2 = [];
                    for z = 1:length(Granger_PriorOnly.grangerspctrm)
                    element_2 = [element_2; Granger_PriorOnly.grangerspctrm(:,:, z)];
                    end

                    element_2 = element_2(element_2(:,2) ~= 0, 2);
                    element_2 = mean(element_2);

                    PriorOnly_Granger_Value_A2P_Sig_ChannelPairs = [PriorOnly_Granger_Value_A2P_Sig_ChannelPairs, element_2];

                    % calculate PriorOnly granger index 

                    granger_index = (element_1 - element_2)/(element_1 + element_2);
                    PriorOnly_Granger_Index_P2A_Sig_ChannelPairs = [PriorOnly_Granger_Index_P2A_Sig_ChannelPairs; granger_index];

                    
                    % calculate PretoneOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              

                    element_1 = [];
                    for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                    element_1 = [element_1; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                    end

                    element_1 = element_1(element_1(:,1) ~= 0, 1);
                    element_1 = mean(element_1);

                    PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs = [PretoneOnly_Granger_Value_P2A_Sig_ChannelPairs, element_1];

                    % calculate PretoneOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    

                    element_2 = [];
                    for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                    element_2 = [element_2; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                    end

                    element_2 = element_2(element_2(:,2) ~= 0, 2);
                    element_2 = mean(element_2);

                    PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs = [PretoneOnly_Granger_Value_A2P_Sig_ChannelPairs, element_2];

                    % calculate PretoneOnly granger index 

                    granger_index = (element_1 - element_2)/(element_1 + element_2);
                    PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs = [PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs; granger_index];

                    end

                    % ***************************************NON-SIG P2A SECTION START******************************* %

                    % get NON-SIG channels  granger index
                   
                    if pvalue_PFC_AC >= 0.05

                    load(datadir,'Granger_PriorOnly');
                    load(datadir, 'Granger_PretoneOnly')
    
                     % calculate NON-SIG PriorOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
    
                     element_1 = [];
                     for z = 1:length(Granger_PriorOnly.grangerspctrm)
                     element_1 = [element_1; Granger_PriorOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_1 = element_1(element_1(:,1) ~= 0, 1);
                     element_1 = mean(element_1); % collapse across frequencies within band
    
                     PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs      = [PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs, element_1];

               
                     % calculate NON-SIG PriorOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
    
                     element_2 = [];
                     for z = 1:length(Granger_PriorOnly.grangerspctrm)
                     element_2 = [element_2; Granger_PriorOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_2 = element_2(element_2(:,2) ~= 0, 2);
                     element_2 = mean(element_2);
                        
                     PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs      = [PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs, element_2];
    
                     % calculate granger index 
    
                     granger_index = (element_1 - element_2)/(element_1 + element_2);
                        
                     PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs = [PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs; granger_index];
    
                     % calculate PretoneOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
    
                     element_1 = [];
                     for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                     element_1 = [element_1; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_1 = element_1(element_1(:,1) ~= 0, 1);
                     element_1 = mean(element_1);
    
                     PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs      = [PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs, element_1];
    
                     % calculate PretoneOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
    
                     element_2 = [];
                     for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                     element_2 = [element_2; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_2 = element_2(element_2(:,2) ~= 0, 2);
                     element_2 = mean(element_2);
    
                     PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs      = [PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs, element_2];
    
                     % calculate granger index 
    
                     granger_index = (element_1 - element_2)/(element_1 + element_2);
                     PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs = [PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs; granger_index];
        
                     end
               end
    
               if strcmp(Direction, 'A2P')
    
                     % ***************************************SIG A2P SECTION START******************************* %

                     if pvalue_AC_PFC <= 0.05
    
                     load(datadir,'Granger_PriorOnly');
                     load(datadir, 'Granger_PretoneOnly')
    
                     % calculate PriorOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
    
                     element_1 = [];
                     for z = 1:length(Granger_PriorOnly.grangerspctrm)
                     element_1 = [element_1; Granger_PriorOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_1 = element_1(element_1(:,1) ~= 0, 1);
                     element_1 = mean(element_1); % collapse across frequencies within band
    
                     % calculate PriorOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
    
                     element_2 = [];
                     for z = 1:length(Granger_PriorOnly.grangerspctrm)
                     element_2 = [element_2; Granger_PriorOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_2 = element_2(element_2(:,2) ~= 0, 2);
                     element_2 = mean(element_2);
    
                     % calculate granger index 
    
                     granger_index = (element_1 - element_2)/(element_1 + element_2);
                     PriorOnly_Granger_Index_A2P_Sig_ChannelPairs = [PriorOnly_Granger_Index_A2P_Sig_ChannelPairs; granger_index];
    
                     % calculate PretoneOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
    
                     element_1 = [];
                     for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                     element_1 = [element_1; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_1 = element_1(element_1(:,1) ~= 0, 1);
                     element_1 = mean(element_1);
    
                     % calculate PretoneOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
    
                     element_2 = [];
                     for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                     element_2 = [element_2; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_2 = element_2(element_2(:,2) ~= 0, 2);
                     element_2 = mean(element_2);
    
                     % calculate granger index 
    
                     granger_index = (element_1 - element_2)/(element_1 + element_2);
                     PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs = [PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs; granger_index];

                     end
                    
                     % ***************************************NON-SIG A2P SECTION START******************************* %

                      % get NON-SIG channels  granger index
                   
                    if pvalue_AC_PFC >= 0.05

                    load(datadir,'Granger_PriorOnly');
                    load(datadir, 'Granger_PretoneOnly')
    
                     % calculate NON-SIG PriorOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
    
                     element_1 = [];
                     for z = 1:length(Granger_PriorOnly.grangerspctrm)
                     element_1 = [element_1; Granger_PriorOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_1 = element_1(element_1(:,1) ~= 0, 1);
                     element_1 = mean(element_1); % collapse across frequencies within band
    
                     PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs      = [PriorOnly_Granger_Value_P2A_NonSig_ChannelPairs, element_1];

               
                     % calculate NON-SIG PriorOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
    
                     element_2 = [];
                     for z = 1:length(Granger_PriorOnly.grangerspctrm)
                     element_2 = [element_2; Granger_PriorOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_2 = element_2(element_2(:,2) ~= 0, 2);
                     element_2 = mean(element_2);
                        
                     PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs      = [PriorOnly_Granger_Value_A2P_NonSig_ChannelPairs, element_2];
    
                     % calculate granger index 
    
                     granger_index = (element_1 - element_2)/(element_1 + element_2);
                        
                     PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs = [PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs; granger_index];
    
                     % calculate PretoneOnly granger value averaged across frequency band for the PFC-to-AC direction, results in single granger value              
    
                     element_1 = [];
                     for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                     element_1 = [element_1; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_1 = element_1(element_1(:,1) ~= 0, 1);
                     element_1 = mean(element_1);
    
                     PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs      = [PretoneOnly_Granger_Value_P2A_NonSig_ChannelPairs, element_1];
    
                     % calculate PretoneOnly granger value averaged across frequency band for the AC-to-PFC direction, results in single granger value    
    
                     element_2 = [];
                     for z = 1:length(Granger_PretoneOnly.grangerspctrm)
                     element_2 = [element_2; Granger_PretoneOnly.grangerspctrm(:,:, z)];
                     end
    
                     element_2 = element_2(element_2(:,2) ~= 0, 2);
                     element_2 = mean(element_2);
    
                     PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs      = [PretoneOnly_Granger_Value_A2P_NonSig_ChannelPairs, element_2];
    
                     % calculate granger index 
    
                     granger_index = (element_1 - element_2)/(element_1 + element_2);
                     PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs = [PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs; granger_index];
    
                    end
                     
                 end

            end
        end
               
    end

end


                    
            %% plot the Granger Index for the P2A direction 

 if strcmp(Direction, 'P2A')

                    % plot the significant channels 

                    figure(1)

                    % scatter plot of prior Condition granger values 
                    hold on; % Retain current plot when adding new plots
                    % Scatter plot
                    scatter(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs);


                    % Add labels
                    xlabel('PriorOnly PFC-2-AC Granger Index');
                    ylabel('PretoneOnly PFC-2-AC Granger Index');
                    % Create the title string using sprintf
                    titleString = sprintf('%s - %s: Granger Index of Significantly Modulated Channel Pairs', Animal, Frequency_Band);
                    % Set the title of the figure
                    title(titleString);

                    % Set axes to be equal
                    axis equal;

                     % Set the origin to (0,0)
                    % Determine the range for the axes
                    xLimits = [min(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs), max(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs)];
                    yLimits = [min(PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs), max(PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs)];
                    
                    % Apply the determined limits
                    xlim(xLimits);
                    ylim(yLimits);

                    % % Add a (1,1) linear line starting at the origin (0,0)
                    % hold on; % Retain current plot when adding new plots
                    % plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Add a (1,1) linear line starting at the origin (0,0)

                    plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Perform a paired t-test
                    [h, p] = ttest(PriorOnly_Granger_Index_P2A_Sig_ChannelPairs, PretoneOnly_Granger_Index_P2A_Sig_ChannelPairs);

                    % Display p-value on the plot
                    text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');

                    hold off; % Release hold on the current plot

                    figure(2)

                    % plot the non-significant channels 

                    % scatter plot of prior Condition granger values 
                    hold on; % Retain current plot when adding new plots
                    % Scatter plot
                    scatter(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs, PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs);


                    % Add labels
                    xlabel('PriorOnly PFC-2-AC Granger Index');
                    ylabel('PretoneOnly PFC-2-AC Granger Index');
                    % Create the title string using sprintf
                    titleString = sprintf('%s - %s: Granger Index of Non-Significantly Modulated Channel Pairs', Animal, Frequency_Band);
                    % Set the title of the figure
                    title(titleString);

                    % Set axes to be equal
                    axis equal;

                     % Set the origin to (0,0)
                    % Determine the range for the axes
                    xLimits = [min(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs), max(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs)];
                    yLimits = [min(PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs), max(PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs)];
                    
                    % Apply the determined limits
                    xlim(xLimits);
                    ylim(yLimits);


                    % % Add a (1,1) linear line starting at the origin (0,0)
                    % hold on; % Retain current plot when adding new plots
                    % plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Add a (1,1) linear line starting at the origin (0,0)

                    plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Perform a paired t-test
                    [h, p] = ttest(PriorOnly_Granger_Index_P2A_NonSig_ChannelPairs, PretoneOnly_Granger_Index_P2A_NonSig_ChannelPairs);

                    % Display p-value on the plot
                    text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');


 end

 %% plot the Granger Index for the A2P direction 


 if strcmp(Direction, 'A2P')

                    % plot the significant channels 

                    figure(1)

                    % scatter plot of prior Condition granger values 
                    hold on; % Retain current plot when adding new plots
                    % Scatter plot
                    scatter(PriorOnly_Granger_Index_A2P_Sig_ChannelPairs, PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs);


                    % Add labels
                    xlabel('PriorOnly AC-2-PFC Granger Index');
                    ylabel('PretoneOnly AC-2-PFC Granger Index');
                   % Create the title string using sprintf
                    titleString = sprintf('%s - %s: Granger Index of Significantly Modulated Channel Pairs', Animal, Frequency_Band);
                    % Set the title of the figure
                    title(titleString);

                    % Set axes to be equal
                    axis equal;

                     % Set the origin to (0,0)
                    % Determine the range for the axes
                    xLimits = [min(PriorOnly_Granger_Index_A2P_Sig_ChannelPairs), max(PriorOnly_Granger_Index_A2P_Sig_ChannelPairs)];
                    yLimits = [min(PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs), max(PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs)];
                    
                    % Apply the determined limits
                    xlim(xLimits);
                    ylim(yLimits);

                    % % Add a (1,1) linear line starting at the origin (0,0)
                    % hold on; % Retain current plot when adding new plots
                    % plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Add a (1,1) linear line starting at the origin (0,0)

                    plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Perform a paired t-test
                    [h, p] = ttest(PriorOnly_Granger_Index_A2P_Sig_ChannelPairs, PretoneOnly_Granger_Index_A2P_Sig_ChannelPairs);

                    % Display p-value on the plot
                    text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');

                    hold off; % Release hold on the current plot

                    figure(2)

                    % plot the non-significant channels 

                    % scatter plot of prior Condition granger values 
                    hold on; % Retain current plot when adding new plots
                    
                    % Scatter plot
                    scatter(PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs, PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs);

                    % Add labels
                    xlabel('PriorOnly AC -2-PFC Granger Index');
                    ylabel('PretoneOnly AC-2-PFC Granger Index');
                    % Create the title string using sprintf
                    titleString = sprintf('%s - %s: Granger Index of Non-Significantly Modulated Channel Pairs', Animal, Frequency_Band);
                    % Set the title of the figure
                    title(titleString);

                    % Set axes to be equal
                    axis equal;

                     % Set the origin to (0,0)
                    % Determine the range for the axes
                    xLimits = [min(PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs), max(PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs)];
                    yLimits = [min(PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs), max(PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs)];
                    
                    % Apply the determined limits
                    xlim(xLimits);
                    ylim(yLimits);


                    % % Add a (1,1) linear line starting at the origin (0,0)
                    % hold on; % Retain current plot when adding new plots
                    % plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Add a (1,1) linear line starting at the origin (0,0)

                    plot([min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], [min([xLimits(1), yLimits(1)]), max([xLimits(2), yLimits(2)])], 'r--'); % Red dashed line with slope 1

                    % Perform a paired t-test
                    [h, p] = ttest(PriorOnly_Granger_Index_A2P_NonSig_ChannelPairs, PretoneOnly_Granger_Index_A2P_NonSig_ChannelPairs);

                    % Display p-value on the plot
                    text(max(xlim)*0.05, max(ylim)*0.95, sprintf('p-value = %.4f', p), 'FontSize', 12, 'Color', 'b');

 end

 %% Plot granger values for prioronly Condition 


%prior_only 
