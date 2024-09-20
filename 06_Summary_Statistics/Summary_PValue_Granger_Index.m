
% This code takes the granger spectra from channel pairs that are
% significantly modulated by context, and tests the correlation between
% directionalities (PFC-to-AC) or (AC-to-PFC) within context condition (OnlyPrior or OnlyPretone) or
% between context conditions within a single direction. 

%Note: there is also code that calculates a simple granger index and uses
%that as a metric to compare across conditions. 

%This code also carries out a paired t-test between the two conditions

%% Define data

Frequency_Band  = 'theta';               % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'Coherence';           % 'Granger' or 'Coherence'
animals         = {'MrM'};         % 'MrCassius' and/or 'MrM'
Direction       = 'P2A';                 % 'P2A' or 'A2P'      

Animal  = 'MrM';       % 'both' or 'MrCassius' or 'MrM'
freq_band_list  = {'theta'};



%rootdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis';
rootdir = 'D:\2024_07_24_Analysis';
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


%% Pull all the necessary granger values and calculate granger index 

for i = 1:length(sessions) 

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        
        Epoch  = 'testToneOnset';

        if exist(fullfile(rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
            
            if strcmp(Statistic, 'Granger')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), ['*' Frequency_Band '*Granger_data.mat']));
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

                    load(datadir,'Granger_PriorOnly');
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
