%% Define data

Frequency_Band  = 'alpha';                     % Options: 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'PSI';                       % Options: 'Granger' or 'Coherence'
animals         = {'MrCassius', 'MrM'};        % Options: 'MrCassius' and/or 'MrM'
Epoch           = 'preCueOnset';             % Options: 'testToneOnset', 'preCueOnset', 'moveOnset'
SNR             = 'AllSNR';                    % Options: 'HighSNR' or 'LowSNR'or 'AllSNR'
Behavior        = 'Wrong';                     % Options: 'Correct' or 'Wrong'
Condition       = 'PriorOnly';               % Options: 'PriorOnly' or 'PretoneOnly'
Congruency_1    = 'Congruent';                 % Options: 'Congruent' or 'Incongruent'
Congruency_2    = 'Incongruent';               % Options: 'Congruent' or 'Incongruent
Direction_1     = 'PFC_2_AC';                  % Options: 'PFC_2_AC' or 'AC_2_PFC'
Direction_2     = 'AC_2_PFC';                  % Options: 'PFC_2_AC' or 'AC_2_PFC'

% Select PSI data directory
% as of 1/14/25 the most important input is the directory because it
% determines what two conditions are stored in the file 

rootdir  = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Analysis_Alpha\PriorOnly\PriorOnly_PreCue_Wrong_Congruency_PSI';
savedir  = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Histogram_Summary\Alpha';                                        % path to the parent directory where all new data will be stored (save structure: RecDate >> Animal >> Epoch >> all figures/files)
sessions = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Condition_1_PSI_grid = cell(20, 20); 
Condition_2_PSI_grid = cell(20, 20);

%% Loop through all sessions and sort PSI values by channel-pair 

% Preallocate cell arrays for storing results
num_sessions   = length(sessions);
monte_Results  = cell(num_sessions, 1);
Wilcox_Results = cell(num_sessions, 1);

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};

        if exist(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct'), 'dir')
            if strcmp(Statistic, 'PSI')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', ['*' Frequency_Band '*PSI_data.mat']));
            end
        else
        continue;
        end
        
        if ~isempty(files)
            fprintf('Processing session %s: ', num2str(RecDate));
        end

            for k = 1:length(files) % runs if there are any processed channel pairs in given frequency band
        
            file_name = files(k).name;

            % pull this sessions, PSI data 
            if strcmp(Statistic, 'PSI')
            channel_files = dir(fullfile(rootdir, RecDate, Animal, 'preCue', 'Correct', [extractBefore(file_name, 'PSI') '*']));
            end

            % 
            datadir = fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', file_name);
            load(datadir);
           
            dataMatrix_1  = PSI_Condition_1.psispctrm; 
            labels_1      = PSI_Condition_1.labelcmb;         
            dataMatrix_2  = PSI_Condition_2.psispctrm; 
            labels_2      = PSI_Condition_2.labelcmb;         
    
           % Populate Condition_1_PSI_grid with data
        
           for rowIdx_1 = 1:size(dataMatrix_1, 1)

            % Extract channel pair from label
            pfcLabel = labels_1{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = labels_1{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
            
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
        
            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx = find(strcmp(channels, strrep(acChannel, 'AC_', '')));
        
                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Condition_1_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Condition_1_PSI_grid{pfcIdx, acIdx} = dataMatrix_1(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Condition_1_PSI_grid{pfcIdx, acIdx} = [Condition_1_PSI_grid{pfcIdx, acIdx}; dataMatrix_1(rowIdx_1, :)];
                    end

                end
           end

           % Populate Condition_2_PSI_grid with data 

           for rowIdx_2 = 1:size(dataMatrix_2, 1)

            % Extract channel pair from label
            pfcLabel = labels_2{rowIdx_2, 1}; % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = labels_2{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')
            
            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel
        
            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));
        
                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Condition_2_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Condition_2_PSI_grid{pfcIdx, acIdx} = dataMatrix_2(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Condition_2_PSI_grid{pfcIdx, acIdx} = [Condition_2_PSI_grid{pfcIdx, acIdx}; dataMatrix_2(rowIdx_2, :)];
                    end                 
                end
           end

           end
    end
end


%% Concatonate across channels 

% Initialize an empty array to store all trials
Condition_1_all_trials = [];

% Iterate through each cell in Condition_1_PSI_grid
for i = 1:numel(Condition_1_PSI_grid)
    % Concatenate the trials from the current cell
    Condition_1_all_trials = [Condition_1_all_trials; Condition_1_PSI_grid{i}];
end

% Initialize an empty array to store all trials
Condition_2_all_trials = [];

% Iterate through each cell in Condition_1_PSI_grid
for b = 1:numel(Condition_2_PSI_grid)
    % Concatenate the trials from the current cell
    Condition_2_all_trials = [Condition_2_all_trials; Condition_2_PSI_grid{i}];
end



%% Band-wise collect all values for each channel pair, and collapse into column

if strcmp(Frequency_Band,'theta') == 1

% Extract columns 5 to 8 and reshape into a single column
Condition_1 = Condition_1_all_trials(:, 5:8);
Condition_1 = Condition_1(:);

Condition_2 = Condition_2_all_trials(:, 5:8);
Condition_2 = Condition_2(:);

Condition_1_PFC_2_AC = Condition_1(Condition_1 > 0);
Condition_1_AC_2_PFC = Condition_1(Condition_1 < 0);

Condition_2_PFC_2_AC = Condition_2(Condition_2 > 0);
Condition_2_AC_2_PFC = Condition_2(Condition_2 < 0);

% Define the file names based on the variables, including the directory
fn_1 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_1, Direction_1, Frequency_Band));
fn_2 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_1, Direction_2, Frequency_Band));
fn_3 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_2, Direction_1, Frequency_Band));
fn_4 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_2, Direction_2, Frequency_Band));

% Save the variables into separate .mat files in the specified directory
save(fn_1, 'Condition_1_PFC_2_AC');  % Replace 'Variable1' with your actual variable name
save(fn_2, 'Condition_1_AC_2_PFC');  % Replace 'Variable2' with your actual variable name
save(fn_3, 'Condition_2_PFC_2_AC');  % Replace 'Variable3' with your actual variable name
save(fn_4, 'Condition_2_AC_2_PFC');  % Replace 'Variable4' with your actual variable name

end

if strcmp(Frequency_Band,'alpha') == 1

% Extract columns 5 to 8 and reshape into a single column
Condition_1 = Condition_1_all_trials(:, 9:14);
Condition_1 = Condition_1(:);

Condition_2 = Condition_2_all_trials(:, 9:14);
Condition_2 = Condition_2(:);

Condition_1_PFC_2_AC = Condition_1(Condition_1 > 0);
Condition_1_AC_2_PFC = Condition_1(Condition_1 < 0);

Condition_2_PFC_2_AC = Condition_2(Condition_2 > 0);
Condition_2_AC_2_PFC = Condition_2(Condition_2 < 0);

% Define the file names based on the variables, including the directory
fn_1 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_1, Direction_1, Frequency_Band));
fn_2 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_1, Direction_2, Frequency_Band));
fn_3 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_2, Direction_1, Frequency_Band));
fn_4 = fullfile(savedir, sprintf('%s_%s_%s_%s_%s_%s_%s.mat', Condition, Epoch, SNR, Behavior, Congruency_2, Direction_2, Frequency_Band));

% Save the variables into separate .mat files in the specified directory
save(fn_1, 'Condition_1_PFC_2_AC');  % Replace 'Variable1' with your actual variable name
save(fn_2, 'Condition_1_AC_2_PFC');  % Replace 'Variable2' with your actual variable name
save(fn_3, 'Condition_2_PFC_2_AC');  % Replace 'Variable3' with your actual variable name
save(fn_4, 'Condition_2_AC_2_PFC');  % Replace 'Variable4' with your actual variable name

end

if strcmp(Frequency_Band,'beta') == 1
% Extract columns 5 to 8 and reshape into a single column
Condition_1_beta = Condition_1_all_trials(:, 15:30);
Condition_1_beta = Condition_1_beta(:);

Condition_2_beta = Condition_2_all_trials(:, 15:30);
Condition_2_beta = Condition_2_beta(:);
end

%% Comparison: Prior versus Pretone 

% Inefficiency_Flag.

Frequency_Band = 'alpha';

if strcmp(Frequency_Band,'theta') == 1

datadir = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Histogram_Summary';

Set_1 = '\PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC_theta';
Set_2 = '\PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC_theta';
Set_3 = '\PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC_theta';
Set_4 = '\PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC_theta';

Set_5 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC_theta';
Set_6 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC_theta';
Set_7 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC_theta';
Set_8 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC_theta';

Set_9 =  '\PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC_theta';
Set_10 = '\PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC_theta';
Set_11 = '\PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC_theta';
Set_12 = '\PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC_theta';

Set_13 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC_theta';
Set_14 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC_theta';
Set_15 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC_theta';
Set_16 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC_theta';


Set_17 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC_theta';
Set_18 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC_theta';
Set_19 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC_theta';
Set_20 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC_theta';

Set_21 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC_theta';
Set_22 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC_theta';
Set_23 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC_theta';
Set_24 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC_theta';

Set_25 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC_theta';
Set_26 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC_theta';
Set_27 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC_theta';
Set_28 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC_theta';

Set_29 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC_theta';
Set_30 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC_theta';
Set_31 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC_theta';
Set_32 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC_theta';

Set_33 = '\PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC_theta';
Set_34 = '\PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC_theta';
Set_35 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC_theta';
Set_36 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC_theta';

Set_37 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC_theta';
Set_38 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC_theta';
Set_39 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC_theta';
Set_40 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC_theta';

Set_41 = '\PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC_theta';
Set_42 = '\PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC_theta';
Set_43 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC_theta';
Set_44 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC_theta';

Set_45 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC_theta';
Set_46 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC_theta';
Set_47 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC_theta';
Set_48 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC_theta';

data_fn_1 = sprintf('%s%s.mat', datadir, Set_1);
data_fn_2 = sprintf('%s%s.mat', datadir, Set_2);
data_fn_3 = sprintf('%s%s.mat', datadir, Set_3);
data_fn_4 = sprintf('%s%s.mat', datadir, Set_4);
data_fn_5 = sprintf('%s%s.mat', datadir, Set_5);
data_fn_6 = sprintf('%s%s.mat', datadir, Set_6);
data_fn_7 = sprintf('%s%s.mat', datadir, Set_7);
data_fn_8 = sprintf('%s%s.mat', datadir, Set_8);
data_fn_9 = sprintf('%s%s.mat', datadir, Set_9);
data_fn_10 = sprintf('%s%s.mat', datadir, Set_10);
data_fn_11 = sprintf('%s%s.mat', datadir, Set_11);
data_fn_12 = sprintf('%s%s.mat', datadir, Set_12);
data_fn_13 = sprintf('%s%s.mat', datadir, Set_13);
data_fn_14 = sprintf('%s%s.mat', datadir, Set_14);
data_fn_15 = sprintf('%s%s.mat', datadir, Set_15);
data_fn_16 = sprintf('%s%s.mat', datadir, Set_16);
data_fn_17 = sprintf('%s%s.mat', datadir, Set_17);
data_fn_18 = sprintf('%s%s.mat', datadir, Set_18);
data_fn_19 = sprintf('%s%s.mat', datadir, Set_19);
data_fn_20 = sprintf('%s%s.mat', datadir, Set_20);
data_fn_21 = sprintf('%s%s.mat', datadir, Set_21);
data_fn_22 = sprintf('%s%s.mat', datadir, Set_22);
data_fn_23 = sprintf('%s%s.mat', datadir, Set_23);
data_fn_24 = sprintf('%s%s.mat', datadir, Set_24);
data_fn_25 = sprintf('%s%s.mat', datadir, Set_25);
data_fn_26 = sprintf('%s%s.mat', datadir, Set_26);
data_fn_27 = sprintf('%s%s.mat', datadir, Set_27);
data_fn_28 = sprintf('%s%s.mat', datadir, Set_28);
data_fn_29 = sprintf('%s%s.mat', datadir, Set_29);
data_fn_30 = sprintf('%s%s.mat', datadir, Set_30);
data_fn_31 = sprintf('%s%s.mat', datadir, Set_31);
data_fn_32 = sprintf('%s%s.mat', datadir, Set_32);
data_fn_33 = sprintf('%s%s.mat', datadir, Set_33);
data_fn_34 = sprintf('%s%s.mat', datadir, Set_34);
data_fn_35 = sprintf('%s%s.mat', datadir, Set_35);
data_fn_36 = sprintf('%s%s.mat', datadir, Set_36);
data_fn_37 = sprintf('%s%s.mat', datadir, Set_37);
data_fn_38 = sprintf('%s%s.mat', datadir, Set_38);
data_fn_39 = sprintf('%s%s.mat', datadir, Set_39);
data_fn_40 = sprintf('%s%s.mat', datadir, Set_40);
data_fn_41 = sprintf('%s%s.mat', datadir, Set_41);
data_fn_42 = sprintf('%s%s.mat', datadir, Set_42);
data_fn_43 = sprintf('%s%s.mat', datadir, Set_43);
data_fn_44 = sprintf('%s%s.mat', datadir, Set_44);
data_fn_45 = sprintf('%s%s.mat', datadir, Set_45);
data_fn_46 = sprintf('%s%s.mat', datadir, Set_46);
data_fn_47 = sprintf('%s%s.mat', datadir, Set_47);
data_fn_48 = sprintf('%s%s.mat', datadir, Set_48);

load(data_fn_1);
PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_2);
PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_3);
PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_4);
PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_5);
PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_6);
PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_7);
PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_8);
PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_9);
PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_10);
PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_11);
PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_12);
PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_13);
PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_14);
PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_15);
PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_16);
PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_17);
PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_18);
PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_19);
PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_20);
PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_21);
PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_22);
PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_23);
PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_24);
PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_25);
PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_26);
PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_27);
PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_28);
PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_29);
PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_30);
PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_31);
PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_32);
PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_33);
PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_34);
PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_35);
PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_36);
PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_37);
PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_38);
PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_39);
PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_40);
PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_41);
PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_42);
PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_43);
PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_44);
PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_45);
PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_46);
PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_47);
PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_48);
PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

end



if strcmp(Frequency_Band,'alpha') == 1

datadir = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Histogram_Summary\Alpha';

Set_1 = '\PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC_alpha';
Set_2 = '\PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC_alpha';
Set_3 = '\PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC_alpha';
Set_4 = '\PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC_alpha';

Set_5 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC_alpha';
Set_6 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC_alpha';
Set_7 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC_alpha';
Set_8 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC_alpha';

Set_9 =  '\PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC_alpha';
Set_10 = '\PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC_alpha';
Set_11 = '\PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC_alpha';
Set_12 = '\PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC_alpha';

Set_13 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC_alpha';
Set_14 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC_alpha';
Set_15 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC_alpha';
Set_16 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC_alpha';


Set_17 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC_alpha';
Set_18 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC_alpha';
Set_19 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC_alpha';
Set_20 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC_alpha';

Set_21 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC_alpha';
Set_22 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC_alpha';
Set_23 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC_alpha';
Set_24 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC_alpha';

Set_25 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC_alpha';
Set_26 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC_alpha';
Set_27 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC_alpha';
Set_28 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC_alpha';

Set_29 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC_alpha';
Set_30 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC_alpha';
Set_31 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC_alpha';
Set_32 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC_alpha';

Set_33 = '\PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC_alpha';
Set_34 = '\PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC_alpha';
Set_35 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC_alpha';
Set_36 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC_alpha';

Set_37 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC_alpha';
Set_38 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC_alpha';
Set_39 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC_alpha';
Set_40 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC_alpha';

Set_41 = '\PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC_alpha';
Set_42 = '\PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC_alpha';
Set_43 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC_alpha';
Set_44 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC_alpha';

Set_45 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC_alpha';
Set_46 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC_alpha';
Set_47 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC_alpha';
Set_48 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC_alpha';

data_fn_1 = sprintf('%s%s.mat', datadir, Set_1);
data_fn_2 = sprintf('%s%s.mat', datadir, Set_2);
data_fn_3 = sprintf('%s%s.mat', datadir, Set_3);
data_fn_4 = sprintf('%s%s.mat', datadir, Set_4);
data_fn_5 = sprintf('%s%s.mat', datadir, Set_5);
data_fn_6 = sprintf('%s%s.mat', datadir, Set_6);
data_fn_7 = sprintf('%s%s.mat', datadir, Set_7);
data_fn_8 = sprintf('%s%s.mat', datadir, Set_8);
data_fn_9 = sprintf('%s%s.mat', datadir, Set_9);
data_fn_10 = sprintf('%s%s.mat', datadir, Set_10);
data_fn_11 = sprintf('%s%s.mat', datadir, Set_11);
data_fn_12 = sprintf('%s%s.mat', datadir, Set_12);
data_fn_13 = sprintf('%s%s.mat', datadir, Set_13);
data_fn_14 = sprintf('%s%s.mat', datadir, Set_14);
data_fn_15 = sprintf('%s%s.mat', datadir, Set_15);
data_fn_16 = sprintf('%s%s.mat', datadir, Set_16);
data_fn_17 = sprintf('%s%s.mat', datadir, Set_17);
data_fn_18 = sprintf('%s%s.mat', datadir, Set_18);
data_fn_19 = sprintf('%s%s.mat', datadir, Set_19);
data_fn_20 = sprintf('%s%s.mat', datadir, Set_20);
data_fn_21 = sprintf('%s%s.mat', datadir, Set_21);
data_fn_22 = sprintf('%s%s.mat', datadir, Set_22);
data_fn_23 = sprintf('%s%s.mat', datadir, Set_23);
data_fn_24 = sprintf('%s%s.mat', datadir, Set_24);
data_fn_25 = sprintf('%s%s.mat', datadir, Set_25);
data_fn_26 = sprintf('%s%s.mat', datadir, Set_26);
data_fn_27 = sprintf('%s%s.mat', datadir, Set_27);
data_fn_28 = sprintf('%s%s.mat', datadir, Set_28);
data_fn_29 = sprintf('%s%s.mat', datadir, Set_29);
data_fn_30 = sprintf('%s%s.mat', datadir, Set_30);
data_fn_31 = sprintf('%s%s.mat', datadir, Set_31);
data_fn_32 = sprintf('%s%s.mat', datadir, Set_32);
data_fn_33 = sprintf('%s%s.mat', datadir, Set_33);
data_fn_34 = sprintf('%s%s.mat', datadir, Set_34);
data_fn_35 = sprintf('%s%s.mat', datadir, Set_35);
data_fn_36 = sprintf('%s%s.mat', datadir, Set_36);
data_fn_37 = sprintf('%s%s.mat', datadir, Set_37);
data_fn_38 = sprintf('%s%s.mat', datadir, Set_38);
data_fn_39 = sprintf('%s%s.mat', datadir, Set_39);
data_fn_40 = sprintf('%s%s.mat', datadir, Set_40);
data_fn_41 = sprintf('%s%s.mat', datadir, Set_41);
data_fn_42 = sprintf('%s%s.mat', datadir, Set_42);
data_fn_43 = sprintf('%s%s.mat', datadir, Set_43);
data_fn_44 = sprintf('%s%s.mat', datadir, Set_44);
data_fn_45 = sprintf('%s%s.mat', datadir, Set_45);
data_fn_46 = sprintf('%s%s.mat', datadir, Set_46);
data_fn_47 = sprintf('%s%s.mat', datadir, Set_47);
data_fn_48 = sprintf('%s%s.mat', datadir, Set_48);

load(data_fn_1);
PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_2);
PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_3);
PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_4);
PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_5);
PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_6);
PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_7);
PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_8);
PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_9);
PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_10);
PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_11);
PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_12);
PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_13);
PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_14);
PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_15);
PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_16);
PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_17);
PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_18);
PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_19);
PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_20);
PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_21);
PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_22);
PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_23);
PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_24);
PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_25);
PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_26);
PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_27);
PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_28);
PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_29);
PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_30);
PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_31);
PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_32);
PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_33);
PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_34);
PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_35);
PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_36);
PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_37);
PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_38);
PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_39);
PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_40);
PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_41);
PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_42);
PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_43);
PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_44);
PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_45);
PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_46);
PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_47);
PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_48);
PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

end




if strcmp(Frequency_Band,'beta') == 1

datadir = 'C:\Users\Corey Roach\Documents\00_DATA\PSI_Histogram_Summary';

Set_1 = '\PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC_beta';
Set_2 = '\PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC_beta';
Set_3 = '\PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC_beta';
Set_4 = '\PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC_beta';

Set_5 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC_beta';
Set_6 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC_beta';
Set_7 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC_beta';
Set_8 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC_beta';

Set_9 =  '\PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC_beta';
Set_10 = '\PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC_beta';
Set_11 = '\PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC_beta';
Set_12 = '\PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC_beta';

Set_13 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC_beta';
Set_14 = '\PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC_beta';
Set_15 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC_beta';
Set_16 = '\PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC_beta';


Set_17 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC_beta';
Set_18 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC_beta';
Set_19 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC_beta';
Set_20 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC_beta';

Set_21 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC_beta';
Set_22 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC_beta';
Set_23 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC_beta';
Set_24 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC_beta';

Set_25 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC_beta';
Set_26 = '\PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC_beta';
Set_27 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC_beta';
Set_28 = '\PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC_beta';

Set_29 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC_beta';
Set_30 = '\PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC_beta';
Set_31 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC_beta';
Set_32 = '\PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC_beta';

Set_33 = '\PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC_beta';
Set_34 = '\PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC_beta';
Set_35 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC_beta';
Set_36 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC_beta';

Set_37 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC_beta';
Set_38 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC_beta';
Set_39 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC_beta';
Set_40 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC_beta';

Set_41 = '\PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC_beta';
Set_42 = '\PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC_beta';
Set_43 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC_beta';
Set_44 = '\PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC_beta';

Set_45 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC_beta';
Set_46 = '\PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC_beta';
Set_47 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC_beta';
Set_48 = '\PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC_beta';

data_fn_1 = sprintf('%s%s.mat', datadir, Set_1);
data_fn_2 = sprintf('%s%s.mat', datadir, Set_2);
data_fn_3 = sprintf('%s%s.mat', datadir, Set_3);
data_fn_4 = sprintf('%s%s.mat', datadir, Set_4);
data_fn_5 = sprintf('%s%s.mat', datadir, Set_5);
data_fn_6 = sprintf('%s%s.mat', datadir, Set_6);
data_fn_7 = sprintf('%s%s.mat', datadir, Set_7);
data_fn_8 = sprintf('%s%s.mat', datadir, Set_8);
data_fn_9 = sprintf('%s%s.mat', datadir, Set_9);
data_fn_10 = sprintf('%s%s.mat', datadir, Set_10);
data_fn_11 = sprintf('%s%s.mat', datadir, Set_11);
data_fn_12 = sprintf('%s%s.mat', datadir, Set_12);
data_fn_13 = sprintf('%s%s.mat', datadir, Set_13);
data_fn_14 = sprintf('%s%s.mat', datadir, Set_14);
data_fn_15 = sprintf('%s%s.mat', datadir, Set_15);
data_fn_16 = sprintf('%s%s.mat', datadir, Set_16);
data_fn_17 = sprintf('%s%s.mat', datadir, Set_17);
data_fn_18 = sprintf('%s%s.mat', datadir, Set_18);
data_fn_19 = sprintf('%s%s.mat', datadir, Set_19);
data_fn_20 = sprintf('%s%s.mat', datadir, Set_20);
data_fn_21 = sprintf('%s%s.mat', datadir, Set_21);
data_fn_22 = sprintf('%s%s.mat', datadir, Set_22);
data_fn_23 = sprintf('%s%s.mat', datadir, Set_23);
data_fn_24 = sprintf('%s%s.mat', datadir, Set_24);
data_fn_25 = sprintf('%s%s.mat', datadir, Set_25);
data_fn_26 = sprintf('%s%s.mat', datadir, Set_26);
data_fn_27 = sprintf('%s%s.mat', datadir, Set_27);
data_fn_28 = sprintf('%s%s.mat', datadir, Set_28);
data_fn_29 = sprintf('%s%s.mat', datadir, Set_29);
data_fn_30 = sprintf('%s%s.mat', datadir, Set_30);
data_fn_31 = sprintf('%s%s.mat', datadir, Set_31);
data_fn_32 = sprintf('%s%s.mat', datadir, Set_32);
data_fn_33 = sprintf('%s%s.mat', datadir, Set_33);
data_fn_34 = sprintf('%s%s.mat', datadir, Set_34);
data_fn_35 = sprintf('%s%s.mat', datadir, Set_35);
data_fn_36 = sprintf('%s%s.mat', datadir, Set_36);
data_fn_37 = sprintf('%s%s.mat', datadir, Set_37);
data_fn_38 = sprintf('%s%s.mat', datadir, Set_38);
data_fn_39 = sprintf('%s%s.mat', datadir, Set_39);
data_fn_40 = sprintf('%s%s.mat', datadir, Set_40);
data_fn_41 = sprintf('%s%s.mat', datadir, Set_41);
data_fn_42 = sprintf('%s%s.mat', datadir, Set_42);
data_fn_43 = sprintf('%s%s.mat', datadir, Set_43);
data_fn_44 = sprintf('%s%s.mat', datadir, Set_44);
data_fn_45 = sprintf('%s%s.mat', datadir, Set_45);
data_fn_46 = sprintf('%s%s.mat', datadir, Set_46);
data_fn_47 = sprintf('%s%s.mat', datadir, Set_47);
data_fn_48 = sprintf('%s%s.mat', datadir, Set_48);

load(data_fn_1);
PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_2);
PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_3);
PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_4);
PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_5);
PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_6);
PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_7);
PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_8);
PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_9);
PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_10);
PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_11);
PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_12);
PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_13);
PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_14);
PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_15);
PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_16);
PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_17);
PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_18);
PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_19);
PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_20);
PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_21);
PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_22);
PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_23);
PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_24);
PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_25);
PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_26);
PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_27);
PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_28);
PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_29);
PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_30);
PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_31);
PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_32);
PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_33);
PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_34);
PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_35);
PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_36);
PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_37);
PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_38);
PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_39);
PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC = Condition_1_PFC_2_AC;
clear Condition_1_PFC_2_AC

load(data_fn_40);
PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC = Condition_1_AC_2_PFC;
clear Condition_1_AC_2_PFC

load(data_fn_41);
PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_42);
PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_43);
PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_44);
PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_45);
PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_46);
PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

load(data_fn_47);
PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC = Condition_2_PFC_2_AC;
clear Condition_2_PFC_2_AC

load(data_fn_48);
PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC = Condition_2_AC_2_PFC;
clear Condition_2_AC_2_PFC

end




%% TestTone Pretone versus Prior collapsed across, Behavior, SNR, and Congruency.

% Concatonate PriorOnly PFC to AC
PriorOnly_PFC_to_AC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;...
                           PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;  PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC;...
                           PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];
% Concatonate PriorOnly AC to PFC
PriorOnly_AC_to_PFC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;...
                           PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;  PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC;...
                           PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];

% Concatonate PretoneOnly PFC to AC
PretoneOnly_PFC_to_AC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;...
                             PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;  PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC;...
                             PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
PretoneOnly_AC_to_PFC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;...
                             PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;  PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC;...
                             PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];


%[~, bins] = hist([PriorOnly_PFC_to_AC_cat; abs(PriorOnly_AC_to_PFC_cat); PretoneOnly_PFC_to_AC_cat; abs(PretoneOnly_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(PriorOnly_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(PriorOnly_AC_to_PFC_cat), bins);
posCounts_2 = hist(PretoneOnly_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(PretoneOnly_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges; 
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges 


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.1 .1 -65 65]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.1 .1 -55 55]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

%% Correct versus Wrong collapsed across, Condition, SNR, and Congruency.

% Concatonate Correct Trials PFC to AC
Correct_PFC_to_AC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC;...
                         PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC;PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC;...
                         PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC];


% Concatonate Correct Trials AC to PFC
Correct_AC_to_PFC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC;...
                         PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC;PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC;...
                         PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC];


% Concatonate PretoneOnly PFC to AC
Wrong_PFC_to_AC_cat =   [PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC;...
                         PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;...
                         PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
Wrong_AC_to_PFC_cat = [PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC;...
                       PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;...
                       PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];


[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(Correct_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(Correct_AC_to_PFC_cat), bins);
posCounts_2 = hist(Wrong_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(Wrong_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0, 0], 'EdgeColor', 'k', 'BarWidth',    1);                     % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0, 0, 1], 'EdgeColor', 'k', 'BarWidth', 1, 'FaceAlpha', 0.5);      % Blue with black edges  


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.15 .15 -65 65]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;


figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0, 0], 'EdgeColor', 'k', 'BarWidth',    1);                     % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0, 0, 1], 'EdgeColor', 'k', 'BarWidth', 1, 'FaceAlpha', 0.5);      % Blue with black edges 

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.15 .15 -65 65]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

%% Congruent versus Incongruent, OnlyPrior Trials collapsed across Behavior, and SNR.

% Concatonate Correct Trials PFC to AC
PriorOnly_Congruent_PFC_to_AC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;...
                                     PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC];


% Concatonate Correct Trials AC to PFC
PriorOnly_Congruent_AC_to_PFC_cat  = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;...
                                      PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC];


% Concatonate PretoneOnly PFC to AC
PriorOnly_Incongruent_PFC_to_AC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;...
                             PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
PriorOnly_Incongruent_AC_to_PFC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;...
                             PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];

%[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(PriorOnly_Congruent_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(PriorOnly_Congruent_AC_to_PFC_cat), bins);
posCounts_2 = hist(PriorOnly_Incongruent_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(PriorOnly_Incongruent_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.1 .1 -55 55]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;


figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.1 .1 -55 55]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


%% Congruent versus Incongruent, OnlyPretone Trials collapsed across Behavior, and SNR.

% Concatonate Correct Trials PFC to AC
Pretone_Congruent_PFC_to_AC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;...
                           PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC];


% Concatonate Correct Trials AC to PFC
Pretone_Congruent_AC_to_PFC_cat  = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;...
                            PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC];


% Concatonate PretoneOnly PFC to AC
Pretone_Incongruent_PFC_to_AC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;...
                             PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
Pretone_Incongruent_AC_to_PFC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;...
                             PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];

[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(Pretone_Congruent_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(Pretone_Congruent_AC_to_PFC_cat), bins);
posCounts_2 = hist(Pretone_Incongruent_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(Pretone_Incongruent_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.1 .1 -55 55]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;


figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.1 .1 -55 55]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

%% LED Epoch versus TestTone Epoch, OnlyPrior

% PreCue Epoch, PriorOnly Trials,  
% testTone Epoch, PriorOnly Trials, Behavior, SNR, Congruency 

% Concatonate PriorOnly PFC to AC
TestTone_PriorOnly_PFC_to_AC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;...
                                    PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;  PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC;...
                                    PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];
% Concatonate PriorOnly AC to PFC
TestTone_PriorOnly_AC_to_PFC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;...
                                    PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;  PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC;...
                                    PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];


% Concatonate PretoneOnly PFC to AC
LED_PriorOnly_PFC_to_AC_cat = [PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC; PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC;...
                               PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC; PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
LED_PriorOnly_AC_to_PFC_cat = [PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC; PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC;...
                               PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC; PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC];

%[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(LED_PriorOnly_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(LED_PriorOnly_AC_to_PFC_cat), bins);
posCounts_2 = hist(TestTone_PriorOnly_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(TestTone_PriorOnly_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.15 .15 -65 65]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;


figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.15 .15 -65 65]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


%% LED Epoch versus TestTone Epoch, OnlyPretone Trials 

% PreCue Epoch, PriorOnly Trials,  
% testTone Epoch, PriorOnly Trials, Behavior, SNR, Congruency 

% Concatonate PretoneOnly PFC to AC
TestTone_PretoneOnly_PFC_to_AC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC; PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;...
                             PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;  PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC;...
                             PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC; PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
TestTone_PretoneOnly_AC_to_PFC_cat = [PretoneOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC; PretoneOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;...
                             PretoneOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;  PretoneOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC;...
                             PretoneOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC; PretoneOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];


% Concatonate PretoneOnly PFC to AC
LED_PretoneOnly_PFC_to_AC_cat = [PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC;...
                                 PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC; PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
LED_PretoneOnly_AC_to_PFC_cat = [PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC;...
                               PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC; PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC];

%[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(LED_PretoneOnly_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(LED_PretoneOnly_AC_to_PFC_cat), bins);
posCounts_2 = hist(TestTone_PretoneOnly_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(TestTone_PretoneOnly_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.15 .15 -65 65]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.15 .15 -65 65]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

%% Informative LED versus Non-informative LED 

% Concatonate PretoneOnly PFC to AC
LED_PriorOnly_PFC_to_AC_cat = [PriorOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC; PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC;...
                               PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC; PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
LED_PriorOnly_AC_to_PFC_cat = [PriorOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC; PriorOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC;...
                               PriorOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC; PriorOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC];

% Concatonate PretoneOnly PFC to AC
LED_PretoneOnly_PFC_to_AC_cat = [PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_PFC_2_AC; PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_PFC_2_AC;...
                                 PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_PFC_2_AC; PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
LED_PretoneOnly_AC_to_PFC_cat = [PretoneOnly_preCueOnset_AllSNR_Correct_Congruent_AC_2_PFC; PretoneOnly_preCueOnset_AllSNR_Correct_Incongruent_AC_2_PFC;...
                               PretoneOnly_preCueOnset_AllSNR_Wrong_Congruent_AC_2_PFC; PretoneOnly_preCueOnset_AllSNR_Wrong_Incongruent_AC_2_PFC];

%[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(LED_PriorOnly_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(LED_PriorOnly_AC_to_PFC_cat), bins);
posCounts_2 = hist(LED_PretoneOnly_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(LED_PretoneOnly_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.15 .15 -65 65]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.15 .15 -65 65]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges


%% Correct versus Wrong collapsed across, Condition, SNR, and Congruency.

% Concatonate Correct Trials PFC to AC
Correct_PFC_to_AC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_PFC_2_AC;PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_PFC_2_AC;...
                         PriorOnly_testToneOnset_LowSNR_Correct_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_PFC_2_AC];


% Concatonate Correct Trials AC to PFC
Correct_AC_to_PFC_cat = [PriorOnly_testToneOnset_HighSNR_Correct_Congruent_AC_2_PFC;PriorOnly_testToneOnset_HighSNR_Correct_Incongruent_AC_2_PFC;...
                         PriorOnly_testToneOnset_LowSNR_Correct_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Correct_Incongruent_AC_2_PFC];


% Concatonate PretoneOnly PFC to AC
Wrong_PFC_to_AC_cat =   [PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_PFC_2_AC;PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_PFC_2_AC;...
                         PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_PFC_2_AC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_PFC_2_AC];

% Concatonate PretoneOnly AC to PFC
Wrong_AC_to_PFC_cat = [PriorOnly_testToneOnset_HighSNR_Wrong_Congruent_AC_2_PFC;PriorOnly_testToneOnset_HighSNR_Wrong_Incongruent_AC_2_PFC;...
                       PriorOnly_testToneOnset_LowSNR_Wrong_Congruent_AC_2_PFC; PriorOnly_testToneOnset_LowSNR_Wrong_Incongruent_AC_2_PFC];


[~, bins] = hist([Correct_PFC_to_AC_cat; abs(Correct_AC_to_PFC_cat); Wrong_PFC_to_AC_cat; abs(Wrong_AC_to_PFC_cat)], 100);

% Compute histograms for each dataset
posCounts_1 = hist(Correct_PFC_to_AC_cat, bins);
negCounts_1 = hist(abs(Correct_AC_to_PFC_cat), bins);
posCounts_2 = hist(Wrong_PFC_to_AC_cat, bins);
negCounts_2 = hist(abs(Wrong_AC_to_PFC_cat), bins);

% Normalize counts to percentages
posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

% Plot histograms

figure();

hold on;
% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   posCounts_1,   'FaceColor', [1, 0, 0], 'EdgeColor', 'k', 'BarWidth',    1);                     % Red with black edges
bar(-bins,  -negCounts_1,  'FaceColor', [0, 0, 1], 'EdgeColor', 'k', 'BarWidth', 1, 'FaceAlpha', 0.5);      % Blue with black edges  


% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
% xlim([-.15 .15])
axis([-.15 .15 -50 50]) 
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;


figure()

hold on 
bar(bins,   posCounts_2,   'FaceColor', [1, 0, 0], 'EdgeColor', 'k', 'BarWidth',    1);                     % Red with black edges
bar(-bins,  -negCounts_2,  'FaceColor', [0, 0, 1], 'EdgeColor', 'k', 'BarWidth', 1, 'FaceAlpha', 0.5);      % Blue with black edges 

% Labels and Legend
xlabel('PSI Value', 'FontSize', 18);
ylabel('Proportion of Sites', 'FontSize', 18);
axis([-.15 .15 -50 50]) 

% Set legend and adjust its font size
leg = legend({'Top-down  Direction', 'Bottom-up Direction'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size

% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;

%'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;
%FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
