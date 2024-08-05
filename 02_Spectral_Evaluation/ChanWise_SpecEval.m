%% Housekeeping

DATA_DIR    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed\MrCassius\preCueOnset';
SAVE_DIR    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_04_18_Analysis'; % save file directory
addpath(genpath(DATA_DIR));
isSave = 1;

%% Scramble the OnlyPrior data 
% this chunk takes about 30 seconds - 1 minute

            % define data to analyze
            
            Animal         = 'MrCassius';        % Options: 'MrCassius', 'MrM'; 
            RecDate        = '190330';        
            Epoch          = 'preCueOnset';      % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
            Condition      = 'OnlyPrior';        % Options: prior condition
            Behavior       = 'Correct';          % Options:  'Correct', Wrong, due to memory constraints 
        
            % load fieldtrip data formatted by FormatLFP_ft_v2.m
            
            fName = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
            load(fName);
        
            % make a phase scrambled data set 
        
            data_phase_scrambled = data; 
            
            % Define your row-wise function
            function_handle = @(x) randphasespec(x); 
            
            processed_data = cell(size(data.trial));
        
                for i = 1:numel(data.trial)
                    % Get the current cell
                    current_cell = data.trial{i};
                
                    % Apply the function row-wise to the current cell's double array
                    processed_data_array = rowfun(function_handle, table(current_cell));
                
                    % Convert the table back to an array and assign it to the corresponding cell in the new array
                    processed_data{i} = table2array(processed_data_array);
                
                end
        
            data_phase_scrambled.trial = processed_data;
            
            clear processed_data
            clear processed_data_array

            %% For OnlyPrior Trials, get the time-freq series for both the data and the scrambled data

            % this chunk takes about 5-10 minutes

            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
            [ tfreq_c, tfreq_scrambled_c] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
            end

            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
            [ tfreq_w, tfreq_scrambled_w] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
            end

           
%% For OnlyPrior Trials, statistically compare OnlyPrior data with phase-scrambled match 

% make a list of frequencies 

%freq_band_list = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};

if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
% ScrambleTest(freq_band_list, Condition, Behavior, RecDate, Epoch, data, tfreq_c, tfreq_scrambled_c, SAVE_DIR)
ScrambleTest(Epoch, data, tfreq_c, tfreq_scrambled_c)
end

if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
ScrambleTest(freq_band_list, Condition, Behavior, RecDate, Epoch, data, tfreq_w, tfreq_scrambled_w, SAVE_DIR)
end

%% clear everything from the onlyPrior condition 
KeepList = {'freq_band_list', 'Animal', 'RecDate', 'Epoch', 'Behavior', 'choice', 'err', 'pretone', 'pretoneLength', 'prior', 'SNR', 'Condition', 'data', 'data_phase_scrambled', 'SAVE_DIR'};
allVars = who;
varstoClear = setdiff(allVars, KeepList);
clear(varstoClear{:});

%% Repeat for the Pretone Condition 

Condition      = 'OnlyPretone';        % Options: prior condition


%% For OnlyPretone Trials, get the time-freq series for both the data and the scrambled data

if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
[ tfreq_c, tfreq_scrambled_c] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
end

if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
[ tfreq_w, tfreq_scrambled_w] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
end

    %%   For OnlyPretone Trials, statistically compare OnlyPrior data with phase-scrambled match and save 


if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
ScrambleTest(freq_band_list, Condition, Behavior, RecDate, Epoch, data, tfreq_c, tfreq_scrambled_c, SAVE_DIR)
end

if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
ScrambleTest(freq_band_list, Condition, Behavior, RecDate, Epoch, data, tfreq_w, tfreq_scrambled_w, SAVE_DIR)
end