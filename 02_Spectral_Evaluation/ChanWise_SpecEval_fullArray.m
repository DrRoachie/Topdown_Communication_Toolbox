%% Housekeeping

DATA_DIR    = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\02b_Preprocessed';
SAVE_DIR    = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\05b_Signficant_Channels_epoch_fullArray'; 

%% Scramble the OnlyPrior data 
% this chunk takes about 30 seconds - 1 minute

Epoch     = 'testToneOnset';
Behavior  = 'Correct';

animal_list     = {'MrCassius', 'MrM'};
condition_list  = {'OnlyPrior', 'OnlyPretone'};

for a = 1:length(animal_list)
    Animal = animal_list{a};

    for c = 1:length(condition_list)
        Condition = condition_list{c};

        % --- get session list ---
        session_dir = fullfile(DATA_DIR, Animal, Epoch);
        sessions = dir(session_dir);
        sessions = sessions([sessions.isdir] & ~startsWith({sessions.name}, '.'));
        recdates = {sessions.name};

        for s = 1:length(recdates)
            RecDate = recdates{s};

        % Build expected .mat file name
        fName = fullfile(session_dir, RecDate, sprintf('%s-%s_bdLFP_%s_ft.mat', Animal, RecDate, Epoch));
        if ~exist(fName, 'file')
            warning('Skipping missing file: %s', fName);
            continue;
        end

        fprintf('\nProcessing: %s | %s | %s\n', Animal, RecDate, Epoch);
        load(fName);  % loads variable 'data'

        % make a phase scrambled data set 
        
            data_phase_scrambled = data; 
            
            % Define your row-wise function
            function_handle = @(x) randphasespec_fullArray(x); 
            
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

            %% For OnlyPrior & OnlyPretone Trials, get the time-freq series for both the data and the scrambled data

            % this chunk takes about 5-10 minutes

            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
            [ tfreq_c, tfreq_scrambled_c] = TimeFreq_Estimation_fullArray(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
            end

            if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
            [ tfreq_w, tfreq_scrambled_w] = TimeFreq_Estimation_fullArray(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
            end

           
            %% For OnlyPrior Trials, statistically compare OnlyPrior data with phase-scrambled match 
            
            % make a list of frequencies 
            
            %freq_band_list = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};
            
            if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
            significant_channels = ScrambleTest_fullArray(Condition, Behavior, RecDate, Epoch, data, tfreq_c, tfreq_scrambled_c);
        
            % Save results
            Frequency_Band = 'theta';  
            filename = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
            save_folder = fullfile(SAVE_DIR, Animal, Epoch, RecDate);
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            save(fullfile(save_folder, filename), 'significant_channels');
        end
        
        if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
            significant_channels = ScrambleTest_fullArray(Condition, Behavior, RecDate, Epoch, data, tfreq_w, tfreq_scrambled_w);
        
            % Save results
            Frequency_Band = 'theta';  
            filename = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
            save_folder = fullfile(SAVE_DIR, Animal, Epoch, RecDate);
            if ~exist(save_folder, 'dir')
                mkdir(save_folder);
            end
            save(fullfile(save_folder, filename), 'significant_channels');
        end

        % clear everything from the onlyPrior condition 
        KeepList = {'freq_band_list', 'Animal', 'RecDate', 'Epoch', 'Behavior', 'choice', 'err', 'pretone', 'pretoneLength', 'prior', 'SNR', 'Condition', 'data', 'data_phase_scrambled', 'SAVE_DIR'};
        allVars = who;
        varstoClear = setdiff(allVars, KeepList);
        clear(varstoClear{:});
        end
    end
end


% %% For OnlyPretone Trials, get the time-freq series for both the data and the scrambled data
% 
% if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
% [ tfreq_c, tfreq_scrambled_c] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
% end
% 
% if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
% [ tfreq_w, tfreq_scrambled_w] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
% end
% 
%     %%   For OnlyPretone Trials, statistically compare OnlyPrior data with phase-scrambled match and save 
% 
% 
% if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
% ScrambleTest(freq_band_list, Condition, Behavior, RecDate, Epoch, data, tfreq_c, tfreq_scrambled_c, SAVE_DIR)
% end
% 
% if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
% ScrambleTest(freq_band_list, Condition, Behavior, RecDate, Epoch, data, tfreq_w, tfreq_scrambled_w, SAVE_DIR)
% end
% 
% 
% 
% 
% 
% significant_channels = ScrambleTest_fullArray(Condition, Behavior, RecDate, Epoch, data, tfreq_c, tfreq_scrambled_c);
% 
% % Save results
% Frequency_Band = 'theta';  % or whichever band you're processing
% filename = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
% save_folder = fullfile(SAVE_DIR, Animal, Epoch, RecDate);
% if ~exist(save_folder, 'dir')
%     mkdir(save_folder);
% end
% save(fullfile(save_folder, filename), 'significant_channels');
% 
% 
