%% Set up alert sound when chunks are finish running

%% Set directory and get a list of all files in the folder with the desired file name pattern.

datadir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed\MrMiyagi\testToneOnset';
savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2\MrM\testTone';
sessions = dir(fullfile(datadir,'*.mat'));
addpath(genpath(datadir));



%% Define data to analyze

Animal    = 'MrM';      
Epoch     = 'testToneOnset';        % options: 'preCueOnset', 'moveOnset', & 'testToneOnset';

%% This code runs through every  

for k = 1:length(sessions)

baseFileName = sessions(k).name;
fullFileName = fullfile(sessions(k).folder, baseFileName);
fprintf(1, 'Now reading %s\n', fullFileName);
load(baseFileName)
pat2 = digitsPattern;
RecDate = extract(baseFileName, pat2);

        % Cuts the tesetToneOnset Epoc to 200 ms
        if strcmp(Epoch,'testToneOnset') == 1 
            numTrial = length(data.time);
            
            for t = 1:numTrial
            
             % Define the new start and end times for the first trial
             new_first_trial_start = 901;
             new_first_trial_end = 1101;

             newtimeframe = data.time{t}(:, new_first_trial_start:new_first_trial_end);
             data.time{t} = newtimeframe;
           
             newtrialframe = data.trial{t}(:, new_first_trial_start:new_first_trial_end);
             data.trial{t} = newtrialframe;

            end

             % Initialize the adjusted sampleinfo with the new first trial timings
             adjusted_sampleinfo = zeros(size(data.sampleinfo));
             adjusted_sampleinfo(1, :) = [new_first_trial_start, new_first_trial_end];
            
            % Calculate the duration of the first trial
            first_trial_duration = new_first_trial_end - new_first_trial_start;

            % Adjust the timings for the remaining trials
            for i = 2:size(data.sampleinfo, 1)
                % Calculate the duration of the current trial
                trial_duration =  first_trial_duration;
                
                % Set the start time to be 5 ms after the end of the previous trial
                new_start_time = adjusted_sampleinfo(i-1, 2) + 5;
                
                % Set the end time to be new_start_time plus the duration of the current trial
                new_end_time = new_start_time + trial_duration;
                
                % Update the adjusted sampleinfo
                adjusted_sampleinfo(i, :) = [new_start_time, new_end_time];
            end
            
            % Update the original data.sampleinfo with the adjusted timings
            data.sampleinfo = adjusted_sampleinfo;
          
             save_file_name = baseFileName;
                 save(fullfile(savedir,save_file_name), 'choice', 'data', 'err','info','pretone', 'pretoneLength', 'prior', 'SNR', 'stim','trial_id'); 
        
        end

        
        if strcmp(Epoch,'preCueOnset') == 1 
            numTrial = length(data.time);
            
            for t = 1:numTrial
            
             % Define the new start and end times for the first trial
             new_first_trial_start = 101;
             new_first_trial_end = 301;

             newtimeframe = data.time{t}(:, new_first_trial_start:new_first_trial_end);
             data.time{t} = newtimeframe;
           
             newtrialframe = data.trial{t}(:, new_first_trial_start:new_first_trial_end);
             data.trial{t} = newtrialframe;

            end

             % Initialize the adjusted sampleinfo with the new first trial timings
             adjusted_sampleinfo = zeros(size(data.sampleinfo));
             adjusted_sampleinfo(1, :) = [new_first_trial_start, new_first_trial_end];
            
            % Calculate the duration of the first trial
            first_trial_duration = new_first_trial_end - new_first_trial_start;

            % Adjust the timings for the remaining trials
            for i = 2:size(data.sampleinfo, 1)
                % Calculate the duration of the current trial
                trial_duration =  first_trial_duration;
                
                % Set the start time to be 5 ms after the end of the previous trial
                new_start_time = adjusted_sampleinfo(i-1, 2) + 5;
                
                % Set the end time to be new_start_time plus the duration of the current trial
                new_end_time = new_start_time + trial_duration;
                
                % Update the adjusted sampleinfo
                adjusted_sampleinfo(i, :) = [new_start_time, new_end_time];
            end
            
            % Update the original data.sampleinfo with the adjusted timings
            data.sampleinfo = adjusted_sampleinfo;
          
             save_file_name = baseFileName;
                 save(fullfile(savedir,save_file_name), 'choice', 'data', 'err','info','pretone', 'pretoneLength', 'prior', 'SNR', 'stim','trial_id'); 
        
        end
        
        % if strcmp(Epoch,'moveOnset') == 1 
        % 
        % end

end

