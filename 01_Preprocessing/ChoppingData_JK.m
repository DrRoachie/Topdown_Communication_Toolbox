%% Set up alert sound when chunks are finish running

%% Set directory and get a list of all files in the folder with the desired file name pattern.

%original directory... want to keep available
%datadir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed\MrMiyagi\testToneOnset';
%savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2\MrM\testTone';
%sessions = dir(fullfile(datadir,'*.mat'));
%addpath(genpath(datadir));

datadir = 'C:\Users\auditory research la\Desktop\Communication_Toolbox\01_Preprocessing\PP_MrCassius_JK\preCueOnset';
savedir = 'C:\Users\auditory research la\Desktop\Communication_Toolbox\01_Preprocessing\EC_MrCassius_JK\preCue';
sessions = dir(fullfile(datadir,'*.mat'));
addpath(genpath(datadir));


%% Define data to analyze

Animal    = 'MrM';      
Epoch     = 'preCueOnset';        % options: 'preCueOnset', 'moveOnset', & 'testToneOnset';

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

%% New One
% Data can now be sorted without calling upon animals and epochs
% seperately, using just a folder that contains the unsorted PP data
datadir = 'C:\Users\auditory research la\Desktop\Communication_Toolbox\01_Preprocessing\PP_DataMixed_JK';  % Non-specific folder
savedir = 'C:\Users\auditory research la\Desktop\Communication_Toolbox\01_Preprocessing\EC_Processed_Sorted';
sessions = dir(fullfile(datadir, '**', '*.mat'));  
addpath(genpath(datadir));

% Loop through each session provided 
for k = 1:length(sessions)

    baseFileName = sessions(k).name;
    fullFileName = fullfile(sessions(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName);

    % detects the animal and epoch by just looking at the base name
    if contains(baseFileName, 'MrCassius')
        Animal = 'MrCassius';
    elseif contains(baseFileName, 'MrM')
        Animal = 'MrM';
    else
        % if no epoch info found
        warning('Animal not found in file name: %s', baseFileName);
        continue;
    end

    if contains(baseFileName, 'preCueOnset')
        Epoch = 'preCueOnset';
    elseif contains(baseFileName, 'testToneOnset')
        Epoch = 'testToneOnset';
    else
        warning('Epoch not found in file name: %s', baseFileName);
        continue;
    end

    %creates the save folder with the animal and epoch name
    saveFolder = fullfile(savedir, Animal, Epoch);

    % Create directory if it doesn't exist
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

   % Cuts the tesetToneOnset Epoc to 200 ms
    if strcmp(Epoch, 'testToneOnset')
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

        % Save the processed data
        save_file_name = baseFileName;
        save(fullfile(saveFolder, save_file_name), 'choice', 'data', 'err', 'info', 'pretone', 'pretoneLength', 'prior', 'SNR', 'stim', 'trial_id');
    end
    % Cuts the preCueOnset Epoc to 200 ms
    if strcmp(Epoch, 'preCueOnset')
        numTrial = length(data.time);

        for t = 1:numTrial
            % Define the new start and end times for the first trial
            new_first_trial_start = 101;
            new_first_trial_end = 301;
disp(['Trial ' num2str(t) ' has ' num2str(size(data.time{t}, 2)) ' time points']);
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

        % Save the processed data
        save_file_name = baseFileName;
        save(fullfile(saveFolder, save_file_name), 'choice', 'data', 'err', 'info', 'pretone', 'pretoneLength', 'prior', 'SNR', 'stim', 'trial_id');
    end

    % if strcmp(Epoch, 'moveOnset')
    %   
    % end
end