%% Code that takes testToneOnset epoch from rereferencing step, and cuts out the 400 ms pretone period and saves. 

% Data can now be sorted without calling upon animals and epochs
% seperately, using just a folder that contains the unsorted PP data
epoch    = 'testToneOnset';  % Options: 'preCueOnset'
datadir  = '\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\02_Preprocessed_bipolar';
savedir  = '\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\04_Epoc_Cut_bipolar';

animals = { 'MrCassius', 'MrM'};   % 'MrCassius',

% Loop through each session provided 
for a = 1:length(animals)
    animal = animals{a};
    file_animal = animal;
        if strcmp(animal, 'MrM')
            file_animal = 'MrMiyagi';
        end

    epoch_dir = fullfile(datadir, animal, epoch);
    session_dirs = dir(epoch_dir);
    session_dirs = session_dirs([session_dirs.isdir] & ~startsWith({session_dirs.name}, '.'));

    for k = 1:length(session_dirs)
        RecDate = session_dirs(k).name;
        session_path = fullfile(epoch_dir, RecDate);

        baseFileName = sprintf('%s-%s_bdLFP_%s_ft.mat', file_animal, RecDate, epoch);
        fullFileName = fullfile(session_path, baseFileName);
        if ~exist(fullFileName, 'file')
            warning('Missing file: %s', fullFileName);
            continue;
        end

        fprintf(1, 'Now reading %s\n', fullFileName);
        load(fullFileName);


    %creates the save folder with the animal and epoch name
    saveFolder = fullfile(savedir, animal, 'PretoneOnset', RecDate);
    
    % Create directory if it doesn't exist
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

   % Cuts the tesetToneOnset Epoc to 200 ms
    if strcmp(Epoch, 'testToneOnset')
        numTrial = length(data.time);
        
        for t = 1:numTrial
            % Define the new start and end times for the first trial
            % new_first_trial_start = 901; % first column of target tone 
            % new_first_trial_end = 1101;  % last column of target tone 

            new_first_trial_start = 301;  % first column of pretone period
            new_first_trial_end   = 701;  % last column of pretone period

            newtimeframe = data.time{t}(:, new_first_trial_start:new_first_trial_end);
            data.time{t} = newtimeframe;

            newtrialframe = data.trial{t}(:, new_first_trial_start:new_first_trial_end);
            data.trial{t} = newtrialframe;
        end

        % Initialize the adjusted sampleinfo with the new first trial timings
        adjusted_sampleinfo       = zeros(size(data.sampleinfo));
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
        baseFileName_Preteone = sprintf('%s-%s_bdLFP_%s_ft.mat', file_animal, RecDate, 'PretoneOnset');
        save_file_name = baseFileName_Pretone;
        save(fullfile(saveFolder, save_file_name), 'choice', 'data', 'err', 'info', 'pretone', 'pretoneLength', 'prior', 'SNR', 'stim', 'trial_id');
    end
    
  

    end
end