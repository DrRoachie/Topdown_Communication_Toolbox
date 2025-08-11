addpath 'C:\Users\ARL\Documents\fieldtrip-20230613';
ft_defaults
addpath(genpath('D:\Top_Down_Coherence_Project\00_DATA'));


%% Housekeeping

DATA_DIR    = 'D:\Top_Down_Coherence_Project\00_DATA\02b_Preprocessed';
SAVE_DIR    = 'D:\Top_Down_Coherence_Project\00_DATA\05b_Signficant_Channels_epoch_fullArray'; 

session_info = {
    'MrMiyagi', '190417';
    'MrMiyagi', '190422';
    'MrMiyagi', '190425';
    'MrMiyagi', '190427';
    'MrMiyagi', '190502';
    'MrMiyagi', '190514';
    'MrMiyagi', '190525';
    'MrMiyagi', '190527';
    'MrMiyagi', '190530';
    'MrMiyagi', '190601';
    'MrMiyagi', '190604';
    'MrMiyagi', '190704';
    'MrMiyagi', '190709';
    'MrMiyagi', '190717';
    'MrMiyagi', '190719';
    'MrMiyagi', '190722';
};

% session_info = {
%     'MrCassius', '190330';
%     'MrCassius', '190404';
%     'MrCassius', '190413';
%     'MrCassius', '190416';
%     'MrM', '190417';
%     'MrCassius', '190418';
%     'MrCassius', '190419';
%     'MrCassius', '190421';
%     'MrM', '190422';
%     'MrCassius', '190423';
%     'MrM', '190425';
%     'MrM', '190427';
%     'MrM', '190502';
%     'MrM', '190514';
%     'MrCassius', '190515';
%     'MrCassius', '190517';
%     'MrM', '190525';
%     'MrM', '190527';
%     'MrM', '190530';
%     'MrCassius', '190531';
%     'MrM', '190601';
%     'MrCassius', '190603';
%     'MrM', '190604';
%     'MrCassius', '190605';
%     'MrCassius', '190703';
%     'MrM', '190704';
%     'MrM', '190709';
%     'MrCassius', '190711';
%     'MrCassius', '190713';
%     'MrM', '190717';
%     'MrM', '190719';
%     'MrCassius', '190720';
%     'MrM', '190722';
%     'MrCassius', '190723';
%     'MrCassius', '190725'
% };

Epoch = {'testToneOnset' ; 'preCueOnset'};          % options: 'preCueOnset'
behavior_list = {'Correct'};                        % options: 'Wrong' 
condition_list = {'OnlyPrior' ; 'OnlyPretone'};     % options: 'OnlyPretone'
freq_band_list = {'theta' ;  'beta'};               % options: 'alpha'; 'beta'; 'gamma'; 'highGamma'
%%

for i = 1:size(session_info,1)
    Animal = session_info{i,1};
    RecDate = session_info{i,2};
    
    for e = 1:length(Epoch)
        epochName = Epoch{e};

        for c = 1:length(condition_list)
            Condition = condition_list{c};

            for b = 1:length(behavior_list)
                Behavior = behavior_list{b};

                fName = fullfile(DATA_DIR, 'MrM', epochName, RecDate, sprintf('%s-%s_bdLFP_%s_ft.mat', Animal, RecDate, epochName));
                if ~exist(fName, 'file')
                    warning('Skipping missing file: %s', fName);
                    continue;
                end

                fprintf('\nProcessing: %s | %s | %s\n', Animal, RecDate, epochName);
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

                        % For OnlyPrior & OnlyPretone Trials, get the time-freq series for both the data and the scrambled data

                        % this chunk takes about 5-10 minutes

                        if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
                        [ tfreq_c, tfreq_scrambled_c] = TimeFreq_Estimation_fullArray(epochName, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
                        end

                        if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
                        [ tfreq_w, tfreq_scrambled_w] = TimeFreq_Estimation_fullArray(epochName, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled);
                        end


                        % For OnlyPrior Trials, statistically compare OnlyPrior data with phase-scrambled match 

                        if strcmp(Behavior, 'Correct') == 1 || strcmp(Behavior, 'Both') == 1  
                            ScrambleTest_fullArray(freq_band_list, Condition, Behavior, RecDate, epochName, Animal, data, tfreq_c, tfreq_scrambled_c, SAVE_DIR);
                        end

                        if strcmp(Behavior, 'Wrong') == 1 || strcmp(Behavior, 'Both') == 1  
                            ScrambleTest_fullArray(freq_band_list, Condition, Behavior, RecDate, epochName, Animal, data, tfreq_c, tfreq_scrambled_c, SAVE_DIR);
                        end

                    % use this clear instead
                    clear tfreq_c tfreq_scrambled_c tfreq_w tfreq_scrambled_w data_phase_scrambled processed_data processed_data_array

            end
        end
    end
end


