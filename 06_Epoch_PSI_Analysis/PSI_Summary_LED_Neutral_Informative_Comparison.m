%% Epoch Analysis
% The epoch analysis compares PSI values during the LED presentation and
% the testTone presentation 
% [OnlyPrior = informative cues and testtone; OnlyPretone = noninformative
% LED/pretones and then testTone]. There are two separate files one for
% OnlyPrior and one for OnlyPretone for the user convienence, could be
% consolidated into a single script. 

%% Step 1: Book-Keeping 

Condition         = 'OnlyPrior';                                                                                             % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];                                              % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrCassius'};                                                                                           % Options: 'MrCassius' or 'MrM'
rootdir           = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis_OnlyPrior_All_Congruency';                                    
sessions          = dir(fullfile(rootdir, '19*'));
                                                                                      % Options: 'yes' or 'no'

%% Step 2: Loop through all sessions and sort PSI values by channel-pair 

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
preCueOnset_PSI_grid   = cell(20, 20); 
testToneOnset_PSI_grid = cell(20, 20);


for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            preCue_file   = dir(fullfile(session_path, ['*', 'preCueOnset_' Frequency_Band '*PSI_data.mat']));
            testTone_file = dir(fullfile(session_path, ['*', 'testToneOnset_' Frequency_Band '*PSI_data.mat']));

            if ~isempty(preCue_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(preCue_file)

        preCue_filename   = preCue_file(k).name;
        testTone_filename = testTone_file(k).name;

        temp_preCue_dir   = fullfile(session_path, preCue_filename);
        temp_testTone_dir = fullfile(session_path, testTone_filename);

        % Load the .mat file safely
         load(temp_preCue_dir,'PSI_preCueOnset' );
         load(temp_testTone_dir,'PSI_testToneOnset');

        preCue_Matrix     = PSI_preCueOnset.psispctrm; 
        preCue_labels     = PSI_preCueOnset.labelcmb;         
        testTone_Matrix   = PSI_testToneOnset.psispctrm; 
        testTone_labels   = PSI_testToneOnset.labelcmb;       

      end    

           % Populate preCueOnset_PSI_grid with data 

           for rowIdx_1 = 1:size(preCue_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = preCue_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = preCue_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    preCueOnset_PSI_grid{pfcIdx, acIdx} = preCue_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    preCueOnset_PSI_grid{pfcIdx, acIdx} = [preCueOnset_PSI_grid{pfcIdx, acIdx}; preCue_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate testToneOnset_PSI_grid with data 

           for rowIdx_2 = 1:size(testTone_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = testTone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = testTone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    testToneOnset_PSI_grid{pfcIdx, acIdx} = testTone_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    testToneOnset_PSI_grid{pfcIdx, acIdx} = [testToneOnset_PSI_grid{pfcIdx, acIdx}; testTone_Matrix(rowIdx_2, :)];
                    end                 
                end
           end

         end
     end
end

%% Step 3: Average selected bands to output a single PSI value per channel pair 

% PreCue Epoch
if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
preCue_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = preCueOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            preCue_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end


% TestTone Epoch 

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;

% Initialize a new grid to store the averaged values
testTone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            dataposition = testToneOnset_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            columnAverage = mean(dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            testTone_averagedGrid(pfcIdx, acIdx) = columnAverage;
        end
    end
end
end



%% Step 6: Plot with Mean of Means + 95% CI [SUMMARY OF HEAT MAP]

    abs_preCue   = abs(preCue_averagedGrid);
    abs_testTone = abs(testTone_averagedGrid);
    
    % Mask out zeros (placeholder values)
    mask_nonzero_pre  = preCue_averagedGrid     ~= 0;
    mask_nonzero_test = testTone_averagedGrid   ~= 0;
    
    % Extract directional values
    PFC_2_AC_PreCue   = abs_preCue(  mask_nonzero_pre  & preCue_averagedGrid     > 0);
    AC_2_PFC_PreCue   = abs_preCue(  mask_nonzero_pre  & preCue_averagedGrid     < 0);
    PFC_2_AC_TestTone = abs_testTone(mask_nonzero_test & testTone_averagedGrid   > 0);
    AC_2_PFC_TestTone = abs_testTone(mask_nonzero_test & testTone_averagedGrid   < 0);
    
    group_data = {PFC_2_AC_PreCue, AC_2_PFC_PreCue, PFC_2_AC_TestTone, AC_2_PFC_TestTone};
    means = cellfun(@mean, group_data);
    SEMs  = cellfun(@(x) std(x)/sqrt(length(x)), group_data);
    CIs   = SEMs * 1.96;  % 95% CI

