%% SpecBalanced Pooled Z-Scored PSI & Coherence Visualization

% This code visualizes connectivity between PFC and AC based on SpecBalanced channel pairs
% using both Phase Slope Index (PSI) and Coherence metrics.
% Visualize:
%    - Z-scored PSI and Coherence as heatmaps (Prior, Pretone)
%    - Spectra plots across frequency with 95% CI bands for Prior and Pretone (PSI and Coherence)
%    - Histograms of PSI values (signed → no abs) to visualize PFC→AC vs AC→PFC directionality
%    - Histograms of Coherence values (prior vs pretone) 


% 
Condition         = 'OnlyPrior';                                                 % Options: 'OnlyPrior' or 'OnlyPretone' or 'Both'
Frequency_Band    = 'theta';
Behavior          = 'correct';
Congruency        = 'congruent';
SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000];  % Options: any combination of  -11/6, -5/3, -1.5000, -1.2500, 0, 1.2500, 1.5000, 5/3, 11/6; 
animals           = {'MrM'};                                                     % Options: 'MrCassius' and/or 'MrM'
rootdir           = 'C:\Users\auditory research la\OneDrive\Desktop\Connectivity_Results';                                 
sessions          = dir(fullfile(rootdir, '19*'));

% Define the PFC and AC channel labels and generate temporary grids

channels               = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);
Mean_OnlyPrior_PSI_grid   = cell(20, 20); 
Mean_OnlyPretone_PSI_grid = cell(20, 20);
Z_OnlyPrior_PSI_grid   = cell(20, 20); 
Z_OnlyPretone_PSI_grid = cell(20, 20);
ZDiff_PSI_grid     = cell(20, 20); 

% Adding grids for coherence 
Mean_OnlyPrior_Coh_grid   = cell(20, 20); 
Mean_OnlyPretone_Coh_grid = cell(20, 20);
Z_OnlyPrior_Coh_grid   = cell(20, 20); 
Z_OnlyPretone_Coh_grid = cell(20, 20);
ZDiff_Coh_grid     = cell(20, 20); 



%% Loop through all sessions and sort PSI values by channel-pair 

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if exist(session_path, 'dir') % Check if the directory exists
            PSI_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));
            PSI_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_PSI.mat']));

            Coh_OnlyPrior_file   = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPrior_',  Behavior,'_', Congruency, '_', Frequency_Band,'_Coh.mat']));
            Coh_OnlyPretone_file = dir(fullfile(session_path, [Animal, '_', RecDate, '_', 'OnlyPretone_', Behavior,'_', Congruency, '_', Frequency_Band,'_Coh.mat']));

            if ~isempty(PSI_OnlyPrior_file)
                fprintf('Processing session %s:\n', RecDate);
            else
                continue;
            end

      for k = 1:length(PSI_OnlyPrior_file)

        PSI_OnlyPrior_filename   = PSI_OnlyPrior_file(k).name;
        PSI_OnlyPretone_filename = PSI_OnlyPretone_file(k).name;

        temp_OnlyPrior_dir   = fullfile(session_path, PSI_OnlyPrior_filename);
        temp_OnlyPretone_dir = fullfile(session_path, PSI_OnlyPretone_filename);

        % Load the .mat file safely
        load(temp_OnlyPrior_dir,'PSI_OnlyPrior' );
        load(temp_OnlyPretone_dir,'PSI_OnlyPretone');

        load(fullfile(session_path, Coh_OnlyPrior_file.name), 'Coh_OnlyPrior');
        load(fullfile(session_path, Coh_OnlyPretone_file.name), 'Coh_OnlyPretone');

        Mean_OnlyPrior_Matrix   = PSI_OnlyPrior.psispctrm; 
        Mean_OnlyPretone_Matrix = PSI_OnlyPretone.psispctrm; 
        Z_OnlyPrior_Matrix      = PSI_OnlyPrior.psispctrm; 
        Z_OnlyPrior_labels      = PSI_OnlyPrior.labelcmb;         
        Z_OnlyPretone_Matrix    = PSI_OnlyPretone.psispctrm; 
        Z_OnlyPretone_labels    = PSI_OnlyPretone.labelcmb;  

        % z-score the data --> pooled (both conditions)
        % Step 1: Pool all values from both matrices

        pooledData = [Z_OnlyPrior_Matrix(:); Z_OnlyPretone_Matrix(:)];

        % Step 2: Compute global mean and std
        globalMean = mean(pooledData);
        globalStd  = std(pooledData);


        % Step 3: Z-score both matrices
        Z_OnlyPrior_Matrix   = (Z_OnlyPrior_Matrix   - globalMean) / globalStd;
        Z_OnlyPretone_Matrix = (Z_OnlyPretone_Matrix - globalMean) / globalStd;

        % z-score the data --> within condition 
        % Z-score OnlyPrior within condition
        globalMean_Prior = mean(Z_OnlyPrior_Matrix(:));
        globalStd_Prior  = std(Z_OnlyPrior_Matrix(:));
        Z_OnlyPrior_Matrix = (Z_OnlyPrior_Matrix - globalMean_Prior) / globalStd_Prior;
        
        % Z-score OnlyPretone within condition
        globalMean_Pretone = mean(Z_OnlyPretone_Matrix(:));
        globalStd_Pretone  = std(Z_OnlyPretone_Matrix(:));
        Z_OnlyPretone_Matrix = (Z_OnlyPretone_Matrix - globalMean_Pretone) / globalStd_Pretone;


        % Match label combinations between OnlyPrior and OnlyPretone
        % Outputs:
        %    - commonLabels = the shared pairs
        %    - idxPrior = indices of those pairs in Prior labels 
        %    - idxPretone = indices of those pairs in Pretone labels 
        [commonLabels, idxPrior, idxPretone] = intersect( ...               % first argument gives all Prior pairs as strings
            strcat(Z_OnlyPrior_labels(:,1), Z_OnlyPrior_labels(:,2)), ...   % second argument gives all Pretone pairs as strings 
            strcat(Z_OnlyPretone_labels(:,1), Z_OnlyPretone_labels(:,2)) ...
        );
        
        % Subset matrices to only include matched channel pairs
        Z_OnlyPrior_Matrix = Z_OnlyPrior_Matrix(idxPrior, :);
        Z_OnlyPretone_Matrix = Z_OnlyPretone_Matrix(idxPretone, :);
        
        Z_OnlyPrior_labels = Z_OnlyPrior_labels(idxPrior, :);
        Z_OnlyPretone_labels = Z_OnlyPretone_labels(idxPretone, :);

        zDiff_Matrix         = Z_OnlyPrior_Matrix - Z_OnlyPretone_Matrix;   % since # of prior and pretone pairs are different --> error. 
        
        % === Coherence Z-scoring ===
        
        % Match label combinations between OnlyPrior and OnlyPretone for Coherence
        [commonLabels_Coh, idxPrior_Coh, idxPretone_Coh] = intersect( ...
            strcat(Coh_OnlyPrior.labelcmb(:,1), Coh_OnlyPrior.labelcmb(:,2)), ...
            strcat(Coh_OnlyPretone.labelcmb(:,1), Coh_OnlyPretone.labelcmb(:,2)) ...
        );
        
        % Subset matrices to only include matched channel pairs
        Coh_OnlyPrior_Matrix = Coh_OnlyPrior.cohspctrm(idxPrior_Coh, :);
        Coh_OnlyPretone_Matrix = Coh_OnlyPretone.cohspctrm(idxPretone_Coh, :);
        
        % Align labelcmb too for plotting later if needed
        Coh_OnlyPrior.labelcmb = Coh_OnlyPrior.labelcmb(idxPrior_Coh, :);
        Coh_OnlyPretone.labelcmb = Coh_OnlyPretone.labelcmb(idxPretone_Coh, :);

        % === Coherence Z-scoring with pooled conditions ===

        % Step 1: Pool all values from both matrices
        pooledCohData = [Coh_OnlyPrior.cohspctrm(:); Coh_OnlyPretone.cohspctrm(:)];
        
        % Step 2: Compute global mean and std
        globalMean_Coh = mean(pooledCohData);
        globalStd_Coh  = std(pooledCohData);
        
        % Step 3: Z-score both matrices
        Z_OnlyPrior_Coh_Matrix   = (Coh_OnlyPrior.cohspctrm - globalMean_Coh) / globalStd_Coh;
        Z_OnlyPretone_Coh_Matrix = (Coh_OnlyPretone.cohspctrm - globalMean_Coh) / globalStd_Coh;
        
        % compute ZDiff if needed
        % zDiff_Coh_Matrix = Z_OnlyPrior_Coh_Matrix - Z_OnlyPretone_Coh_Matrix;

        
        % % Z-score OnlyPrior within condition
        % globalMean_Coh_Prior = mean(Coh_OnlyPrior_Matrix(:));
        % globalStd_Coh_Prior  = std(Coh_OnlyPrior_Matrix(:));
        % Z_OnlyPrior_Coh_Matrix = (Coh_OnlyPrior_Matrix - globalMean_Coh_Prior) / globalStd_Coh_Prior;
        % 
        % % Z-score OnlyPretone within condition
        % globalMean_Coh_Pretone = mean(Coh_OnlyPretone_Matrix(:));
        % globalStd_Coh_Pretone  = std(Coh_OnlyPretone_Matrix(:));
        % Z_OnlyPretone_Coh_Matrix = (Coh_OnlyPretone_Matrix - globalMean_Coh_Pretone) / globalStd_Coh_Pretone;
        
        
        % zDiff_Coh_Matrix = Z_OnlyPrior_Coh_Matrix - Z_OnlyPretone_Coh_Matrix;


      end    

           % Populate Z_OnlyPrior_Coherence_grid with data 

           for rowIdx_1 = 1:size(Z_OnlyPrior_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Z_OnlyPrior_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx} = Mean_OnlyPrior_Matrix(rowIdx_1, :);
                    Z_OnlyPrior_PSI_grid{pfcIdx, acIdx} = Z_OnlyPrior_Matrix(rowIdx_1, :);
                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = Coh_OnlyPrior.cohspctrm(rowIdx_1, :);
                    Z_OnlyPrior_Coh_grid{pfcIdx, acIdx} = Z_OnlyPrior_Coh_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx}; Mean_OnlyPrior_Matrix(rowIdx_1, :)];
                    Z_OnlyPrior_PSI_grid{pfcIdx, acIdx} = [Z_OnlyPrior_PSI_grid{pfcIdx, acIdx}; Z_OnlyPrior_Matrix(rowIdx_1, :)];
                    Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx}; Coh_OnlyPrior.cohspctrm(rowIdx_1, :)];
                    Z_OnlyPrior_Coh_grid{pfcIdx, acIdx} = [Z_OnlyPrior_Coh_grid{pfcIdx, acIdx}; Z_OnlyPrior_Coh_Matrix(rowIdx_1, :)];
                    end

                end
           end

           % Populate OnlyPretone_PSI_grid with data 

           for rowIdx_2 = 1:size(Z_OnlyPretone_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_OnlyPretone_labels{rowIdx_2, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_OnlyPretone_labels{rowIdx_2, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once'); % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels, strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(Z_OnlyPretone_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx} = Mean_OnlyPretone_Matrix(rowIdx_2, :);
                    Z_OnlyPretone_PSI_grid{pfcIdx, acIdx} = Z_OnlyPretone_Matrix(rowIdx_2, :);
                    Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx} = Coh_OnlyPretone.cohspctrm(rowIdx_2, :);
                    Z_OnlyPretone_Coh_grid{pfcIdx, acIdx} = Z_OnlyPretone_Coh_Matrix(rowIdx_2, :);
                    else
                    % Append the new row of data to the existing data
                    Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx} = [Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx}; Mean_OnlyPretone_Matrix(rowIdx_2, :)];
                    Z_OnlyPretone_PSI_grid{pfcIdx, acIdx} = [Z_OnlyPretone_PSI_grid{pfcIdx, acIdx}; Z_OnlyPretone_Matrix(rowIdx_2, :)];
                    Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx} = [Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx}; Coh_OnlyPretone.cohspctrm(rowIdx_2, :)];
                    Z_OnlyPretone_Coh_grid{pfcIdx, acIdx} = [Z_OnlyPretone_Coh_grid{pfcIdx, acIdx}; Z_OnlyPretone_Coh_Matrix(rowIdx_2, :)];
                    end                 
                end
           end
           
           % Populate ZDiff_Coherence_PSI_grid with data 

           for rowIdx_1 = 1:size(zDiff_Matrix, 1)

            % Extract channel pair from label
            pfcLabel = Z_OnlyPrior_labels{rowIdx_1, 1};  % First column: PFC label (e.g., 'D1_PFC_ch05')
            acLabel  = Z_OnlyPrior_labels{rowIdx_1, 2};  % Second column: AC label (e.g., 'D3_AC_ch04')

            % Remove the 'D' component and isolate the channel labels
            pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');  % Extract PFC channel
            acChannel  = regexp(acLabel, 'AC_ch\d+', 'match', 'once');    % Extract AC channel

            % Find indices in the grid
            pfcIdx = find(strcmp(channels, strrep(pfcChannel, 'PFC_', '')));
            acIdx  = find(strcmp(channels,  strrep(acChannel, 'AC_', '')));

                % Assign the data row to the grid position
                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(ZDiff_PSI_grid{pfcIdx, acIdx})
                    % Initialize the cell with the first row of data
                    ZDiff_PSI_grid{pfcIdx, acIdx} = zDiff_Matrix(rowIdx_1, :);
                    else
                    % Append the new row of data to the existing data
                    ZDiff_PSI_grid{pfcIdx, acIdx} = [ZDiff_PSI_grid{pfcIdx, acIdx}; zDiff_Matrix(rowIdx_1, :)];
                    end

                end
           end
            

         end
     end
end

%% Average selected bands to output a single Coherence value per channel pair (OnlyPrior)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;
end

% Initialize a new grid to store the averaged values
Z_OnlyPrior_PSI_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions
Mean_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions
std_OnlyPrior_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

Mean_OnlyPrior_Coh_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions
std_OnlyPrior_Coh_averagedGrid     = nan(20, 20); % Optional: for std
Z_OnlyPrior_Coh_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions



% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_OnlyPrior_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition    = Z_OnlyPrior_PSI_grid{pfcIdx, acIdx};
            Mean_dataposition = Mean_OnlyPrior_PSI_grid{pfcIdx, acIdx};
            
            % Compute the average of the selected columns
            Z_columnAverage    = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage = mean(Mean_dataposition(:, columnRange), 'all', 'omitnan');
            std_columnAverage = std(Mean_dataposition(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_OnlyPrior_PSI_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage;
            Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage;
            std_OnlyPrior_averagedGrid(pfcIdx, acIdx)  = std_columnAverage;

        end
    end
end

% Loop through each position in the grid for Coherence
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_OnlyPrior_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition_Coh    = Z_OnlyPrior_Coh_grid{pfcIdx, acIdx};
            Mean_dataposition_Coh = Mean_OnlyPrior_Coh_grid{pfcIdx, acIdx};
            
            % Compute the average of the selected columns
            Z_columnAverage_Coh    = mean(Z_dataposition_Coh(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage_Coh = mean(Mean_dataposition_Coh(:, columnRange), 'all', 'omitnan');
            std_columnAverage_Coh  = std(Mean_dataposition_Coh(:, columnRange), 0, 'all', 'omitnan');

            % Store the average in the new grid
            Z_OnlyPrior_Coh_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage_Coh;
            Mean_OnlyPrior_Coh_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage_Coh;
            std_OnlyPrior_Coh_averagedGrid(pfcIdx, acIdx)  = std_columnAverage_Coh;
        end
    end
end


%% Average selected bands to output a single PSI value per channel pair (OnlyPretone)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;
end

% Initialize a new grid to store the averaged values
Z_OnlyPretone_PSI_averagedGrid    = nan(20, 20); % Use NaN to handle empty grid positions
Mean_OnlyPretone_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions
std_OnlyPretone_averagedGrid  = nan(20, 20); % Use NaN to handle empty grid positions

Mean_OnlyPretone_Coh_averagedGrid  = nan(20, 20); % Use NaN to handle empty grid positions
std_OnlyPretone_Coh_averagedGrid   = nan(20, 20); % Optional: for std
Z_OnlyPretone_Coh_averagedGrid       = nan(20, 20); % Use NaN to handle empty grid positions



% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_OnlyPretone_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition    = Z_OnlyPretone_PSI_grid{pfcIdx, acIdx};
            Mean_dataposition = Mean_OnlyPretone_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage    = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage = mean(Mean_dataposition(:, columnRange), 'all', 'omitnan');
            std_columnAverage = std(Mean_dataposition(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_OnlyPretone_PSI_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage;
            Mean_OnlyPretone_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage;
            std_OnlyPretone_averagedGrid(pfcIdx, acIdx)  = std_columnAverage;
        end
    end
end

% Loop through each position in the grid for Coherence (OnlyPretone)
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(Z_OnlyPretone_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition_Coh    = Z_OnlyPretone_Coh_grid{pfcIdx, acIdx};
            Mean_dataposition_Coh = Mean_OnlyPretone_Coh_grid{pfcIdx, acIdx};
            
            % Compute the average of the selected columns
            Z_columnAverage_Coh    = mean(Z_dataposition_Coh(:, columnRange), 'all', 'omitnan');
            Mean_columnAverage_Coh = mean(Mean_dataposition_Coh(:, columnRange), 'all', 'omitnan');
            std_columnAverage_Coh  = std(Mean_dataposition_Coh(:, columnRange), 0, 'all', 'omitnan');   

            % Store the average in the new grid
            Z_OnlyPretone_Coh_averagedGrid(pfcIdx, acIdx)    = Z_columnAverage_Coh;
            Mean_OnlyPretone_Coh_averagedGrid(pfcIdx, acIdx) = Mean_columnAverage_Coh;
            std_OnlyPretone_Coh_averagedGrid(pfcIdx, acIdx)  = std_columnAverage_Coh;

        end
    end
end


%% Average selected bands to output a single PSI value per channel pair (ZDiff)

if strcmp(Frequency_Band,'theta') == 1
% Define the range of columns to average
columnRange = 5:8;
end

if strcmp(Frequency_Band,'alpha') == 1
% Define the range of columns to average
columnRange = 9:14;
end

if strcmp(Frequency_Band,'beta') == 1
% Define the range of columns to average
columnRange = 15:30;
end

% Initialize a new grid to store the averaged values
ZDiff_averagedGrid = nan(20, 20); % Use NaN to handle empty grid positions

% Loop through each position in the grid
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(ZDiff_PSI_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition = ZDiff_PSI_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage = mean(Z_dataposition(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            ZDiff_averagedGrid(pfcIdx, acIdx) = Z_columnAverage;
        end
    end
end


% Coherence 
ZDiff_Coh_averagedGrid = nan(20, 20);
for pfcIdx = 1:20
    for acIdx = 1:20
        % Check if the grid position contains data
        if ~isempty(ZDiff_Coh_grid{pfcIdx, acIdx})
            % Extract the data matrix for the current position
            Z_dataposition_Coh = ZDiff_Coh_grid{pfcIdx, acIdx};

            % Compute the average of the selected columns
            Z_columnAverage_Coh = mean(Z_dataposition_Coh(:, columnRange), 'all', 'omitnan');

            % Store the average in the new grid
            ZDiff_Coh_averagedGrid(pfcIdx, acIdx) = Z_columnAverage_Coh;
        end
    end
end



% %% Plot the Z-scored PSI  
% 
% % Replace NaNs with zero
% Z_OnlyPrior_PSI_averagedGrid(isnan(Z_OnlyPrior_PSI_averagedGrid)) = 0;
% Z_OnlyPretone_PSI_averagedGrid(isnan(Z_OnlyPretone_PSI_averagedGrid)) = 0;
% ZDiff_averagedGrid(isnan(Z_OnlyPrior_PSI_averagedGrid)) = 0;
% 
% % Determine the global color range
% globalMin = -5.0;
% globalMax = 5.0;
% 
% % Clamp values to the global color range
% averagedGrid_1_clamped = max(globalMin, min(globalMax, Z_OnlyPrior_PSI_averagedGrid));
% averagedGrid_2_clamped = max(globalMin, min(globalMax, Z_OnlyPretone_PSI_averagedGrid));
% averagedGrid_3_clamped = max(globalMin, min(globalMax, ZDiff_averagedGrid));
% 
% % Plot the first heatmap
% figure();
% subplot(1, 3, 1);
% imagesc(averagedGrid_1_clamped);
% colormap(jet(256)); % 256 instead of default 64
% clim([globalMin, globalMax]);
% colorbar;
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title(['Zscore-OnlyPrior-',  Frequency_Band, '-', Animal]);
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% axis square;
% 
% % Plot the second heatmap
% subplot(1, 3, 2);
% imagesc(averagedGrid_2_clamped);
% colormap(jet(256)); % 256 instead of default 64
% clim([globalMin, globalMax]);
% colorbar;
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title(['Zscore-OnlyPretone-',  Frequency_Band, '-', Animal]);
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% axis square;
% 
% % Plot the third heatmap
% subplot(1, 3, 3);
% imagesc(averagedGrid_3_clamped);
% colormap(jet(256)); % 256 instead of default 64
% clim([globalMin, globalMax]);
% colorbar;
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title(['ZDiff',  Frequency_Band, '-', Animal]);
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% axis square;

% %% Plot the Z-scored Coherence   
% 
% % Replace NaNs with zero
% Z_OnlyPrior_Coh_averagedGrid(isnan(Z_OnlyPrior_Coh_averagedGrid)) = 0;
% Z_OnlyPretone_Coh_averagedGrid(isnan(Z_OnlyPretone_Coh_averagedGrid)) = 0;
% ZDiff_Coh_averagedGrid(isnan(ZDiff_Coh_averagedGrid)) = 0;
% 
% % Determine the global color range (same as PSI for consistency)
% globalMin = -5.0;
% globalMax = 5.0;
% 
% % Clamp values to the global color range
% averagedGrid_Coh_1_clamped = max(globalMin, min(globalMax, Z_OnlyPrior_Coh_averagedGrid));
% averagedGrid_Coh_2_clamped = max(globalMin, min(globalMax, Z_OnlyPretone_Coh_averagedGrid));
% averagedGrid_Coh_3_clamped = max(globalMin, min(globalMax, ZDiff_Coh_averagedGrid));
% 
% % Plot the first heatmap (OnlyPrior Coherence)
% figure();
% subplot(1, 3, 1);
% imagesc(averagedGrid_Coh_1_clamped);
% colormap(jet(256)); % 256 instead of default 64
% clim([globalMin, globalMax]);
% colorbar;
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title(['Zscore Coherence - OnlyPrior - ', Frequency_Band, ' - ', Animal]);
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% axis square;
% 
% % Plot the second heatmap (OnlyPretone Coherence)
% subplot(1, 3, 2);
% imagesc(averagedGrid_Coh_2_clamped);
% colormap(jet(256)); % 256 instead of default 64
% clim([globalMin, globalMax]);
% colorbar;
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title(['Zscore Coherence - OnlyPretone - ', Frequency_Band, ' - ', Animal]);
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% axis square;
% 
% % Plot the third heatmap (ZDiff Coherence)
% subplot(1, 3, 3);
% imagesc(averagedGrid_Coh_3_clamped);
% colormap(jet(256)); % 256 instead of default 64
% clim([globalMin, globalMax]);
% colorbar;
% xlabel('AC Channels');
% ylabel('PFC Channels');
% title(['ZDiff Coherence - ', Frequency_Band, ' - ', Animal]);
% xticks(1:20);
% yticks(1:20);
% xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
% axis square;



%% Plot RAW PSI and RAW Coherence heatmaps

% % Replace NaNs with 0 
% Mean_OnlyPrior_averagedGrid(isnan(Mean_OnlyPrior_averagedGrid)) = 0;
% Mean_OnlyPretone_averagedGrid(isnan(Mean_OnlyPretone_averagedGrid)) = 0;
% 

globalMin_PSI = -0.025;
globalMax_PSI = 0.025;

figure;

% Subplot 1: Prior
subplot(1, 2, 1);
imagesc(Mean_OnlyPrior_averagedGrid);
colormap(jet(256));
clim([globalMin_PSI, globalMax_PSI]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPrior_averagedGrid);

xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW PSI - OnlyPrior - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%  plot boxes around individual NaN cells ---
hold on;

for pfcIdx = 1:20
    for acIdx = 1:20
        if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
            % Plot a black square around this cell
            x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
            y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
            plot(x, y, 'k-', 'LineWidth', 1.5);
        end
    end
end

hold off;


% Subplot 2: Pretone
subplot(1, 2, 2);
imagesc(Mean_OnlyPretone_averagedGrid);
colormap(jet(256));
clim([globalMin_PSI, globalMax_PSI]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPrior_averagedGrid);

xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW PSI - OnlyPretone - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%  plot boxes around individual NaN cells ---
hold on;

for pfcIdx = 1:20
    for acIdx = 1:20
        if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
            % Plot a black square around this cell
            x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
            y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
            plot(x, y, 'k-', 'LineWidth', 1.5);
        end
    end
end

hold off;



% Plot RAW Coherence heatmaps: Prior + Pretone side by side

% % Replace NaNs with 0
% Mean_OnlyPrior_Coh_averagedGrid(isnan(Mean_OnlyPrior_Coh_averagedGrid)) = 0;
% Mean_OnlyPretone_Coh_averagedGrid(isnan(Mean_OnlyPretone_Coh_averagedGrid)) = 0;


globalMin_Coh = 0.005;
globalMax_Coh = max([Mean_OnlyPrior_Coh_averagedGrid(:); Mean_OnlyPretone_Coh_averagedGrid(:)]);

figure;

% Subplot 1: Prior
subplot(1, 2, 1);
imagesc(Mean_OnlyPrior_Coh_averagedGrid);
colormap(jet(256));
clim([globalMin_Coh, globalMax_Coh]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPrior_averagedGrid);

xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW Coherence - OnlyPrior - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;



%  plot boxes around individual NaN cells ---
hold on;

for pfcIdx = 1:20
    for acIdx = 1:20
        if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
            % Plot a black square around this cell
            x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
            y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
            plot(x, y, 'k-', 'LineWidth', 1.5);
        end
    end
end

hold off;


% Subplot 2: Pretone
subplot(1, 2, 2);
imagesc(Mean_OnlyPretone_Coh_averagedGrid);
colormap(jet(256));
clim([globalMin_Coh, globalMax_Coh]);
colorbar;
% Make NaNs transparent → this is the KEY step
set(gca, 'Color', 'w');  % Optional, doesn't hurt
h = gca;
h.Children.AlphaData = ~isnan(Mean_OnlyPrior_averagedGrid);

xlabel('AC Channels');
ylabel('PFC Channels');
title(['RAW Coherence - OnlyPretone - ', Frequency_Band, ' - ', Animal]);
xticks(1:20);
yticks(1:20);
xticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false));
axis square;


%  plot boxes around individual NaN cells ---
hold on;

for pfcIdx = 1:20
    for acIdx = 1:20
        if isnan(Mean_OnlyPrior_averagedGrid(pfcIdx, acIdx))
            % Plot a black square around this cell
            x = [acIdx-0.5, acIdx+0.5, acIdx+0.5, acIdx-0.5, acIdx-0.5];
            y = [pfcIdx-0.5, pfcIdx-0.5, pfcIdx+0.5, pfcIdx+0.5, pfcIdx-0.5];
            plot(x, y, 'k-', 'LineWidth', 1.5);
        end
    end
end

hold off;


%% plot the average full PSI spectra across sessions with 95% confidence interval around them.

% Initialize an empty array to store all trials
PSI_OnlyPrior_all_trials = [];

% Iterate through each cell in preCueOnset_PSI_grid
for i = 1:numel(Z_OnlyPrior_PSI_grid)
    % Concatenate the trials from the current cell
    PSI_OnlyPrior_all_trials = [PSI_OnlyPrior_all_trials; Z_OnlyPrior_PSI_grid{i}];
end

% Initialize an empty array to store all trials
PSI_OnlyPretone_all_trials = [];

% Iterate through each cell in testToneOnset_PSI_grid
for b = 1:numel(Z_OnlyPretone_PSI_grid)
    % Concatenate the trials from the current cell
    PSI_OnlyPretone_all_trials = [PSI_OnlyPretone_all_trials; Z_OnlyPretone_PSI_grid{b}];
end

% 

% Compute mean and standard error across trials
mean_data_Prior = mean(PSI_OnlyPrior_all_trials, 1);
mean_data_Pretone = mean(PSI_OnlyPretone_all_trials, 1);

% Apply absolute values to the means
mean_data_Prior = abs(mean_data_Prior);
mean_data_Pretone = abs(mean_data_Pretone);

% Calculate standard error
SE_Prior = std(PSI_OnlyPrior_all_trials, 0, 1) / sqrt(size(PSI_OnlyPrior_all_trials, 1));
SE_Pretone = std(PSI_OnlyPretone_all_trials, 0, 1) / sqrt(size(PSI_OnlyPretone_all_trials, 1));

% Calculate 95% confidence intervals
CI_Prior = tinv(0.975, size(PSI_OnlyPrior_all_trials, 1) - 1) * SE_Prior;
CI_Pretone = tinv(0.975, size(PSI_OnlyPretone_all_trials, 1) - 1) * SE_Pretone;

% % Plotting z-score 
% freq = 1:101;  % Assuming 101 frequencies
% 
% figure;
% hold on;
% 
% % Plot mean data with error bars (absolute values)
% errorbar(freq, mean_data_Prior, CI_Prior, 'b', 'LineWidth', 1.5);
% errorbar(freq, mean_data_Pretone, CI_Pretone, 'r', 'LineWidth', 1.5);
% 
% % Customize plot labels and legend
% xlabel('Frequency');
% ylabel('|PSI|');
% xlim([4 35])
% title('PSI Comparison: Prior vs Pretone');
% legend('Prior', 'Pretone', 'Location', 'best');
% 
% % Set plot aesthetics
% set(gca, 'FontSize', 12);
% grid on;
% box on;
% hold off;

% % RAW PSI Spectra
% freq = 1:101;
% figure;
% hold on;
% 
% errorbar(freq, mean(PSI_OnlyPrior_all_trials, 1), ...
%          std(PSI_OnlyPrior_all_trials, 0, 1) / sqrt(size(PSI_OnlyPrior_all_trials,1)), ...
%          'b', 'LineWidth', 1.5);
% 
% errorbar(freq, mean(PSI_OnlyPretone_all_trials, 1), ...
%          std(PSI_OnlyPretone_all_trials, 0, 1) / sqrt(size(PSI_OnlyPretone_all_trials,1)), ...
%          'r', 'LineWidth', 1.5);
% 
% xlabel('Frequency (Hz)');
% ylabel('RAW PSI');
% xlim([4 35]);
% title(['RAW PSI Spectra: Prior vs Pretone']);
% legend('Prior', 'Pretone');
% grid on;
% box on;
% hold off;
% % 
% % RAW Coherence Spectra
% figure;
% hold on;
% 
% errorbar(freq, mean(Coh_OnlyPrior_all_trials, 1), ...
%          std(Coh_OnlyPrior_all_trials, 0, 1) / sqrt(size(Coh_OnlyPrior_all_trials,1)), ...
%          'b', 'LineWidth', 1.5);
% 
% errorbar(freq, mean(Coh_OnlyPretone_all_trials, 1), ...
%          std(Coh_OnlyPretone_all_trials, 0, 1) / sqrt(size(Coh_OnlyPretone_all_trials,1)), ...
%          'r', 'LineWidth', 1.5);
% 
% xlabel('Frequency (Hz)');
% ylabel('RAW Coherence');
% xlim([4 35]);
% title(['RAW Coherence Spectra: Prior vs Pretone']);
% legend('Prior', 'Pretone');
% grid on;
% box on;
% hold off;


% Count true zeros in OnlyPrior_theta
num_zeros_Prior = sum(PSI_OnlyPrior_all_trials == 0);

% Count true zeros in OnlyPretone_theta
num_zeros_Pretone = sum(PSI_OnlyPretone_all_trials == 0);

% Plot the average full Coherence spectra across sessions with 95%
% confidence interval around them --> z-score

% % Initialize an empty array to store all trials
% Coh_OnlyPrior_all_trials = [];
% for i = 1:numel(Z_OnlyPrior_Coh_grid)
%     Coh_OnlyPrior_all_trials = [Coh_OnlyPrior_all_trials; Z_OnlyPrior_Coh_grid{i}];
% end
% 
% Coh_OnlyPretone_all_trials = [];
% for b = 1:numel(Z_OnlyPretone_Coh_grid)
%     Coh_OnlyPretone_all_trials = [Coh_OnlyPretone_all_trials; Z_OnlyPretone_Coh_grid{b}];
% end
% 
% % Compute mean and standard error across trials
% mean_data_Coh_Prior = mean(Coh_OnlyPrior_all_trials, 1);
% mean_data_Coh_Pretone = mean(Coh_OnlyPretone_all_trials, 1);
% 
% % Calculate standard error
% SE_Coh_Prior = std(Coh_OnlyPrior_all_trials, 0, 1) / sqrt(size(Coh_OnlyPrior_all_trials, 1));
% SE_Coh_Pretone = std(Coh_OnlyPretone_all_trials, 0, 1) / sqrt(size(Coh_OnlyPretone_all_trials, 1));
% 
% % Calculate 95% confidence intervals
% CI_Coh_Prior = tinv(0.975, size(Coh_OnlyPrior_all_trials, 1) - 1) * SE_Coh_Prior;
% CI_Coh_Pretone = tinv(0.975, size(Coh_OnlyPretone_all_trials, 1) - 1) * SE_Coh_Pretone;
% 
% % Plotting
% freq = 1:101;  % Assuming 101 frequencies
% 
% figure;
% hold on;
% 
% % Plot mean data with error bars
% errorbar(freq, mean_data_Coh_Prior, CI_Coh_Prior, 'b', 'LineWidth', 1.5);
% errorbar(freq, mean_data_Coh_Pretone, CI_Coh_Pretone, 'r', 'LineWidth', 1.5);
% 
% % Customize plot labels and legend
% xlabel('Frequency');
% ylabel('Coherence');
% xlim([4 35]);  % You can adjust this to your band of interest
% title('Coherence Comparison: Prior vs Pretone');
% legend('Prior', 'Pretone', 'Location', 'best');
% 
% % Set plot aesthetics
% set(gca, 'FontSize', 12);
% grid on;
% box on;
% hold off;


%% Band-wise collect all values for each channel pair, and collapse into column

% if strcmp(Frequency_Band,'theta') == 1
% 
%     % Extract columns 5 to 8 and reshape into a single column
%     OnlyPrior_theta   = PSI_OnlyPrior_all_trials(:, 5:8);
%     OnlyPrior_theta   = OnlyPrior_theta(:);
%     OnlyPretone_theta = PSI_OnlyPretone_all_trials(:, 5:8);
%     OnlyPretone_theta = OnlyPretone_theta(:);
% 
%     % Split by direction
%     OnlyPrior_PFC_2_AC   = OnlyPrior_theta(OnlyPrior_theta > 0);
%     OnlyPretone_PFC_2_AC = OnlyPretone_theta(OnlyPretone_theta > 0);
%     OnlyPrior_AC_2_PFC   = OnlyPrior_theta(OnlyPrior_theta < 0);
%     OnlyPretone_AC_2_PFC = OnlyPretone_theta(OnlyPretone_theta < 0);
% 
%     % Define bin edges from -6 to +6 (centered)
%     binEdges = linspace(-10, 10, 50);
% 
%     % Plot all histograms on the same figure
%     figure;
%     hold on;
% 
%     histogram(OnlyPrior_theta, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], ...
%               'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Blue
%     histogram(OnlyPretone_theta, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], ...
%               'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Light red
% 
%     xlabel('Phase Slope Index (PSI)');
%     ylabel('Count');
%     xlim([-10 10])
%     legend({'OnlyPrior', 'OnlyPretone'});
%     title('Z-Scored PSI Values');
%     box on;
%     set(gca, 'FontSize', 14);
%     hold off;
% 
% end

%% RAW PSI Histogram
if strcmp(Frequency_Band,'theta') == 1
    OnlyPrior_theta = PSI_OnlyPrior_all_trials(:, 5:8);
    OnlyPretone_theta = PSI_OnlyPretone_all_trials(:, 5:8);

    OnlyPrior_theta = OnlyPrior_theta(:);
    OnlyPretone_theta = OnlyPretone_theta(:);

    binEdges = linspace(-10, 10, 50);

    figure;
    hold on;

    histogram(OnlyPrior_theta, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    histogram(OnlyPretone_theta, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    xlabel('RAW PSI Value');
    ylabel('Count');
    xlim([-10 10]);
    legend({'OnlyPrior', 'OnlyPretone'});
    title('RAW PSI: Prior vs Pretone');
    box on;
    set(gca, 'FontSize', 14);
    hold off;
end

% RAW Coherence Histogram
if strcmp(Frequency_Band,'theta') == 1
    OnlyPrior_Coh_theta = Coh_OnlyPrior_all_trials(:, 5:8);
    OnlyPretone_Coh_theta = Coh_OnlyPretone_all_trials(:, 5:8);

    OnlyPrior_Coh_theta = OnlyPrior_Coh_theta(:);
    OnlyPretone_Coh_theta = OnlyPretone_Coh_theta(:);

    binEdges = linspace(0, 1, 50); % Coherence is positive, range [0,1]

    figure;
    hold on;

    histogram(OnlyPrior_Coh_theta, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    histogram(OnlyPretone_Coh_theta, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    xlabel('RAW Coherence Value');
    ylabel('Count');
    xlim([0 1]);
    legend({'OnlyPrior', 'OnlyPretone'});
    title('RAW Coherence: Prior vs Pretone');
    box on;
    set(gca, 'FontSize', 14);
    hold off;
end

%% plotting histogram 

[~, bins] = hist([OnlyPrior_PFC_2_AC; abs(OnlyPrior_AC_2_PFC); OnlyPretone_PFC_2_AC; abs(OnlyPretone_AC_2_PFC)], 100);
edges = linspace(0, 10, 10);                % 100 bins from 0 to 0.1
bins  = edges(1:end-1) + diff(edges)/2;       % Bin centers


% Compute histograms for each dataset
OnlyPrior_PFC_2_AC_counts   = hist(OnlyPrior_PFC_2_AC, bins);
OnlyPrior_AC_2_PFC_counts   = hist(abs(OnlyPrior_AC_2_PFC), bins);
OnlyPretone_PFC_2_AC_counts = hist(OnlyPretone_PFC_2_AC, bins);
OnlyPretone_AC_2_PFC_counts = hist(abs(OnlyPretone_AC_2_PFC), bins);

% Normalize counts to percentages
OnlyPrior_PFC_2_AC_counts   = OnlyPrior_PFC_2_AC_counts / sum(OnlyPrior_PFC_2_AC_counts) * 100;
OnlyPrior_AC_2_PFC_counts   = OnlyPrior_AC_2_PFC_counts / sum(OnlyPrior_AC_2_PFC_counts) * 100;
OnlyPretone_PFC_2_AC_counts = OnlyPretone_PFC_2_AC_counts / sum(OnlyPretone_PFC_2_AC_counts) * 100;
OnlyPretone_AC_2_PFC_counts = OnlyPretone_AC_2_PFC_counts / sum(OnlyPretone_AC_2_PFC_counts) * 100;

% Plot histograms

figure();

hold on;

% Condition 1: Red (Positive) & Soft Red (Negative)
bar(bins,   OnlyPrior_PFC_2_AC_counts,   'FaceColor', [1, 0.6, 0.6],   'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins, -OnlyPrior_AC_2_PFC_counts,   'FaceColor', [0.6, 0.6, 1],   'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges
bar(bins,   OnlyPretone_PFC_2_AC_counts, 'FaceColor', [0.8, 0.1, 0.1], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1);  % Soft Red with black edges;                    % Red with black edges
bar(-bins, -OnlyPretone_AC_2_PFC_counts, 'FaceColor', [0.2, 0.4, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

% Labels and Legend
xlabel('PSI Value');
ylabel('Proportion of Sites');
% ylim([-65 65])
xlim([-10 10])
%axis([-.025 .025 -65 65]) 
% Adjust legend to include all conditions
leg = legend({'LED PFC→AC', 'LED AC→PFC', 'Pretone PFC→AC', 'Pretone AC→PFC'}, 'Location', 'SouthEast');
set(leg, 'FontSize', 12); % Reduce legend font size
% Improve aesthetics
set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Keep axis font size at 18
box on;
hold off;



% 
% % Z-score coherence
% if strcmp(Frequency_Band, 'theta') == 1
% 
%     % Extract columns 5 to 8 and reshape into a single column for Coherence
%     OnlyPrior_theta_Coh   = Coh_OnlyPrior_all_trials(:, 5:8);
%     OnlyPrior_theta_Coh   = OnlyPrior_theta_Coh(:);
% 
%     OnlyPretone_theta_Coh = Coh_OnlyPretone_all_trials(:, 5:8);
%     OnlyPretone_theta_Coh = OnlyPretone_theta_Coh(:);
% 
%     % Plot histogram
%     binEdges = linspace(0, 2, 50); % Coherence values typically in [0,1] — I give up to 2 in case of z-scoring
% 
%     figure;
%     hold on;
% 
%     histogram(OnlyPrior_theta_Coh, 'BinEdges', binEdges, 'FaceColor', [0.2 0.6 1], ...
%               'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Blue
%     histogram(OnlyPretone_theta_Coh, 'BinEdges', binEdges, 'FaceColor', [1 0.4 0.4], ...
%               'FaceAlpha', 0.5, 'EdgeColor', 'k');  % Light red
% 
%     xlabel('Z-scored Coherence');
%     ylabel('Count');
%     legend({'OnlyPrior', 'OnlyPretone'});
%     title('Z-Scored Coherence Values');
%     box on;
%     set(gca, 'FontSize', 14);
%     hold off;
% 
% end

%% 


% %% Statistics
% 
% % Assuming preCue and testTone are both 20x20 cell arrays
% % and each cell contains a struct with numeric fields
% 
% 
% subtracted_PSI_grid = cell(size(OnlyPrior_PSI_grid)); % Preallocate
% 
% for i = 1:size(OnlyPrior_PSI_grid, 1)
%     for j = 1:size(OnlyPrior_PSI_grid, 2)
%         subtracted_PSI_grid{i, j} = OnlyPretone_PSI_grid{i, j} - OnlyPrior_PSI_grid{i, j};
%     end
% end
% 
% 
% subtracted_all_trials = [];
% 
% % Iterate through each cell in preCueOnset_PSI_grid
% for i = 1:numel(OnlyPrior_PSI_grid)
%     % Concatenate the trials from the current cell
%    subtracted_all_trials = [subtracted_all_trials; subtracted_PSI_grid{i}];
% end
% 
% % Extract columns 5 to 8 and reshape into a single column
% diff_theta = subtracted_all_trials(:, 5:8);
% diff_theta = diff_theta(:);
% 
% diff_PFC_2_AC = OnlyPrior_theta(OnlyPrior_theta > 0);
% diff_AC_2_PFC = OnlyPrior_theta(OnlyPrior_theta < 0);

% %% Jamie Kuminski Core/Belt Analysis
% 
% Core_OnlyPrior_averagedGrid   = OnlyPrior_averagedGrid;
% Core_OnlyPretone_averagedGrid = OnlyPretone_averagedGrid;
% save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Core_Spatial_PSI_Values.mat'), 'Core_OnlyPrior_averagedGrid', 'Core_OnlyPretone_averagedGrid');
% 
% Belt_OnlyPrior_averagedGrid   = OnlyPrior_averagedGrid;
% Belt_OnlyPretone_averagedGrid = OnlyPretone_averagedGrid;
% save(fullfile('C:\Users\Corey Roach\Desktop\Expectation_Analysis', 'Belt_Spatial_PSI_Values.mat'), 'Belt_OnlyPrior_averagedGrid', 'Belt_OnlyPretone_averagedGrid');
% 
% Subtract the session-wise means from each other and plot as a histogram.
% 
% % Subtract the conditions
% Belt_Subtracted = Belt_OnlyPrior_averagedGrid - Belt_OnlyPretone_averagedGrid;
% Core_Subtracted = Core_OnlyPrior_averagedGrid - Core_OnlyPretone_averagedGrid; 
% 
% % Reshape the matrices into column vectors
% Belt_Subtracted_Column = Belt_Subtracted(:);
% Core_Subtracted_Column = Core_Subtracted(:);
% 
% % Subtract the conditions
% Belt_Subtracted = Belt_OnlyPrior_averagedGrid - Belt_OnlyPretone_averagedGrid;
% Core_Subtracted = Core_OnlyPrior_averagedGrid - Core_OnlyPretone_averagedGrid; 
% 
% % Reshape the matrices into column vectors
% Belt_Subtracted_Column = Belt_Subtracted(:);
% Core_Subtracted_Column = Core_Subtracted(:);
% 
% % Set the same number of bins for both histograms
% numBins = 20;
% 
% % Calculate the edges for both histograms using the same bin range
% edges = linspace(min([Belt_Subtracted_Column; Core_Subtracted_Column]), max([Belt_Subtracted_Column; Core_Subtracted_Column]), numBins+1);
% 
% % Calculate the histogram counts for each condition using the same bin edges
% [counts_belt, ~] = histcounts(Belt_Subtracted_Column, edges);
% [counts_core, ~] = histcounts(Core_Subtracted_Column, edges);
% 
% % Plot the histograms with transparency
% hold on;
% 
% % Plot the first histogram for Belt_Subtracted with a specific color and transparency
% h1 = bar(edges(1:end-1), counts_belt, 'FaceColor', [0.4, 1.0, 0.4], 'EdgeColor', 'k', 'BarWidth', 1);
% set(h1, 'FaceAlpha', 0.5);  % Set transparency for Belt_Subtracted
% 
% % Plot the second histogram for Core_Subtracted with a different color and transparency
% h2 = bar(edges(1:end-1), counts_core, 'FaceColor', [1.0, 0.4, 0.4], 'EdgeColor', 'k', 'BarWidth', 1);
% set(h2, 'FaceAlpha', 0.5);  % Set transparency for Core_Subtracted
% 
% hold off;
% 
% % Customize the plot
% set(gca, 'XTick', linspace(min(edges), max(edges), numBins), 'XTickLabel', round(linspace(min(edges), max(edges), numBins), 2));  % Set x-axis labels
% ylabel('Frequency');
% xlabel('Subtracted mean PSI values');% 
% % Optionally, add a legend to differentiate the histograms
% legend({'Belt Subtracted', 'Core Subtracted'});
% 
% 
% 
% %
% 
% % Perform K-S test for normality on Core_Subtracted and Belt_Subtracted
% % Null Hypothesis: Data follows a normal distribution
% 
% % Test on Core_Subtracted (shifted to zero mean)
% [h_core, p_core] = kstest(Core_Subtracted_Column - mean(Core_Subtracted_Column));  % Center around zero
% fprintf('K-S Test for Core Subtracted: h = %d, p = %.4f\n', h_core, p_core);
% 
% % Test on Belt_Subtracted (shifted to zero mean)
% [h_belt, p_belt] = kstest(Belt_Subtracted_Column - mean(Belt_Subtracted_Column));  % Center around zero
% fprintf('K-S Test for Belt Subtracted: h = %d, p = %.4f\n', h_belt, p_belt);
% 
% % Perform K-S test to compare the distributions of Core and Belt Subtracted
% % Null Hypothesis: The two samples are from the same distribution
% 
% % Perform K-S test between the two distributions
% [h_ks, p_ks] = kstest2(Core_Subtracted_Column, Belt_Subtracted_Column);  % Compare distributions
% fprintf('K-S Test between Core and Belt Subtracted: h = %d, p = %.4f\n', h_ks, p_ks);
% 
% 
% % Skewness and Kurtosis for Core_Subtracted
% skew_core = skewness(Core_Subtracted_Column);
% kurtosis_core = kurtosis(Core_Subtracted_Column);
% 
% % Skewness and Kurtosis for Belt_Subtracted
% skew_belt = skewness(Belt_Subtracted_Column);
% kurtosis_belt = kurtosis(Belt_Subtracted_Column);
% 
% % Display the results
% fprintf('Core Subtracted - Skewness: %.4f, Kurtosis: %.4f\n', skew_core, kurtosis_core);
% fprintf('Belt Subtracted - Skewness: %.4f, Kurtosis: %.4f\n', skew_belt, kurtosis_belt);

