
% Finds highest three clusters of p-value counts from the array generated
% in Summary_PValue_Spatial_Figure.m (must run prior to this script)

% Alternatively, a save chunk could be added to
% Summary_PValue_Spatial_Figure.m and variables ChanPair_Array_AC_PFC and
% ChanPair_Array_PFC_AC can be saved and loaded in this script.

%% Define input array to be evaluated 

Animal          = 'MrCassius';                   % 'both' or 'MrCassius' or 'MrM' (for file name ONLY, should match input to Summary_PValue_Spatial_Figure.m) 
Frequency_Band  = 'theta';                 % for file name ONLY, should match input to Summary_PValue_Spatial_Figure.m
Statistic       = 'Coherence';
savedir         = 'D:\2024_09_27_Analysis\XCorr_Histogram_Data\Coherence_preCue_xcorr_histogram_data';

%% Coherence 

% CR checked through the code and inputing a 20x20 coherence array rather
% than the single directionality granger array does not break the logic of
% the cluster identification. 

if strcmp(Statistic, 'Coherence')

        Input_Array = ChanPair_Array_Coh;    % 'ChanPair_Array_PFC_AC' or 'ChanPair_Array_AC_PFC'
        
        % Generate a blurred version of the p-value count array
        
        Blur_Array = Input_Array; % make copy of input array (to repeat blur function)
        figure; heatmap(Blur_Array);
        for i = 1:3 % blur the p-value count array 3 times (determined experimentally, preserves the most defined peaks)
            Blur_Array = blurArray(Blur_Array, 2);
        end

        figure; heatmap(Blur_Array);
        title(['Coherence Smoothed Connectivity Map', ' ', Frequency_Band, ' ', Animal]);
        xlabel('AC');
        ylabel('PFC');
        
        % Find highest 3 clusters of p-values and store in array
        Blur_Array_Copy = Blur_Array; % make copy of blurred array (to be able to remove maximum and find next highest peak)
        Results = []; % array storing results of peak cluster detection (first/second/third max rank, sender channel, receiver channel, max id logical, # of p-values)
        
        for k = 1:3
            % find the next max value
            max_value = max(Blur_Array_Copy, [], 'all');
            [max_send_chan, max_rec_chan] = find(Blur_Array_Copy == max_value);
            max_send_chan = max_send_chan(1); % if there are multiple maxima, arbitrarily take the first
            max_rec_chan = max_rec_chan(1);
            % find the +/-2 range around the max value
            for ii = 0:4
                send_chan = (max_send_chan - 2) + ii; % start from (max_send_chan - 2) and increase to (max_send_chan + 2)
                if (send_chan < 1) || (send_chan > 17) % make sure channel index is within array bounds
                    continue;
                end
                for jj = 0:4
                    rec_chan = (max_rec_chan - 2) + jj; % for each row of senders, start from (max_rec_chan - 2) and increase to (max_rec_chan + 2)
                    if (rec_chan < 1) || (rec_chan > 17) % make sure channel index is within array bounds
                        continue;
                    end
                    max_chan_logical = (send_chan == max_send_chan) & (rec_chan == max_rec_chan); % logical value for identifying which pair is the max value
                    % update Results array
                    New_Row = [k, send_chan + 3, rec_chan + 3, max_chan_logical, Input_Array(send_chan, rec_chan)]; % channel numbers go from 3:22
                    Results = [Results; New_Row];
                    % remove the channel pair once it's been added to Results
                    Blur_Array_Copy(send_chan, rec_chan) = 0;
                end
            end
            
        end

        % Remove overlapping channel pairs
        
        C = unique(Results(:,2:3), 'rows'); % get all unique channel pairs
        
        for i = 1:length(C)
            pair_indices = find((C(i,1) == Results(:,2)) & C(i,2) == Results(:,3)); % get indices of every occurrence of the channel pair
            if length(pair_indices) > 1 % if there are more than one occurrence of the channel pair
                Results(pair_indices(2:end),:) = []; % delete all duplicate occurrences except for the first (highest priority goes to higher max rank)
            end
        end

        max_results = Results;

        fName   = sprintf('pvalue_max_channels_%s_%s_Coh.mat', Animal, Frequency_Band);
        
        save(fullfile(savedir,fName), 'max_results');

end


%% Granger Section 


if strcmp(Statistic, 'Granger')
       
    
        % PFC to AC Direction

        Direction       = 'PFC_AC';                % relevant when 'Statistic' = 'Granger'...Options: 'PFC_AC' or 'AC_PFC'

        % this is the chunk that you use if you ran Summary_Value_Spatial_Figure.m as Granger and want yo use the repectively strongest connection for each direction

        if strcmp(Direction, 'AC_PFC')
            Input_Array = ChanPair_Array_AC_PFC;    % 'ChanPair_Array_PFC_AC' or 'ChanPair_Array_AC_PFC'
        elseif strcmp(Direction, 'PFC_AC')
            Input_Array = ChanPair_Array_PFC_AC;
        end

        % Generate a blurred version of the p-value count array
        
        Blur_Array = Input_Array; % make copy of input array (to repeat blur function)
        figure; heatmap(Blur_Array);
        for i = 1:3 % blur the p-value count array 3 times (determined experimentally, preserves the most defined peaks)
            Blur_Array = blurArray(Blur_Array, 2);
        end

       if strcmp(Direction, 'PFC_AC')
        figure; heatmap(Blur_Array);
        title(['PFC2AC',' ', 'Granger Smoothed Connectivity Map', ' ', Frequency_Band, ' ', Animal]);
        xlabel('AC');
        ylabel('PFC');
       end

       if strcmp(Direction, 'AC_PFC')
        figure; heatmap(Blur_Array);
        title(['AC2PFC',' ', 'Granger Smoothed Connectivity Map', ' ', Frequency_Band, ' ', Animal]);
        xlabel('AC');
        ylabel('PFC');
       end

        % Find highest 3 clusters of p-values and store in array
        Blur_Array_Copy = Blur_Array; % make copy of blurred array (to be able to remove maximum and find next highest peak)
        Results = []; % array storing results of peak cluster detection (first/second/third max rank, sender channel, receiver channel, max id logical, # of p-values)
        for k = 1:3
            % find the next max value
            max_value = max(Blur_Array_Copy, [], 'all');
            [max_send_chan, max_rec_chan] = find(Blur_Array_Copy == max_value);
            max_send_chan = max_send_chan(1); % if there are multiple maxima, arbitrarily take the first
            max_rec_chan = max_rec_chan(1);
            % find the +/-2 range around the max value
            for ii = 0:4
                send_chan = (max_send_chan - 2) + ii; % start from (max_send_chan - 2) and increase to (max_send_chan + 2)
                if (send_chan < 1) || (send_chan > 17) % make sure channel index is within array bounds
                    continue;
                end
                for jj = 0:4
                    rec_chan = (max_rec_chan - 2) + jj; % for each row of senders, start from (max_rec_chan - 2) and increase to (max_rec_chan + 2)
                    if (rec_chan < 1) || (rec_chan > 17) % make sure channel index is within array bounds
                        continue;
                    end
                    max_chan_logical = (send_chan == max_send_chan) & (rec_chan == max_rec_chan); % logical value for identifying which pair is the max value
                    % update Results array
                    New_Row = [k, send_chan + 3, rec_chan + 3, max_chan_logical, Input_Array(send_chan, rec_chan)]; % channel numbers go from 3:22
                    Results = [Results; New_Row];
                    % remove the channel pair once it's been added to Results
                    Blur_Array_Copy(send_chan, rec_chan) = 0;
                end
            end
            
        end

        % Remove overlapping channel pairs
        
        C = unique(Results(:,2:3), 'rows'); % get all unique channel pairs
        
        for i = 1:length(C)
            pair_indices = find((C(i,1) == Results(:,2)) & C(i,2) == Results(:,3)); % get indices of every occurrence of the channel pair
            if length(pair_indices) > 1 % if there are more than one occurrence of the channel pair
                Results(pair_indices(2:end),:) = []; % delete all duplicate occurrences except for the first (highest priority goes to higher max rank)
            end
        end

        max_results = Results;
        fName   = sprintf('pvalue_max_channels_%s_%s_%s_%s.mat', Statistic, Animal, Frequency_Band, Direction);
        save(fullfile(savedir,fName), 'max_results');

        
        % AC to PFC Direction
         
        Direction       = 'AC_PFC';                % relevant when 'Statistic' = 'Granger'...Options: 'PFC_AC' or 'AC_PFC'

        % this is the chunk that you use if you ran Summary_Value_Spatial_Figure.m as Granger and want yo use the repectively strongest connection for each direction

        if strcmp(Direction, 'AC_PFC')
            Input_Array = ChanPair_Array_AC_PFC;    % 'ChanPair_Array_PFC_AC' or 'ChanPair_Array_AC_PFC'
        elseif strcmp(Direction, 'PFC_AC')
            Input_Array = ChanPair_Array_PFC_AC;
        end

        % Generate a blurred version of the p-value count array
        
        Blur_Array = Input_Array; % make copy of input array (to repeat blur function)
        figure; heatmap(Blur_Array);
        for i = 1:3 % blur the p-value count array 3 times (determined experimentally, preserves the most defined peaks)
            Blur_Array = blurArray(Blur_Array, 2);
        end

       if strcmp(Direction, 'PFC_AC')
        figure; heatmap(Blur_Array);
        title(['PFC2AC',' ', 'Granger Smoothed Connectivity Map', ' ', Frequency_Band, ' ', Animal]);
        xlabel('AC');
        ylabel('PFC');
       end

       if strcmp(Direction, 'AC_PFC')
        figure; heatmap(Blur_Array);
        title(['AC2PFC',' ', 'Granger Smoothed Connectivity Map', ' ', Frequency_Band, ' ', Animal]);
        xlabel('AC');
        ylabel('PFC');
       end
        % Find highest 3 clusters of p-values and store in array
        Blur_Array_Copy = Blur_Array; % make copy of blurred array (to be able to remove maximum and find next highest peak)
        Results = []; % array storing results of peak cluster detection (first/second/third max rank, sender channel, receiver channel, max id logical, # of p-values)
        for k = 1:3
            % find the next max value
            max_value = max(Blur_Array_Copy, [], 'all');
            [max_send_chan, max_rec_chan] = find(Blur_Array_Copy == max_value);
            max_send_chan = max_send_chan(1); % if there are multiple maxima, arbitrarily take the first
            max_rec_chan = max_rec_chan(1);
            % find the +/-2 range around the max value
            for ii = 0:4
                send_chan = (max_send_chan - 2) + ii; % start from (max_send_chan - 2) and increase to (max_send_chan + 2)
                if (send_chan < 1) || (send_chan > 17) % make sure channel index is within array bounds
                    continue;
                end
                for jj = 0:4
                    rec_chan = (max_rec_chan - 2) + jj; % for each row of senders, start from (max_rec_chan - 2) and increase to (max_rec_chan + 2)
                    if (rec_chan < 1) || (rec_chan > 17) % make sure channel index is within array bounds
                        continue;
                    end
                    max_chan_logical = (send_chan == max_send_chan) & (rec_chan == max_rec_chan); % logical value for identifying which pair is the max value
                    % update Results array
                    New_Row = [k, send_chan + 3, rec_chan + 3, max_chan_logical, Input_Array(send_chan, rec_chan)]; % channel numbers go from 3:22
                    Results = [Results; New_Row];
                    % remove the channel pair once it's been added to Results
                    Blur_Array_Copy(send_chan, rec_chan) = 0;
                end
            end
            
        end

        % Remove overlapping channel pairs
        
        C = unique(Results(:,2:3), 'rows'); % get all unique channel pairs
        
        for i = 1:length(C)
            pair_indices = find((C(i,1) == Results(:,2)) & C(i,2) == Results(:,3)); % get indices of every occurrence of the channel pair
            if length(pair_indices) > 1 % if there are more than one occurrence of the channel pair
                Results(pair_indices(2:end),:) = []; % delete all duplicate occurrences except for the first (highest priority goes to higher max rank)
            end
        end

        max_results = Results;
        fName   = sprintf('pvalue_max_channels_%s_%s_%s_%s.mat', Statistic, Animal, Frequency_Band, Direction);
        save(fullfile(savedir,fName), 'max_results');

end


