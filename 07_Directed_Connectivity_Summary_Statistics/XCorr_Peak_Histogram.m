
%% create and fill arrays with xcorr peak values and time bins

% Define data

rootdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2';
rootchandir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Significant_Channels';

animals         = {'MrCassius'};
Epoch           = 'preCueOnset';
conditions      = {'OnlyPrior'; 'OnlyPretone'};
Behavior        = 'Correct';
frequency_bands = {'theta'};


for c = 1:length(conditions)

Condition           = conditions{c};

for fb = 1:length(frequency_bands)

Frequency_Band      = frequency_bands{fb};


% create time bin array and populate with peak xcorrs

Peak_Time_Bin_Array = [];   % [session date, PFC channel, AC channel, lag time bin, passing/nonoverlapping, peak xcorr value]

for j = 1:length(animals)

    Animal = animals{j};

    sessions = dir(fullfile(rootdir, Animal, extractBefore(Epoch, 'Onset'), '19*'));

    for i = 1:length(sessions)
        
        RecDate = sessions(i).name;
        chandir = fullfile(rootchandir, Animal, extractBefore(Epoch, 'Onset'), RecDate);

        disp(['Now processing ' fullfile(rootdir, Animal, extractBefore(Epoch, 'Onset'),RecDate)]);


        % Fetch the Significant Channels (output of ChanWise_SpecEval)
    
        SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
        files = dir(fullfile(chandir, SigChans_fn));
        
        if ~isempty(files)
            SigChans = load(fullfile(chandir, SigChans_fn), 'significant_channels');
            modifiedCellArray   = regexprep(SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
    
            % Separate elements that start with *PFC and *AC
            PFC_chans  =  modifiedCellArray(startsWith(modifiedCellArray, '*PFC'));
            AC_chans   =  modifiedCellArray(startsWith(modifiedCellArray, '*AC'));

            % Get all pairwise combinations
            [PFC_comb, AC_comb] = ndgrid(PFC_chans, AC_chans);
            Significant_Channels = [PFC_comb(:), AC_comb(:)];
        else
            continue;
        end

        num_channels = length(Significant_Channels(:,1));
        progress_msg = sprintf('0/%d\r', num_channels);
        fprintf('processing channel pair %s', progress_msg);

        for k = 1:num_channels

            Current_ChanPair = Significant_Channels(k,:);

            % updating progress message
            num_backspaces = length(progress_msg);
            fprintf(repmat('\b', 1, num_backspaces));

            progress_msg = sprintf('%d/%d\r', k, num_channels);
            fprintf('%s', progress_msg);

            % calculate xcorr
            [shuffled_xcorr_result, shuffled_lagsResults, xcorr_result, lagsResults] = TrialbyTrialXcorr(rootdir, Animal, RecDate, Epoch, Condition, Current_ChanPair);
    
            % find peak value and time bin
            [~, Peak_Time_Index]    = max(abs(xcorr_result.avg));
            Peak_Value              = xcorr_result.avg(Peak_Time_Index);
            Peak_Time_Bin           = lagsResults(Peak_Time_Index);
    
            % determine if confidence intervals are overlapping
            if xcorr_result.avg(Peak_Time_Index) > shuffled_xcorr_result.avg(Peak_Time_Index)
                Peak_Diff = xcorr_result.lower(Peak_Time_Index) - shuffled_xcorr_result.upper(Peak_Time_Index);
            else
                Peak_Diff = shuffled_xcorr_result.lower(Peak_Time_Index) - xcorr_result.upper(Peak_Time_Index);
            end
    
            is_nonoverlapping = Peak_Diff > 0;
    
            % add to array
            next_row = [str2double(RecDate), str2double(extractAfter(Current_ChanPair{1,1}, 'ch')), str2double(extractAfter(Current_ChanPair{1,2}, 'ch')), Peak_Time_Bin, is_nonoverlapping, Peak_Value];
            Peak_Time_Bin_Array = [Peak_Time_Bin_Array; next_row];
            
            % % plot xcorr
            % figure;
            % hold on;
            % plot(lagsResults, xcorr_result.avg, 'b', 'LineWidth', 1.5);
            % fill([lagsResults, fliplr(lagsResults)], [xcorr_result.upper, fliplr(xcorr_result.lower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % 
            % plot(lagsResults, shuffled_xcorr_result.avg, 'r', 'LineWidth', 0.5);
            % fill([lagsResults, fliplr(lagsResults)], [shuffled_xcorr_result.upper, fliplr(shuffled_xcorr_result.lower)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % 
            % close all;

        end
    end
end

% fill time bin array

Passing_Session_Indices = find(Peak_Time_Bin_Array(:,5));
Passing_Time_Bin_Array = Peak_Time_Bin_Array(Passing_Session_Indices,:);

if strcmp(Condition, 'OnlyPrior')
    OnlyPrior_Time_Bins.(Frequency_Band).alltimebins        = Peak_Time_Bin_Array;
    OnlyPrior_Time_Bins.(Frequency_Band).passingtimebins    = Passing_Time_Bin_Array;
elseif strcmp(Condition, 'OnlyPretone')
    OnlyPretone_Time_Bins.(Frequency_Band).alltimebins      = Peak_Time_Bin_Array;
    OnlyPretone_Time_Bins.(Frequency_Band).passingtimebins    = Passing_Time_Bin_Array;
end

end
end

%%

savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_24_Analysis\XCorr_Histogram_Data\preCue_xcorr_histogram_data';
file_name = 'xcorr_histogram_data.mat';
save(fullfile(savedir, file_name), 'OnlyPrior_Time_Bins', 'OnlyPretone_Time_Bins');

%% load histogram data arrays
% OnlyPrior_Time_Bins and OnlyPretone_Time_Bins saved in D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data

path_to_data    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\XCorr_Histogram_Data';
file_name       = 'xcorr_histogram_data.mat';
load(fullfile(path_to_data, file_name));

clear path_to_data;
clear file_name;

%% plot OnlyPrior and OnlyPretone xcorr peak histogram (with or without separation by sign)

freq_band       = 'theta';
filter_for_sign = true;
edge_size       = 100;

if filter_for_sign

    % positive peak
    % only prior
    bins_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    positive_peak_indices = find(bins_to_plot(:,6) > 0);   % find bins with positive peaks
    bins_to_plot = bins_to_plot(positive_peak_indices,4);

    figure; hold on
    edges = linspace(-200, 200, edge_size);
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Time');
    ylabel('Bin Count');
    title(['Positive XCorr Peaks in ' freq_band]);
    legend show

    % only pretone
    bins_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    positive_peak_indices = find(bins_to_plot(:,6) > 0);   % find bins with positive peaks
    bins_to_plot = bins_to_plot(positive_peak_indices,4);       

    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');

    % negative peak
    % only prior
    bins_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    negative_peak_indices = find(bins_to_plot(:,6) < 0);   % find bins with negative peaks
    bins_to_plot = bins_to_plot(negative_peak_indices,4);

    figure; hold on
    edges = linspace(-200, 200, edge_size);
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Time');
    ylabel('Bin Count');
    title(['Negative XCorr Peaks in ' freq_band]);
    legend show

    % only pretone
    bins_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    negative_peak_indices = find(bins_to_plot(:,6) < 0);   % find bins with positive peaks
    bins_to_plot = bins_to_plot(negative_peak_indices,4);

    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');

else

    bins_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins(:,4);

    figure; hold on
    edges = linspace(-200, 200, edge_size);
    histogram(bins_to_plot, 'BinEdges', edges);
    xlim([-200 200]);
    xlabel('Time');
    ylabel('Bin Count');
    title(['OnlyPretone XCorr Peaks in ' freq_band]);

end

%% plot channel pairs by layer group

superficial = 1:5;
uppermid    = 6:10;
lowermid    = 11:15;
deep        = 16:20;

freq_band_list  = {'theta', 'alpha', 'beta'};
pfc_layer_name  = 'Deep';        % 'Superficial' 'Lowmid' 'Uppermid' or 'Deep'
ac_layer_name   = 'Uppermid';         % 'Superficial' 'Lowmid' 'Uppermid' or 'Deep'
edges           = -200:5:200;

if strcmp(pfc_layer_name, 'Superficial') && strcmp(ac_layer_name, 'Superficial')
    pfc_layer_group  = superficial;
    ac_layer_group   = superficial;
elseif strcmp(pfc_layer_name, 'Lowermid') && strcmp(ac_layer_name, 'Uppermid')
    pfc_layer_group  = lowermid;
    ac_layer_group   = uppermid;
elseif strcmp(pfc_layer_name, 'Deep') && strcmp(ac_layer_name, 'Deep')
    pfc_layer_group  = deep;
    ac_layer_group   = deep;
elseif strcmp(pfc_layer_name, 'Deep') && strcmp(ac_layer_name, 'Uppermid')
    pfc_layer_group  = deep;
    ac_layer_group   = uppermid;
else
    error([pfc_layer_name ' to ' ac_layer_name ' layer group combination not supported.'])
end

figure; 

% positive peaks
for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    subplot(2,length(freq_band_list),i); hold on
    
    % only prior
    array_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with positive peaks and within specified layer group
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Lag (ms)');
    ylabel('Bin Count');
    title(['Positive PFC_{' pfc_layer_name '}-AC_{' ac_layer_name '} XCorr Peaks in ' freq_band]);
    legend show
    
    % only pretone
    array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with positive peaks and within specified layer group
    bins_to_plot = array_to_plot(:,4);     
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
end

% negative peaks
for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    subplot(2,length(freq_band_list),i+length(freq_band_list)); hold on
    
    % only prior
    array_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) < 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with negative peaks and within specified layer group
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Lag (ms)');
    ylabel('Bin Count');
    title(['Negative PFC_{' pfc_layer_name '}-AC_{' ac_layer_name '} XCorr Peaks in ' freq_band]);
    legend show
    
    % only pretone
    array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) < 0) & (ismember(array_to_plot(:,2), pfc_layer_group)) & (ismember(array_to_plot(:,3), ac_layer_group)),:);   % find bins with negative peaks and within specified layer group
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
end

%% plot channel pairs by animal

MrCassius   = [190330, 190404, 1907413, 190414, 190416, 190418, 190419, 190421, 190423, 190429, 190515, 190517, 190531, 190603, 190605, 190703, 190711, 190713, 190718, 190720, 190723, 190725];
MrM         = [190417, 190422, 190425, 190427, 190502, 190514, 190516, 190525, 190527, 190530, 190604, 190704, 190709, 190717, 190719, 190722];

freq_band_list  = {'theta', 'alpha', 'beta'};
animal_name     = 'MrCassius';        % 'MrCassius' or 'MrM'
edges           = -200:5:200;

if strcmp(animal_name, 'MrCassius')
    animal_recdates  = MrCassius;
elseif strcmp(animal_name, 'MrM')
    animal_recdates  = MrM;
end

figure;

% positive peaks
for i = 1:length(freq_band_list)

    freq_band = freq_band_list{i};

    subplot(2,length(freq_band_list),i); hold on
    
    % only prior
    array_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,1), animal_recdates)),:);   % find bins with positive peaks and within specified sessions
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Lag (ms)');
    ylabel('Bin Count');
    title(['Positive ' animal_name ' XCorr Peaks in ' freq_band]);
    legend show
    
    % only pretone
    array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) > 0) & (ismember(array_to_plot(:,1), animal_recdates)),:);   % find bins with positive peaks and within specified sessions
    bins_to_plot = array_to_plot(:,4);     
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
end

% negative peaks
for i = 1:length(freq_band_list)
    
    freq_band = freq_band_list{i};

    subplot(2,length(freq_band_list),i+length(freq_band_list)); hold on
    
    % only prior
    array_to_plot = OnlyPrior_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) < 0) & (ismember(array_to_plot(:,1), animal_recdates)),:);   % find bins with negative peaks and within specified sessions
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'DisplayName', 'PriorOnly');
    xlim([-200 200]);
    xlabel('Lag (ms)');
    ylabel('Bin Count');
    title(['Negative ' animal_name ' XCorr Peaks in ' freq_band]);
    legend show
    
    % only pretone
    array_to_plot = OnlyPretone_Time_Bins.(freq_band).passingtimebins;
    array_to_plot = array_to_plot((array_to_plot(:,6) < 0) & (ismember(array_to_plot(:,1), animal_recdates)),:);   % find bins with negative peaks and within specified sessions
    bins_to_plot = array_to_plot(:,4);
    
    histogram(bins_to_plot, 'BinEdges', edges, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'DisplayName', 'PretoneOnly');
end

