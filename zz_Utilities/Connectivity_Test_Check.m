
%% Establish the subset of data that you want to check. 
Animal          = 'MrM';                  % Options: 'MrCassius', 'MrM'; 
RecDate         = '190525';        
Epoch           = 'testToneOnset';                % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
Behavior        = 'Correct';                    % Options:  'Correct', Wrong, due to memory constraint
freq_band_list   = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};

for fb = 1:length(freq_band_list)

    Frequency_Band = freq_band_list{fb};

%% Set directory and get a list of all files in the folder with the desired file name pattern.

datadir = fullfile('D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2', Animal, extractBefore(Epoch, 'Onset'), RecDate); % The cut epoc data 
chandir = fullfile('D:\00_Significant_chans', Animal, extractBefore(Epoch, 'Onset'), RecDate); % the sig channels that are output
specdir = fullfile('D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed', Animal, 'testToneOnset'); % this is where the uncut preprocessed data is stored; you need the full time series for the spectogram
savedir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis\'; % path to RecDate folders (save structure: RecDate >> Animal >> Epoch >> all figures/files)
% savedir = 'E:\2024_07_24_Analysis';

sessions = dir(fullfile(datadir,'*.mat')); % only ever one session at a time
addpath(genpath(datadir));
addpath(genpath(chandir));

% make save folders and directory
savedir = fullfile(savedir, RecDate, Animal, extractBefore(Epoch, 'Onset'));
if ~exist(savedir, 'dir')  % make folders if they don't exist already
    error(['No save directory for ' RecDate '.']);
end


          % Fetch the OnlyPrior Channels (output of ChanWise_SpecEval)
    
            Condition = 'OnlyPrior';
            
            PriorOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
                        files = dir(fullfile(chandir, PriorOnly_SigChans_fn));
                            
                        if ~isempty(files)

                        PriorOnly_SigChans = load(fullfile(chandir, PriorOnly_SigChans_fn), 'significant_channels');
                        PriorOnly_modifiedCellArray   = regexprep(PriorOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

                        % Separate elements that start with *PFC and *AC
                        PriorOnly_PFC_chans  =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*PFC'));
                        PriorOnly_AC_chans   =  PriorOnly_modifiedCellArray(startsWith(PriorOnly_modifiedCellArray, '*AC'));
                        
                        end

             
                        if isempty(files)
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_coh.txt', Animal, RecDate, Epoch, Frequency_Band, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_zthresh.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band, Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_montecarlo.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band, Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, RecDate, Frequency_Band,Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose('all');
                        else
           

          % Fetch the OnlyPretone Channels 

            Condition = 'OnlyPretone';
            
            PretoneOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
                        files = dir(fullfile(chandir, PretoneOnly_SigChans_fn));
                            
                        if ~isempty(files)
                        PretoneOnly_SigChans = load(fullfile(chandir, PretoneOnly_SigChans_fn), 'significant_channels');
                        PretoneOnly_modifiedCellArray   = regexprep(PretoneOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');

                        % Separate elements that start with *PFC and *AC
                        PretoneOnly_PFC_chans  =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*PFC'));
                        PretoneOnly_AC_chans   =  PretoneOnly_modifiedCellArray(startsWith(PretoneOnly_modifiedCellArray, '*AC'));
                        end

                        if isempty(files)
                                        Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_coh.txt', Animal, RecDate, Epoch, Frequency_Band, Behavior));
                                        Fig2_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_zthresh.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
                                        Fig3_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_montecarlo.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
                                        save_file_name = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_%s_%s_data.txt', Animal, RecDate, Epoch, Frequency_Band,  Behavior));
                                        fid1 = fopen(Fig1_name, 'w');
                                        fid2 = fopen(Fig2_name, 'w');
                                        fid3 = fopen(Fig3_name, 'w');
                                        fid4 = fopen(save_file_name, 'w');
                                        fclose('all');
                        else

            % Get the Channel Pairs that are shared by both the OnlyPrior and OnlyPretone Condition 

            % Find common elements
            common_PFC_chans = intersect(PriorOnly_PFC_chans, PretoneOnly_PFC_chans);
            common_AC_chans = intersect(PriorOnly_AC_chans, PretoneOnly_AC_chans);

            % Get all pairwise combinations
            [PFC_comb, AC_comb] = ndgrid(common_PFC_chans, common_AC_chans);
            
            % Reshape into Nx2 cell array
            Shared_ChannelPairs = [PFC_comb(:), AC_comb(:)];

                        end
                        end
clear_list = {'AC_comb', 'common_AC_chans', 'common_PFC_chans', 'files', 'PFC_comb', 'PriorOnly_PFC_chans', 'PretoneOnly_PFC_chans', ...
                          'PriorOnly_AC_chans', 'PretoneOnly_AC_chans', 'PriorOnly_modifiedCellArray', 'PretoneOnly_modifiedCellArray', ...
                          'PriorOnly_SigChans_fn', 'PretoneOnly_SigChans_fn', 'PriorOnly_SigChans', 'PretoneOnly_SigChans'};
clear(clear_list{:});
clear clear_list;


%% Check number of files for each channel pair (should be 15)

ChanPairs_Missing = [];
ChanPairs_No_Maris = [];

all_files = dir(savedir);
file_names = {all_files.name};

for i = 1:length(Shared_ChannelPairs(:,1))

    PFC_chan = extractAfter(Shared_ChannelPairs{i,1}, '*');
    AC_chan  = extractAfter(Shared_ChannelPairs{i,2}, '*');

    matches = contains(file_names, PFC_chan) & contains(file_names, AC_chan) & contains(file_names, Frequency_Band);
    chanpair_files = file_names(matches);

    if length(chanpair_files) ~= 15
        if length(chanpair_files) == 9
            ChanPairs_No_Maris = [ChanPairs_No_Maris; Shared_ChannelPairs(i,:)];
        else
            ChanPairs_Missing = [ChanPairs_Missing; Shared_ChannelPairs(i,:)];
        end
    end
end

if isempty(Shared_ChannelPairs)
    disp(['No shared channel pairs in ' Frequency_Band '.']);
elseif isempty(ChanPairs_Missing) && isempty(ChanPairs_No_Maris)
    disp(['All channel pairs processed in ' Frequency_Band '!']);
elseif length(Shared_ChannelPairs) == length(ChanPairs_Missing)
    disp([Frequency_Band ' not run yet.'])
else
    disp([num2str(length(ChanPairs_Missing)) '/' num2str(length(Shared_ChannelPairs)) ' channel pairs missing and ' ...
        num2str(length(ChanPairs_No_Maris)) '/' num2str(length(Shared_ChannelPairs)) ' channel pairs run without Maris in ' Frequency_Band '.']);
end

end

