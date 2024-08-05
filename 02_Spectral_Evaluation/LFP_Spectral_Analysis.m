
%% Set directory and get a list of all files in the folder with the desired file name pattern.

datadir = 'C:\Users\Corey Roach\Documents\00_DATA\04_Epoc_Cut\MrCassius\LED\190330';
chandir = 'C:\Users\Corey Roach\Documents\00_DATA\2024_04_26_Analysis';
savedir = 'C:\Users\Corey Roach\Documents\00_DATA\2024_04_29_Analysis\MrCassius\LED';
sessions = dir(fullfile(datadir,'*.mat'));
addpath(genpath(datadir));
addpath(genpath(chandir));
isSave = 0; 

%% Establish the subset of data that you want to process. 
Animal         = 'MrCassius';                  % Options: 'MrCassius', 'MrM'; 
RecDate        = '190330';        
Epoch          = 'preCueOnset';                % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
Behavior       = 'Correct';                    % Options:  'Correct', Wrong, due to memory constraints 

%%

freq_band_list = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};
%freq_band_list  = {'highGamma'}; % all the frequency bands we are going to independently test
threshold_value = 0.5;       % the threshold used for the monte carlo test between PriorOnly and PretoneOnly Conditions

if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1

    for fb = 1:length(freq_band_list)

    Frequency_Band = freq_band_list{fb};
  
   % Get the channel groupings with only spectrally significant channels for both PretoneOnly and PriorOnly trials
    PriorOnly_SigChans_c_fn   = sprintf('PriorOnly_GroupedSigChans_%s_%s_%s_%s.mat', RecDate, Epoch, Behavior, Frequency_Band);
    
    files = dir(fullfile(chandir, PriorOnly_SigChans_c_fn));
    
        if ~isempty(files)
           PriorOnly_SigChans_c = load(fullfile(chandir, PriorOnly_SigChans_c_fn), 'PriorOnly_Chans_Filtered');         % for each frequency band, pull the spectrally relevant channels in the PriorOnly Condition 
        else
            disp('File does not exist.');
        end

                % fields = fieldnames(PriorOnly_SigChans_c.PriorOnly_Chans_Filtered);
                % uniqueValues = {};     
                % 
                % for j = 1:numel(fields)
                % currentField = PriorOnly_SigChans_c.PriorOnly_Chans_Filtered.(fields{j});
                % 
                % % Initialize arrays to store unique elements from both columns
                % allUniqueElements = {};
                % 
                %        % Iterate over each row of the cell array
                %         for j = 1:size(currentField, 1)
                %         % Extract elements from the first column and add to allUniqueElements
                %         allUniqueElements = [allUniqueElements, currentField{j, 1}];
                %         % Extract elements from the second column and add to allUniqueElements
                %         allUniqueElements = [allUniqueElements, currentField{j, 2}];
                %         end
                % 
                %        % Get unique elements from both columns
                %        uniqueElements = unique(allUniqueElements);
                %        % Store the result back into the parent structure
                %        PriorOnly_SigChans_st_c.PriorOnly_Chans_Filtered.(fields{j}) = uniqueElements;
               
                end

        PretoneOnly_SigChans_c_fn = sprintf('PretoneOnly_GroupedSigChans_%s_%s_%s_%s.mat', RecDate, Epoch, Behavior, Frequency_Band);
        files = dir(fullfile(chandir, PretoneOnly_SigChans_c_fn));
        
        if ~isempty(files)
            PretoneOnly_SigChans_c = load(fullfile(chandir, PretoneOnly_SigChans_c_fn), 'PretoneOnly_Chans_Filtered');   % for each frequency band, pull the spectrally relevant channels in the PretoneOnly Condition 
        else
            disp('File does not exist.');
        end
        
        % Carry out the Coherence Test for every layer grouping, graph, and save data
            LayerChan_List      =  {'PFClowmid_ACupmid'; 'PFCupmid_ACupmid'; 'PFCdeep_ACupmid'; 'PFCdeep_ACdeep'; 'PFCupper_ACupper'};
            %LayerGrouping_List =  {'PFClowmid_ACupmid'};
            
            for lc = 1:length(LayerGrouping_List)

                LayerGroup = LayerGrouping_List{lc};
                Condition = 'OnlyPrior';                          % options: 'OnlyPrior', 'OnlyPretone', & 'Both';
                Channels  = PriorOnly_SigChans_st_c.PriorOnly_Chans_Filtered.(LayerGroup);

                        if isempty(Channels)
                            Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                            fid1 = fopen(Fig1_name, 'w');
                            fclose(fid1);
                           
                        else

           
                Condition = 'OnlyPretone';                          % options: 'OnlyPrior', 'OnlyPretone', & 'Both';
                Channels  = PretoneOnly_SigChans_st_c.PretoneOnly_Chans_Filtered.(LayerGroup);

                        if isempty(Channels) 
                            Fig1_name      = fullfile(savedir, sprintf('NULL_%s_%s_%s_%s_coh.txt', RecDate, Frequency_Band, LayerGroup, Behavior));
                            fid1 = fopen(Fig1_name, 'w');
                            fclose(fid1);
                           
                        else

                        end

                        end
            end
    end
end


                % load fieldtrip data formatted by bipolarLFP.m
                
                fName = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
                %fName = strcat(Animal,'-',RecDate,'_bipolarLFP_',Epoch);
                load(fName);

%% set parameters required for fieldtrip functions and data selection
%Used to for selectData function

params.choice = choice;
params.err = err;
params.pretone = pretone;
params.pretoneLength = pretoneLength;
params.prior = cell2char(prior); %cell2char is an auxillary function that needs to be in your path
params.SNR = SNR;

%% Choose the data that you want to analyze based on parameters 

iSelect = setStimulusCondition(Condition);              %Condition is established in the 2nd chunk

iSelect.err = 'c'; % choose correct trials
data_c = selectData(data,params,iSelect);

iSelect.err = 'w'; % choose wrong trials
data_w = selectData(data,params,iSelect);

%% non-parametric computation of the cross-spectral density matrix (slow)
cfg           = [];
cfg.method    = 'mtmfft';
% cfg.taper     = 'hanning';
cfg.output    = 'fourier';
cfg.taper     = 'dpss'; % use this line to use >1 tapers                                                                   
cfg.tapsmofrq = 4;      % use this line to use >1 tapers
cfg.pad       = 1.0;
cfg.foilim    = [1 100];
freq_c        = ft_freqanalysis(cfg, data_c);
freq_w        = ft_freqanalysis(cfg, data_w);
fd_c          = ft_freqdescriptives(cfg,freq_c);
fd_w          = ft_freqdescriptives(cfg,freq_w);

nTrial_c = numel(data_c.trial);
nTrial_w = numel(data_w.trial);
nTaper_c = size(freq_c.fourierspctrm,1) / nTrial_c; % number of tapers
nTaper_w = size(freq_w.fourierspctrm,1) / nTrial_w; % number of tapers
nCh = numel(data_c.label); % number of channel
nFreq = numel(freq_c.freq); % number of points in freq

TrialInfo.nTrial.c = nTrial_c;
TrialInfo.nTrial.w = nTrial_w;
TrialInfo.nTaper = nTaper_c;

%% Spectra from PriorOnly Correct Trials

if strcmp(Channels,'PFCdeep_ACdeep') == 1
prior_PFC_LFP_avg_c = mean(fd_c.powspctrm(16:20,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_c  = mean(fd_c.powspctrm(36:40, :));      % average across the AC channels 
prior_PFC_LFP_std_c = std(fd_c.powspctrm(16:20,   :));      % std across the PFC channels 
prior_AC_LFP_std_c  = std(fd_c.powspctrm(36:40,  :));      % average across the AC channels 
prior_PFCupper_c    = prior_PFC_LFP_avg_c + prior_PFC_LFP_std_c;
prior_PFClower_c    = prior_PFC_LFP_avg_c - prior_PFC_LFP_std_c;
prior_ACupper_c     = prior_AC_LFP_avg_c  + prior_AC_LFP_std_c;
prior_AClower_c     = prior_AC_LFP_avg_c  - prior_AC_LFP_std_c;
end

if strcmp(Channels,'PFCupper_ACupper') == 1
prior_PFC_LFP_avg_c = mean(fd_c.powspctrm(1:5,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_c  = mean(fd_c.powspctrm(21:25, :));      % average across the AC channels 
prior_PFC_LFP_std_c = std(fd_c.powspctrm(1:5,   :));      % std across the PFC channels 
prior_AC_LFP_std_c  = std(fd_c.powspctrm(21:25,  :));      % average across the AC channels 
prior_PFCupper_c    = prior_PFC_LFP_avg_c + prior_PFC_LFP_std_c;
prior_PFClower_c    = prior_PFC_LFP_avg_c - prior_PFC_LFP_std_c;
prior_ACupper_c     = prior_AC_LFP_avg_c  + prior_AC_LFP_std_c;
prior_AClower_c     = prior_AC_LFP_avg_c  - prior_AC_LFP_std_c;
end

if strcmp(Channels,'PFClowmid_ACupmid') == 1
prior_PFC_LFP_avg_c = mean(fd_c.powspctrm(11:15,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_c  = mean(fd_c.powspctrm(26:30, :));      % average across the AC channels 
prior_PFC_LFP_std_c = std(fd_c.powspctrm(11:15,   :));      % std across the PFC channels 
prior_AC_LFP_std_c  = std(fd_c.powspctrm(26:30,  :));      % average across the AC channels 
prior_PFCupper_c    = prior_PFC_LFP_avg_c + prior_PFC_LFP_std_c;
prior_PFClower_c    = prior_PFC_LFP_avg_c - prior_PFC_LFP_std_c;
prior_ACupper_c     = prior_AC_LFP_avg_c  + prior_AC_LFP_std_c;
prior_AClower_c     = prior_AC_LFP_avg_c  - prior_AC_LFP_std_c;
end


if strcmp(Channels,'PFCupmid_ACupmid') == 1
prior_PFC_LFP_avg_c = mean(fd_c.powspctrm(6:10,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_c  = mean(fd_c.powspctrm(26:30, :));      % average across the AC channels 
prior_PFC_LFP_std_c = std(fd_c.powspctrm(6:10,   :));      % std across the PFC channels 
prior_AC_LFP_std_c  = std(fd_c.powspctrm(26:30,  :));      % average across the AC channels 
prior_PFCupper_c    = prior_PFC_LFP_avg_c + prior_PFC_LFP_std_c;
prior_PFClower_c    = prior_PFC_LFP_avg_c - prior_PFC_LFP_std_c;
prior_ACupper_c     = prior_AC_LFP_avg_c  + prior_AC_LFP_std_c;
prior_AClower_c     = prior_AC_LFP_avg_c  - prior_AC_LFP_std_c;
end

if strcmp(Channels,'PFCdeep_ACupmid') == 1
prior_PFC_LFP_avg_c = mean(fd_c.powspctrm(16:20,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_c  = mean(fd_c.powspctrm(26:30, :));      % average across the AC channels 
prior_PFC_LFP_std_c = std(fd_c.powspctrm(16:20,   :));      % std across the PFC channels 
prior_AC_LFP_std_c  = std(fd_c.powspctrm(26:30,  :));      % average across the AC channels 
prior_PFCupper_c    = prior_PFC_LFP_avg_c + prior_PFC_LFP_std_c;
prior_PFClower_c    = prior_PFC_LFP_avg_c - prior_PFC_LFP_std_c;
prior_ACupper_c     = prior_AC_LFP_avg_c  + prior_AC_LFP_std_c;
prior_AClower_c     = prior_AC_LFP_avg_c  - prior_AC_LFP_std_c;
end

%% Spectra from Wrong OnlyPrior Trials 

if strcmp(Channels,'PFCdeep_ACdeep') == 1
prior_PFC_LFP_avg_w = mean(fd_w.powspctrm(16:20,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_w  = mean(fd_w.powspctrm(36:40, :));      % average across the AC channels 
prior_PFC_LFP_std_w = std(fd_w.powspctrm(16:20,   :));      % std across the PFC channels 
prior_AC_LFP_std_w  = std(fd_w.powspctrm(36:40,  :));      % average across the AC channels 
prior_PFCupper_w    = prior_PFC_LFP_avg_w + prior_PFC_LFP_std_w;
prior_PFClower_w    = prior_PFC_LFP_avg_w - prior_PFC_LFP_std_w;
prior_ACupper_w     = prior_AC_LFP_avg_w  + prior_AC_LFP_std_w;
prior_AClower_w     = prior_AC_LFP_avg_w  - prior_AC_LFP_std_w;
end

if strcmp(Channels,'PFCupper_ACupper') == 1
prior_PFC_LFP_avg_w = mean(fd_w.powspctrm(1:5,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_w  = mean(fd_w.powspctrm(21:25, :));      % average across the AC channels 
prior_PFC_LFP_std_w = std(fd_w.powspctrm(1:5,   :));      % std across the PFC channels 
prior_AC_LFP_std_w  = std(fd_w.powspctrm(21:25,  :));      % average across the AC channels 
prior_PFCupper_w    = prior_PFC_LFP_avg_w + prior_PFC_LFP_std_w;
prior_PFClower_w    = prior_PFC_LFP_avg_w - prior_PFC_LFP_std_w;
prior_ACupper_w     = prior_AC_LFP_avg_w  + prior_AC_LFP_std_w;
prior_AClower_w     = prior_AC_LFP_avg_w  - prior_AC_LFP_std_w;
end

if strcmp(Channels,'PFClowmid_ACupmid') == 1
prior_PFC_LFP_avg_w = mean(fd_w.powspctrm(11:15,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_w  = mean(fd_w.powspctrm(26:30, :));      % average across the AC channels 
prior_PFC_LFP_std_w = std(fd_w.powspctrm(11:15,   :));      % std across the PFC channels 
prior_AC_LFP_std_w  = std(fd_w.powspctrm(26:30,  :));      % average across the AC channels 
prior_PFCupper_w    = prior_PFC_LFP_avg_w + prior_PFC_LFP_std_w;
prior_PFClower_w    = prior_PFC_LFP_avg_w - prior_PFC_LFP_std_w;
prior_ACupper_w     = prior_AC_LFP_avg_w  + prior_AC_LFP_std_w;
prior_AClower_w     = prior_AC_LFP_avg_w  - prior_AC_LFP_std_w;
end

if strcmp(Channels,'PFCupmid_ACupmid') == 1
prior_PFC_LFP_avg_w = mean(fd_w.powspctrm(6:10,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_w  = mean(fd_w.powspctrm(26:30, :));      % average across the AC channels 
prior_PFC_LFP_std_w = std(fd_w.powspctrm(6:10,   :));      % std across the PFC channels 
prior_AC_LFP_std_w  = std(fd_w.powspctrm(26:30,  :));      % average across the AC channels 
prior_PFCupper_w    = prior_PFC_LFP_avg_w + prior_PFC_LFP_std_w;
prior_PFClower_w    = prior_PFC_LFP_avg_w - prior_PFC_LFP_std_w;
prior_ACupper_w     = prior_AC_LFP_avg_w  + prior_AC_LFP_std_w;
prior_AClower_w     = prior_AC_LFP_avg_w  - prior_AC_LFP_std_w;
end

if strcmp(Channels,'PFCdeep_ACupmid') == 1
prior_PFC_LFP_avg_w = mean(fd_w.powspctrm(16:20,  :));      % average across the PFC channels; 
prior_AC_LFP_avg_w  = mean(fd_w.powspctrm(26:30, :));      % average across the AC channels 
prior_PFC_LFP_std_w = std(fd_w.powspctrm(16:20,   :));      % std across the PFC channels 
prior_AC_LFP_std_w  = std(fd_w.powspctrm(26:30,  :));      % average across the AC channels 
prior_PFCupper_w    = prior_PFC_LFP_avg_w + prior_PFC_LFP_std_w;
prior_PFClower_w    = prior_PFC_LFP_avg_w - prior_PFC_LFP_std_w;
prior_ACupper_w     = prior_AC_LFP_avg_w  + prior_AC_LFP_std_w;
prior_AClower_w     = prior_AC_LFP_avg_w  - prior_AC_LFP_std_w;
end

%% Define data to analyze

Condition = 'OnlyPretone';        % options: 'OnlyPrior', 'OnlyPretone';

%% load fieldtrip data formatted by bipolarLFP.m

fName = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
%fName = strcat(Animal,'-',RecDate,'_bipolarLFP_',Epoch);
load(fName);

%% set parameters required for fieldtrip functions and data selection
%Used to for selectData function

params.choice = choice;
params.err = err;
params.pretone = pretone;
params.pretoneLength = pretoneLength;
params.prior = cell2char(prior); %cell2char is an auxillary function that needs to be in your path
params.SNR = SNR;

%% Choose the data that you want to analyze based on parameters 

iSelect = setStimulusCondition(Condition);              %Condition is established in the 2nd chunk

iSelect.err = 'c'; % choose correct trials
data_c = selectData(data,params,iSelect);

iSelect.err = 'w'; % choose wrong trials
data_w = selectData(data,params,iSelect);

%% non-parametric computation of the cross-spectral density matrix (slow)
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'hanning';
cfg.output    = 'fourier';
%cfg.taper     = 'dpss'; % use this line to use >1 tapers                                                                   
%cfg.tapsmofrq = 4;      % use this line to use >1 tapers
cfg.pad       = 1.0;
cfg.foilim    = [1 100];
freq_c        = ft_freqanalysis(cfg, data_c);
freq_w        = ft_freqanalysis(cfg, data_w);
fd_c          = ft_freqdescriptives(cfg,freq_c);
fd_w          = ft_freqdescriptives(cfg,freq_w);

nTrial_c = numel(data_c.trial);
nTrial_w = numel(data_w.trial);
nTaper_c = size(freq_c.fourierspctrm,1) / nTrial_c; % number of tapers
nTaper_w = size(freq_w.fourierspctrm,1) / nTrial_w; % number of tapers
nCh = numel(data_c.label); % number of channel
nFreq = numel(freq_c.freq); % number of points in freq

TrialInfo.nTrial.c = nTrial_c;
TrialInfo.nTrial.w = nTrial_w;
TrialInfo.nTaper = nTaper_c;

%% Spectra from OnlyPretone Correct Trials

if strcmp(Channels,'PFCdeep_ACdeep') == 1
pretone_PFC_LFP_avg_c = mean(fd_c.powspctrm(16:20,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_c  = mean(fd_c.powspctrm(36:40, :));      % average across the AC channels 
pretone_PFC_LFP_std_c = std(fd_c.powspctrm(16:20,   :));      % std across the PFC channels 
pretone_AC_LFP_std_c  = std(fd_c.powspctrm(36:40,  :));      % average across the AC channels 
pretone_PFCupper_c    = pretone_PFC_LFP_avg_c + pretone_PFC_LFP_std_c;
pretone_PFClower_c    = pretone_PFC_LFP_avg_c - pretone_PFC_LFP_std_c;
pretone_ACupper_c     = pretone_AC_LFP_avg_c  + pretone_AC_LFP_std_c;
pretone_AClower_c     = pretone_AC_LFP_avg_c  - pretone_AC_LFP_std_c;
end

if strcmp(Channels,'PFCupper_ACupper') == 1
pretone_PFC_LFP_avg_c = mean(fd_c.powspctrm(1:5,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_c  = mean(fd_c.powspctrm(21:25, :));      % average across the AC channels 
pretone_PFC_LFP_std_c = std(fd_c.powspctrm(1:5,   :));      % std across the PFC channels 
pretone_AC_LFP_std_c  = std(fd_c.powspctrm(21:25,  :));      % average across the AC channels 
pretone_PFCupper_c    = pretone_PFC_LFP_avg_c + pretone_PFC_LFP_std_c;
pretone_PFClower_c    = pretone_PFC_LFP_avg_c - pretone_PFC_LFP_std_c;
pretone_ACupper_c     = pretone_AC_LFP_avg_c  + pretone_AC_LFP_std_c;
pretone_AClower_c     = pretone_AC_LFP_avg_c  - pretone_AC_LFP_std_c;
end

if strcmp(Channels,'PFClowmid_ACupmid') == 1
pretone_PFC_LFP_avg_c = mean(fd_c.powspctrm(11:15,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_c  = mean(fd_c.powspctrm(26:30, :));      % average across the AC channels 
pretone_PFC_LFP_std_c = std(fd_c.powspctrm(11:15,   :));      % std across the PFC channels 
pretone_AC_LFP_std_c  = std(fd_c.powspctrm(26:30,  :));      % average across the AC channels 
pretone_PFCupper_c    = pretone_PFC_LFP_avg_c + pretone_PFC_LFP_std_c;
pretone_PFClower_c    = pretone_PFC_LFP_avg_c - pretone_PFC_LFP_std_c;
pretone_ACupper_c     = pretone_AC_LFP_avg_c  + pretone_AC_LFP_std_c;
pretone_AClower_c     = pretone_AC_LFP_avg_c  - pretone_AC_LFP_std_c;
end

if strcmp(Channels,'PFCupmid_ACupmid') == 1
pretone_PFC_LFP_avg_c = mean(fd_c.powspctrm(6:10,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_c  = mean(fd_c.powspctrm(26:30, :));      % average across the AC channels 
pretone_PFC_LFP_std_c = std(fd_c.powspctrm(6:10,   :));      % std across the PFC channels 
pretone_AC_LFP_std_c  = std(fd_c.powspctrm(26:30,  :));      % average across the AC channels 
pretone_PFCupper_c    = pretone_PFC_LFP_avg_c + pretone_PFC_LFP_std_c;
pretone_PFClower_c    = pretone_PFC_LFP_avg_c - pretone_PFC_LFP_std_c;
pretone_ACupper_c     = pretone_AC_LFP_avg_c  + pretone_AC_LFP_std_c;
pretone_AClower_c     = pretone_AC_LFP_avg_c  - pretone_AC_LFP_std_c;
end

if strcmp(Channels,'PFCdeep_ACupmid') == 1
pretone_PFC_LFP_avg_c = mean(fd_c.powspctrm(16:20,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_c  = mean(fd_c.powspctrm(26:30, :));      % average across the AC channels 
pretone_PFC_LFP_std_c = std(fd_c.powspctrm(16:20,   :));      % std across the PFC channels 
pretone_AC_LFP_std_c  = std(fd_c.powspctrm(26:30,  :));      % average across the AC channels 
pretone_PFCupper_c    = pretone_PFC_LFP_avg_c + pretone_PFC_LFP_std_c;
pretone_PFClower_c    = pretone_PFC_LFP_avg_c - pretone_PFC_LFP_std_c;
pretone_ACupper_c     = pretone_AC_LFP_avg_c  + pretone_AC_LFP_std_c;
pretone_AClower_c     = pretone_AC_LFP_avg_c  - pretone_AC_LFP_std_c;
end

%% Spectra from Wrong OnlyPretone Trials 

if strcmp(Channels,'PFCdeep_ACdeep') == 1
pretone_PFC_LFP_avg_w = mean(fd_w.powspctrm(16:20,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_w  = mean(fd_w.powspctrm(36:40, :));      % average across the AC channels 
pretone_PFC_LFP_std_w = std(fd_w.powspctrm(16:20,   :));      % std across the PFC channels 
pretone_AC_LFP_std_w  = std(fd_w.powspctrm(36:40,  :));      % average across the AC channels 
pretone_PFCupper_w    = pretone_PFC_LFP_avg_w + pretone_PFC_LFP_std_w;
pretone_PFClower_w    = pretone_PFC_LFP_avg_w - pretone_PFC_LFP_std_w;
pretone_ACupper_w     = pretone_AC_LFP_avg_w  + pretone_AC_LFP_std_w;
pretone_AClower_w     = pretone_AC_LFP_avg_w  - pretone_AC_LFP_std_w;
end

if strcmp(Channels,'PFCupper_ACupper') == 1
pretone_PFC_LFP_avg_w = mean(fd_w.powspctrm(1:5,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_w  = mean(fd_w.powspctrm(21:25, :));      % average across the AC channels 
pretone_PFC_LFP_std_w = std(fd_w.powspctrm(1:5,   :));      % std across the PFC channels 
pretone_AC_LFP_std_w  = std(fd_w.powspctrm(21:25,  :));      % average across the AC channels 
pretone_PFCupper_w    = pretone_PFC_LFP_avg_w + pretone_PFC_LFP_std_w;
pretone_PFClower_w    = pretone_PFC_LFP_avg_w - pretone_PFC_LFP_std_w;
pretone_ACupper_w     = pretone_AC_LFP_avg_w  + pretone_AC_LFP_std_w;
pretone_AClower_w     = pretone_AC_LFP_avg_w  - pretone_AC_LFP_std_w;
end

if strcmp(Channels,'PFClowmid_ACupmid') == 1
pretone_PFC_LFP_avg_w = mean(fd_w.powspctrm(11:15,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_w  = mean(fd_w.powspctrm(26:30, :));      % average across the AC channels 
pretone_PFC_LFP_std_w = std(fd_w.powspctrm(11:15,   :));      % std across the PFC channels 
pretone_AC_LFP_std_w  = std(fd_w.powspctrm(26:30,  :));      % average across the AC channels 
pretone_PFCupper_w    = pretone_PFC_LFP_avg_w + pretone_PFC_LFP_std_w;
pretone_PFClower_w    = pretone_PFC_LFP_avg_w - pretone_PFC_LFP_std_w;
pretone_ACupper_w     = pretone_AC_LFP_avg_w  + pretone_AC_LFP_std_w;
pretone_AClower_w     = pretone_AC_LFP_avg_w  - pretone_AC_LFP_std_w;
end

if strcmp(Channels,'PFCupmid_ACupmid') == 1
pretone_PFC_LFP_avg_w = mean(fd_w.powspctrm(6:10,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_w  = mean(fd_w.powspctrm(26:30, :));      % average across the AC channels 
pretone_PFC_LFP_std_w = std(fd_w.powspctrm(6:10,   :));      % std across the PFC channels 
pretone_AC_LFP_std_w  = std(fd_w.powspctrm(26:30,  :));      % average across the AC channels 
pretone_PFCupper_w    = pretone_PFC_LFP_avg_w + pretone_PFC_LFP_std_w;
pretone_PFClower_w    = pretone_PFC_LFP_avg_w - pretone_PFC_LFP_std_w;
pretone_ACupper_w     = pretone_AC_LFP_avg_w  + pretone_AC_LFP_std_w;
pretone_AClower_w     = pretone_AC_LFP_avg_w  - pretone_AC_LFP_std_w;
end

if strcmp(Channels,'PFCdeep_ACupmid') == 1
pretone_PFC_LFP_avg_w = mean(fd_w.powspctrm(16:20,  :));      % average across the PFC channels; 
pretone_AC_LFP_avg_w  = mean(fd_w.powspctrm(26:30, :));      % average across the AC channels 
pretone_PFC_LFP_std_w = std(fd_w.powspctrm(16:20,   :));      % std across the PFC channels 
pretone_AC_LFP_std_w  = std(fd_w.powspctrm(26:30,  :));      % average across the AC channels 
pretone_PFCupper_w    = pretone_PFC_LFP_avg_w + pretone_PFC_LFP_std_w;
pretone_PFClower_w    = pretone_PFC_LFP_avg_w - pretone_PFC_LFP_std_w;
pretone_ACupper_w     = pretone_AC_LFP_avg_w  + pretone_AC_LFP_std_w;
pretone_AClower_w     = pretone_AC_LFP_avg_w  - pretone_AC_LFP_std_w;
end

%%


% Sample data
freq_range = 1:100;

% Plot mean line
figure;
sgtitle('CORRECT')

subplot(2,1,1)
hold on
ax1 = gca;  % Get the current axes
set(ax1, 'Color', [0.68, 0.85, 0.90]);  % Light Blue background
plot(freq_range, prior_PFC_LFP_avg_c, 'LineWidth', 2, 'Color', 'b');
plot(freq_range, prior_AC_LFP_avg_c, 'LineWidth', 2, 'Color', 'r');
% Fill area between upper and lower bounds
fill([freq_range, fliplr(freq_range)], [prior_PFCupper_c, fliplr(prior_PFClower_c)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([freq_range, fliplr(freq_range)], [prior_ACupper_c,  fliplr(prior_AClower_c)],  'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off
% Customize plot
legend('','', 'PFC' ,'AC')
xlabel('Frequency (Hz)');
ylabel('LFP Mean Power');
title('LED')

subplot(2,1,2)
hold on
ax2 = gca;  % Get the current axes
set(ax2, 'Color', [0.56, 0.93, 0.56]);  % Light Green background
plot(freq_range, pretone_PFC_LFP_avg_c, 'LineWidth', 2, 'Color', 'b');
plot(freq_range, pretone_AC_LFP_avg_c, 'LineWidth', 2, 'Color', 'r');
% Fill area between upper and lower bounds
fill([freq_range, fliplr(freq_range)], [pretone_PFCupper_c, fliplr(pretone_PFClower_c)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([freq_range, fliplr(freq_range)], [pretone_ACupper_c,  fliplr(pretone_AClower_c)],  'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off
% Customize plot
legend('','', 'PFC' ,'AC')
xlabel('Frequency (Hz)');
ylabel('LFP Mean Power');
title('Pretone')


%%

% Sample data
freq_range = 1:100;

% Plot mean line
figure;
sgtitle('WRONG')

subplot(2,1,1)
hold on
ax1 = gca;  % Get the current axes
set(ax1, 'Color', [0.68, 0.85, 0.90]);  % Light Blue background
plot(freq_range, prior_PFC_LFP_avg_w, 'LineWidth', 2, 'Color', 'b');
plot(freq_range, prior_AC_LFP_avg_w, 'LineWidth', 2, 'Color', 'r');
% Fill area between upper and lower bounds
fill([freq_range, fliplr(freq_range)], [prior_PFCupper_w, fliplr(prior_PFClower_w)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([freq_range, fliplr(freq_range)], [prior_ACupper_w,  fliplr(prior_AClower_w)],  'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off
% Customize plot
legend('','', 'PFC' ,'AC')
xlabel('Frequency (Hz)');

ylabel('LFP Mean Power');
title('LED')

subplot(2,1,2)
hold on
ax2 = gca;  % Get the current axes
set(ax2, 'Color', [0.56, 0.93, 0.56]);  % Light Green background
plot(freq_range, pretone_PFC_LFP_avg_w, 'LineWidth', 2, 'Color', 'b');
plot(freq_range, pretone_AC_LFP_avg_w, 'LineWidth', 2, 'Color', 'r');
% Fill area between upper and lower bounds
fill([freq_range, fliplr(freq_range)], [pretone_PFCupper_w, fliplr(pretone_PFClower_w)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([freq_range, fliplr(freq_range)], [pretone_ACupper_w,  fliplr(pretone_AClower_w)],  'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off
% Customize plot
legend('','', 'PFC' ,'AC')
xlabel('Frequency (Hz)');
ylabel('LFP Mean Power');
title('Pretone')