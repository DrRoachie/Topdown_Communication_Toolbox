%% Define data to analyze
Animal     = 'MrMiyagi';                                        % options: 'MrCassius' & 'MrMiyagi';
Epoch      = 'moveOnset';                                       % options: 'LED', 'testTone';

%% Establish file and save directories
datadir = ['/home/arl/Documents/DATA/02_Preprocessed/' Animal '/' Epoch '/'];
savedir = '/media/arl/SHD/07_MoveOnset_Cut/';
addpath(genpath(datadir));

sessions = dir(fullfile(datadir,'*.mat'));

if strcmp(Animal, 'MrMiyagi')
    Animal = 'MrM';
end

for k = 1:length(sessions)
    % Load in data
    RecDate = regexp(sessions(k).name, '\d{6}', 'match');
    RecDate = RecDate{1};
    fName   = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft.mat');
    load(fName);
    
    % Find index of 50ms after onset
    Time_Array  = data.time{1};
    tolerance   = 1e-4;
    Cut_Index   = find(abs(Time_Array - 0.0500) < tolerance);
    
    % Loop through all trials and cut data.trial to 50ms after onset
    for i = 1:length(data.trial)
        data.trial{i} = data.trial{i}(:,1:Cut_Index);
        data.time{i} = data.time{i}(1:Cut_Index);
    end
    
    % Save cut data
    savefilename = sprintf('%s_%s_%s_cut.mat', Animal, Epoch, RecDate);

    if ~exist(fullfile(savedir, Animal, Epoch, RecDate), 'dir')
        mkdir(fullfile(savedir, Animal, Epoch, RecDate));
    end

    save(fullfile(savedir, Animal, Epoch, RecDate, savefilename),'choice','data','err','info','pretone','pretoneLength','prior','SNR','stim','trial_id');
    disp(['File saved: ' savefilename ' cut to index ' num2str(Cut_Index)]);
end