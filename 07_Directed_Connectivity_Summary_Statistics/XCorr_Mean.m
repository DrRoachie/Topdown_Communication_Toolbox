

% Define data

Epoch               = 'testToneOnset';
Current_ChanPair    = {'*PFC_ch13' '*AC_ch03'};

rootdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\05_Epoc_Cut_2';
animals = {'MrCassius', 'MrM'};

OnlyPrior_XCorr_Array       = [];
OnlyPrior_XCorr_Upper_Array = [];
OnlyPrior_XCorr_Lower_Array = [];

for j = 1:length(animals)

    Animal = animals{j};

    sessions = dir(fullfile(rootdir, Animal, extractBefore(Epoch, 'Onset'), '19*'));

    for i = 1:length(sessions)
        
        RecDate = sessions(i).name;
        disp(['Now processing ' fullfile(rootdir, Animal, extractBefore(Epoch, 'Onset'),RecDate)]);

        [OnlyPrior_shuffled_xcorr_result, OnlyPrior_shuffled_lagsResults,OnlyPrior_xcorr_result, OnlyPrior_lagsResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, 'OnlyPrior', Current_ChanPair);
        OnlyPrior_XCorr_Array       = [OnlyPrior_XCorr_Array; OnlyPrior_xcorr_result.avg];
        OnlyPrior_XCorr_Upper_Array = [OnlyPrior_XCorr_Upper_Array; OnlyPrior_xcorr_result.upper];
        OnlyPrior_XCorr_Lower_Array = [OnlyPrior_XCorr_Lower_Array; OnlyPrior_xcorr_result.lower];
    end

end

Mean_XCorr          = mean(OnlyPrior_XCorr_Array);
Mean_XCorr_Upper    = mean(OnlyPrior_XCorr_Upper_Array);
Mean_XCorr_Lower    = mean(OnlyPrior_XCorr_Lower_Array);


%  plot xcorr results  
figure;
hold on;
plot(OnlyPrior_lagsResults, OnlyPrior_xcorr_result.avg, 'g', 'LineWidth', 1.5);
fill([OnlyPrior_lagsResults, fliplr(OnlyPrior_lagsResults)], [OnlyPrior_xcorr_result.upper, fliplr(OnlyPrior_xcorr_result.lower)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');