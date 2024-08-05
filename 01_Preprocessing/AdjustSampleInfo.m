% Get a list of all session directories
parentDir   = 'C:\Users\Corey Roach\Documents\00_DATA\04_Epoc_Cut\MrCassius\testTone';
sessionDirs = dir('C:\Users\Corey Roach\Documents\00_DATA\04_Epoc_Cut\MrCassius\testTone');
SAVE_DIR    = 'C:\Users\Corey Roach\Documents\00_DATA\Updated_Epoc_Cut\Cassius\testTone';
Epoch       = 'testToneOnset';
 

% Loop through each session directory
for t = 1:length(sessionDirs)

    if sessionDirs(t).isdir && ~strcmp(sessionDirs(t).name, '.') && ~strcmp(sessionDirs(t).name, '..')
        sessionPath = fullfile(parentDir, sessionDirs(t).name);
        
        % Load the .mat file in the session directory
        matFiles = dir(fullfile(sessionPath, '*.mat'));
        
            if ~isempty(matFiles)
            % Load the .mat file
            matFileName = matFiles(1).name; % Assuming there's only one .mat file
            filePath = fullfile(sessionPath, matFileName);
            matData = load(filePath);

                if strcmp(Epoch,'testToneOnset')||strcmp(Epoch,'PreToneEpoch') == 1 
                 
                        % Define increment size
                        duration = 200;
                        intertrial_interval = 5;  
                        
                        % Calculate new sampleinfo
                        numTrials     = size(data.sampleinfo, 1);
                        newSampleinfo = zeros(numTrials, 2);
                
                        for t = 1:numTrials
                            if t == 1
                                newSampleinfo(t, 1) = 1;
                            else
                                newSampleinfo(t, 1) = newSampleinfo(t-1, 2) + intertrial_interval + 1;
                            end
                            
                            newSampleinfo(t, 2) = newSampleinfo(t, 1) + duration - 1;

                            data.sampleinfo = newSampleinfo;

                            save_file_name = fName;
                            save(fullfile(SAVE_DIR,save_file_name),  'data');
                            
                            % Clear loaded data to free up memory
                            clear matData;

                        end

                    end
                end  
            end      
    end 

     