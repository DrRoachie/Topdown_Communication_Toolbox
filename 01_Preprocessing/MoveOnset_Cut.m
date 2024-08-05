
Animal = 'MrM';

datadir = 'E:\07_MoveOnset_Cut';
sessions = dir(fullfile(datadir, Animal, 'moveOnset', '19*'));

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    file_name = [Animal '_moveOnset_' RecDate '_cut.mat'];
    disp(['Now processing ' file_name]);
    load(fullfile(datadir, Animal, 'moveOnset', RecDate, file_name), 'data');

    start_idx = 1;
    for j = 1:length(data.sampleinfo)
        data.sampleinfo(j,1) = start_idx;
        data.sampleinfo(j,2) = data.sampleinfo(j,1) + 550;
        start_idx = data.sampleinfo(j,2) + 5;
    end

    disp('Saving...');
    save(fullfile(datadir, Animal, 'moveOnset', RecDate, file_name), 'data', '-append');
end