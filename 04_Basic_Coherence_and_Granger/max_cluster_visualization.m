%%

% Define the parent directory
parentDir = 'C:\Users\Corey Roach\Desktop\2024_01_17_Analysis\testTone\PFCdeep_ACdeep';

% List subdirectories in the parent directory (should be MrCassius and MrC)
subDirs = dir(parentDir);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name}, {'.', '..'}));

% Initialize an empty table
outputTable = table();

% Loop through each primary subdirectory
for i = 1:length(subDirs)
    primarySubDir = subDirs(i).name;
    
    % List secondary subdirectories in the primary subdirectory
    secondarySubDirs = dir(fullfile(parentDir, primarySubDir));
    secondarySubDirs = secondarySubDirs([secondarySubDirs.isdir] & ~ismember({secondarySubDirs.name}, {'.', '..'}));

    % Loop through each secondary subdirectory
    for j = 1:length(secondarySubDirs)
        secondarySubDir = secondarySubDirs(j).name;
        addpath(genpath(secondarySubDir));

        % Construct the full path to the .mat file
        matFilePath = fullfile(parentDir, primarySubDir, secondarySubDir, 'Statistic_Test_Results.mat');

        % Load the .mat file
        loadedData = load(matFilePath);

        % Extract the desired variables (replace 'your_variable_names' accordingly)
        band         = fieldnames(loadedData.pvalue_C);
        pvalue_c     = struct2cell(loadedData.pvalue_C);
        pvalue_c     = cell2mat(pvalue_c);
        pvalue_w     = struct2cell(loadedData.pvalue_W);
        pvalue_w     = cell2mat(pvalue_w);
        clustermax_c = struct2cell(loadedData.data_max_cluster_sum_Maris_C);
        clustermax_c = cell2mat(clustermax_c);
        clustermax_w = struct2cell(loadedData.data_max_cluster_sum_Maris_W);
        clustermax_w = cell2mat(clustermax_w);
        animal       = repmat({subDirs(i).name}, length(band), 1);
        %animal       = categorical(animal);
        session      = repmat({secondarySubDirs(j).name}, length(band), 1);
        session      = cell2mat(session);
        session      = datetime(num2str(session), 'InputFormat', 'yyyyMMdd', 'Format', 'yyyy-MM-dd');
        sig_c        = pvalue_c <= .05; 
        sig_w        = pvalue_w <= .05;

        % Create a row for the output table
        newRows = table(session, animal, band, pvalue_c, clustermax_c, sig_c, pvalue_w, clustermax_w, sig_w);

        % Append the rows to the output table

        % Append the row to the output table
        outputTable = [outputTable; newRows];

    end
end

%% clear all but output table 

% List of variables to retain
retainVars = {'outputTable'};

% Get a list of all variables in the workspace
allVars = evalin('base', 'who');
% Identify variables to clear
varsToClear = setdiff(allVars, retainVars);
% Clear the unwanted variables
for i = 1:length(varsToClear)
    evalin('base', ['clear ' varsToClear{i}]);
end

%%

% Filter for theta band
theta_idx           = strcmp(outputTable.band, 'theta');
theta_animal        = outputTable.animal(theta_idx);
theta_clustermax_c  = outputTable.clustermax_c(theta_idx);
theta_clustermax_w  = outputTable.clustermax_w(theta_idx);
theta_sig_c         = outputTable.sig_c(theta_idx);
theta_sig_w         = outputTable.sig_w(theta_idx);

alpha_idx           = strcmp(outputTable.band, 'alpha');
alpha_animal        = outputTable.animal(alpha_idx);
alpha_clustermax_c  = outputTable.clustermax_c(alpha_idx);
alpha_clustermax_w  = outputTable.clustermax_w(alpha_idx);
alpha_sig_c         = outputTable.sig_c(alpha_idx);
alpha_sig_w         = outputTable.sig_w(alpha_idx);

beta_idx           = strcmp(outputTable.band, 'beta');
beta_animal        = outputTable.animal(beta_idx);
beta_clustermax_c  = outputTable.clustermax_c(beta_idx);
beta_clustermax_w  = outputTable.clustermax_w(beta_idx);
beta_sig_c         = outputTable.sig_c(beta_idx);
beta_sig_w         = outputTable.sig_w(beta_idx);

gamma_idx           = strcmp(outputTable.band, 'gamma');
gamma_animal        = outputTable.animal(gamma_idx);
gamma_clustermax_c  = outputTable.clustermax_c(gamma_idx);
gamma_clustermax_w  = outputTable.clustermax_w(gamma_idx);
gamma_sig_c         = outputTable.sig_c(gamma_idx);
gamma_sig_w         = outputTable.sig_w(gamma_idx);

highGamma_idx           = strcmp(outputTable.band, 'highGamma');
highGamma_animal        = outputTable.animal(highGamma_idx);
highGamma_clustermax_c  = outputTable.clustermax_c(highGamma_idx);
highGamma_clustermax_w  = outputTable.clustermax_w(highGamma_idx);
highGamma_sig_c         = outputTable.sig_c(highGamma_idx);
highGamma_sig_w         = outputTable.sig_w(highGamma_idx);

all_idx           = strcmp(outputTable.band, 'all');
all_animal        = outputTable.animal(all_idx);
all_clustermax_c  = outputTable.clustermax_c(all_idx);
all_clustermax_w  = outputTable.clustermax_w(all_idx);
all_sig_c         = outputTable.sig_c(all_idx);
all_sig_w         = outputTable.sig_w(all_idx);

%% Plotting
% Constants
x_categories = [1, 2]; % Representing clustermax_c and clustermax_w
animal_shapes = {'o', 's', 'd', '^', 'v', '<', '>', 'p', 'h'}; % Different shapes for different animals
unique_animals = unique(outputTable.animal);
marker_size = 100; % Adjustable marker size
line_width = 1; % Adjustable line width for better visibility

% Prepare figure
figure;

subplot(2,3,1)
hold on;
% Overlay individual data points
for i = 1:length(theta_clustermax_c)
    animal_idx = find(strcmp(unique_animals, theta_animal{i})); % Find the index of the animal for shape
    color_c = 'black';
    color_w = 'black';
    if theta_sig_c(i) == 1
        color_c = 'green';
    end
    if theta_sig_w(i) == 1
        color_w = 'green';
    end
    % Plot individual points with filled color
    scatter(x_categories(1), theta_clustermax_c(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_c, 'MarkerFaceColor', color_c, 'MarkerFaceAlpha', 0.35);
    scatter(x_categories(2), theta_clustermax_w(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_w, 'MarkerFaceColor', color_w, 'MarkerFaceAlpha', 0.35);
    
    % Draw lines connecting the points
    line([x_categories(1), x_categories(2)], [theta_clustermax_c(i), theta_clustermax_w(i)], 'Color', 'black', 'LineWidth', line_width);
end
% Customize the plot
xlim([0.5, 2.5]); % Adjust as needed for better spacing
set(gca, 'xtick', x_categories, 'xticklabel', {'CORRECT', 'WRONG'});
ylabel('Z Max Cluster Value');
title('theta');

hold off;

subplot(2,3,2)
hold on;
% Overlay individual data points
for i = 1:length(alpha_clustermax_c)
    animal_idx = find(strcmp(unique_animals, alpha_animal{i})); % Find the index of the animal for shape
    color_c = 'black';
    color_w = 'black';
    if alpha_sig_c(i) == 1
        color_c = 'green';
    end
    if alpha_sig_w(i) == 1
        color_w = 'green';
    end
    % Plot individual points with filled color
    scatter(x_categories(1), alpha_clustermax_c(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_c, 'MarkerFaceColor', color_c, 'MarkerFaceAlpha', 0.35);
    scatter(x_categories(2), alpha_clustermax_w(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_w, 'MarkerFaceColor', color_w, 'MarkerFaceAlpha', 0.35);
    
    % Draw lines connecting the points
    line([x_categories(1), x_categories(2)], [alpha_clustermax_c(i), alpha_clustermax_w(i)], 'Color', 'black', 'LineWidth', line_width);
end
% Customize the plot
xlim([0.5, 2.5]); % Adjust as needed for better spacing
set(gca, 'xtick', x_categories, 'xticklabel', {'CORRECT', 'WRONG'});
ylabel('Z Max Cluster Value');
title('alpha');
hold off;

subplot(2,3,3)
hold on;
% Overlay individual data points
for i = 1:length(beta_clustermax_c)
    animal_idx = find(strcmp(unique_animals, beta_animal{i})); % Find the index of the animal for shape
    color_c = 'black';
    color_w = 'black';
    if beta_sig_c(i) == 1
        color_c = 'green';
    end
    if beta_sig_w(i) == 1
        color_w = 'green';
    end
    % Plot individual points with filled color
    scatter(x_categories(1), beta_clustermax_c(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_c, 'MarkerFaceColor', color_c, 'MarkerFaceAlpha', 0.35);
    scatter(x_categories(2), beta_clustermax_w(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_w, 'MarkerFaceColor', color_w, 'MarkerFaceAlpha', 0.35);
    
    % Draw lines connecting the points
    line([x_categories(1), x_categories(2)], [beta_clustermax_c(i), beta_clustermax_w(i)], 'Color', 'black', 'LineWidth', line_width);
end
% Customize the plot
xlim([0.5, 2.5]); % Adjust as needed for better spacing
set(gca, 'xtick', x_categories, 'xticklabel', {'CORRECT', 'WRONG'});
ylabel('Z Max Cluster Value');
title('beta');
hold off;

subplot(2,3,4)
hold on;
% Overlay individual data points
for i = 1:length(gamma_clustermax_c)
    animal_idx = find(strcmp(unique_animals, gamma_animal{i})); % Find the index of the animal for shape
    color_c = 'black';
    color_w = 'black';
    if gamma_sig_c(i) == 1
        color_c = 'green';
    end
    if gamma_sig_w(i) == 1
        color_w = 'green';
    end
    % Plot individual points with filled color
    scatter(x_categories(1), gamma_clustermax_c(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_c, 'MarkerFaceColor', color_c, 'MarkerFaceAlpha', 0.35);
    scatter(x_categories(2), gamma_clustermax_w(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_w, 'MarkerFaceColor', color_w, 'MarkerFaceAlpha', 0.35);
    
    % Draw lines connecting the points
    line([x_categories(1), x_categories(2)], [gamma_clustermax_c(i), gamma_clustermax_w(i)], 'Color', 'black', 'LineWidth', line_width);
end
% Customize the plot
xlim([0.5, 2.5]); % Adjust as needed for better spacing
set(gca, 'xtick', x_categories, 'xticklabel', {'CORRECT', 'WRONG'});
ylabel('Z Max Cluster Value');
title('gamma');
hold off;

subplot(2,3,5)
hold on;
% Overlay individual data points
for i = 1:length(highGamma_clustermax_c)
    animal_idx = find(strcmp(unique_animals, highGamma_animal{i})); % Find the index of the animal for shape
    color_c = 'black';
    color_w = 'black';
    if highGamma_sig_c(i) == 1
        color_c = 'green';
    end
    if highGamma_sig_w(i) == 1
        color_w = 'green';
    end
    % Plot individual points with filled color
    scatter(x_categories(1), highGamma_clustermax_c(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_c, 'MarkerFaceColor', color_c, 'MarkerFaceAlpha', 0.35);
    scatter(x_categories(2), highGamma_clustermax_w(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_w, 'MarkerFaceColor', color_w, 'MarkerFaceAlpha', 0.35);
    
    % Draw lines connecting the points
    line([x_categories(1), x_categories(2)], [highGamma_clustermax_c(i), highGamma_clustermax_w(i)], 'Color', 'black', 'LineWidth', line_width);
end
% Customize the plot
xlim([0.5, 2.5]); % Adjust as needed for better spacing
set(gca, 'xtick', x_categories, 'xticklabel', {'CORRECT', 'WRONG'});
ylabel('Z Max Cluster Value');
title('highGamma');
hold off;

subplot(2,3,6)
hold on;
% Overlay individual data points
for i = 1:length(all_clustermax_c)
    animal_idx = find(strcmp(unique_animals, all_animal{i})); % Find the index of the animal for shape
    color_c = 'black';
    color_w = 'black';
    if all_sig_c(i) == 1
        color_c = 'green';
    end
    if all_sig_w(i) == 1
        color_w = 'green';
    end
    % Plot individual points with filled color
    scatter(x_categories(1), all_clustermax_c(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_c, 'MarkerFaceColor', color_c, 'MarkerFaceAlpha', 0.35);
    scatter(x_categories(2), all_clustermax_w(i), animal_shapes{animal_idx},  'filled', 'MarkerEdgeColor', color_w, 'MarkerFaceColor', color_w, 'MarkerFaceAlpha', 0.35);
    
    % Draw lines connecting the points
    line([x_categories(1), x_categories(2)], [all_clustermax_c(i), all_clustermax_w(i)], 'Color', 'black', 'LineWidth', line_width);
end
% Customize the plot
xlim([0.5, 2.5]); % Adjust as needed for better spacing
set(gca, 'xtick', x_categories, 'xticklabel', {'CORRECT', 'WRONG'});
ylabel('Z Max Cluster Value');
title('all');
hold off;