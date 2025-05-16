%% Make scatter plots of clustermax values of the correct trials 


% Extract relevant variables
session = outputTable.session;
clustermax_c = outputTable.clustermax_c;
animal = outputTable.animal;
sig_c = outputTable.sig_c;
band = outputTable.band;

% Create theta scatter plot 

figure;

sgtitle('CORRECT')

subplot(2,3,1)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_c(i);
    currentAnimal = animal{i};
    isSigC = sig_c(i);
    currentBand = band{i};
    
    % Plot only 'theta' data points
    if strcmp(currentBand, 'theta')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigC
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('theta');
hold off;

% alpha
subplot(2,3,2)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_c(i);
    currentAnimal = animal{i};
    isSigC = sig_c(i);
    currentBand = band{i};
    
    % Plot only 'alpha' data points
    if strcmp(currentBand, 'alpha')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigC
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('alpha');
hold off;

%beta
subplot(2,3,3)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_c(i);
    currentAnimal = animal{i};
    isSigC = sig_c(i);
    currentBand = band{i};
    
    % Plot only 'beta' data points
    if strcmp(currentBand, 'beta')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigC
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('beta');
hold off;


% gamma
subplot(2,3,4)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_c(i);
    currentAnimal = animal{i};
    isSigC = sig_c(i);
    currentBand = band{i};
    
    % Plot only 'gamma' data points
    if strcmp(currentBand, 'gamma')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigC
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('gamma');
hold off;

% highGamma
subplot(2,3,5)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_c(i);
    currentAnimal = animal{i};
    isSigC = sig_c(i);
    currentBand = band{i};
    
    % Plot only 'theta' data points
    if strcmp(currentBand, 'highGamma')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigC
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('highGamma');
hold off;

% all
subplot(2,3,6)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_c(i);
    currentAnimal = animal{i};
    isSigC = sig_c(i);
    currentBand = band{i};
    
    % Plot only 'theta' data points
    if strcmp(currentBand, 'all')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigC
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('all');
hold off;



%% Make scatter plots of clustermax values of the wrong trials 


% Extract relevant variables
session = outputTable.session;
clustermax_w = outputTable.clustermax_w;
animal = outputTable.animal;
sig_w = outputTable.sig_w;
band = outputTable.band;

% Create theta scatter plot 

figure;

sgtitle('WRONG')

subplot(2,3,1)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_w(i);
    currentAnimal = animal{i};
    isSigW = sig_w(i);
    currentBand = band{i};
    
    % Plot only 'theta' data points
    if strcmp(currentBand, 'theta')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigW
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('theta');
hold off;

% alpha
subplot(2,3,2)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_w(i);
    currentAnimal = animal{i};
    isSigW = sig_w(i);
    currentBand = band{i};
    
    % Plot only 'alpha' data points
    if strcmp(currentBand, 'alpha')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigW
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('alpha');
hold off;

%beta
subplot(2,3,3)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_w(i);
    currentAnimal = animal{i};
    isSigW = sig_w(i);
    currentBand = band{i};
    
    % Plot only 'beta' data points
    if strcmp(currentBand, 'beta')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigW
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('beta');
hold off;


% gamma
subplot(2,3,4)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_w(i);
    currentAnimal = animal{i};
    isSigW = sig_w(i);
    currentBand = band{i};
    
    % Plot only 'gamma' data points
    if strcmp(currentBand, 'gamma')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigW
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('gamma');
hold off;

% highGamma
subplot(2,3,5)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_w(i);
    currentAnimal = animal{i};
    isSigW= sig_w(i);
    currentBand = band{i};
    
    % Plot only 'theta' data points
    if strcmp(currentBand, 'highGamma')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigW
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('highGamma');
hold off;

% all

subplot(2,3,6)
hold on;
% Set marker shapes for each animal
animalMarkers = containers.Map({'MrCassius', 'MrM'}, {'o', '^'});

% Iterate through each data point
for i = 1:length(session)
    x = session(i);
    y = clustermax_w(i);
    currentAnimal = animal{i};
    isSigW = sig_w(i);
    currentBand = band{i};
    
    % Plot only all data points
    if strcmp(currentBand, 'all')
        % Get the marker shape for the current animal
        marker = animalMarkers(currentAnimal);
        
        % Set marker face color based on sig_c
        if isSigW
            faceColor = 'green';
            edgeColor = 'green';
        else
            faceColor = 'black';
            edgeColor = 'black';
        end
        
        % Plot the data point
        scatter(x, y, 'Marker', marker, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', edgeColor);
    end
end

% Customize plot
xlabel('Session');
ylabel('Clustermax Correct');
title('all');
hold off;