
% For a specified frequency band, loops through all sessions and counts the 
% number of times the p value for coherence or granger passes (< 0.05) for 
% all channel pairs. Heatmap is generated to visualize which pairs are 
% passing more.

%% Define data

Frequency_Band  = 'theta';               % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'Granger';            % 'Granger' or 'Coherence'
animals         = {'MrM'};             % 'MrCassius' and/or 'MrM'

rootdir = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis';
sessions = dir(fullfile(rootdir, '19*'));

if strcmp(Statistic, 'Coherence')
    ChanPair_Array_Coh    = zeros(20,20);   % used to store passing p-value count for coherence
elseif strcmp(Statistic, 'Granger')
    ChanPair_Array_PFC_AC = zeros(20,20);   % used to store passing p-value count for granger (bidirectional)
    ChanPair_Array_AC_PFC = zeros(20,20);
end

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        Epoch  = 'testToneOnset';

        if exist(fullfile(rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
            if strcmp(Statistic, 'Coherence')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), ['*' Frequency_Band '*Coh_data.mat']));
            elseif strcmp(Statistic, 'Granger')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), ['*' Frequency_Band '*Granger_data.mat']));
            end
        else
            continue;
        end

        if ~isempty(files)
            fprintf('Processing session %s: ', num2str(RecDate));
        end
        p_count = 0;
        for k = 1:length(files) % runs if there are any processed channel pairs in given frequency band
        
            file_name = files(k).name;

            % pull channel pair info from file name
            tokens = regexp(file_name, '(.*?)_', 'tokens');
            if isempty(tokens)
                continue;
            end

            Send_Cort   = tokens{4}{1}; 
            Send_Num    = tokens{5}{1};
            Rec_Cort    = tokens{6}{1};
            Rec_Num     = tokens{7}{1};

            Send_Chan   = [Send_Cort '_' Send_Num];
            Rec_Chan    = [Rec_Cort '_' Rec_Num];

            % pull all present data files for channel pair
            if strcmp(Statistic, 'Coherence')
                channel_files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), [extractBefore(file_name, 'Coh') '*']));
            elseif strcmp(Statistic, 'Granger')
                channel_files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), [extractBefore(file_name, 'Granger') '*']));
            end

            % check if Maris statistical test has been run on the channel pair yet (pair can have Coh_data/Granger_data file without
            % having statistical test results stored in that file since RunMaris_v3.m appends data from statistical test)
            if length(channel_files) ~= 15      % number of files for a channel pair is 15 after Maris runs
                continue;
            end
            
            % load p value of frequency band
            datadir = fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), file_name);

            if strcmp(Statistic, 'Coherence')
                % load in p-value
                load(datadir,'pvalue_c');
                pvalue_c = pvalue_c.(Frequency_Band);
    
                % increase passing p-value count at appropriate index if p-value < 0.05
                pfc_idx    = str2double(extractAfter(Send_Num, 'ch')) - 2;    % channels start at ch03
                ac_idx     = str2double(extractAfter(Rec_Num, 'ch')) - 2;     % adjust index to fit all pairs in 20x20 array
    
                if pvalue_c <= 0.05
                    ChanPair_Array_Coh(pfc_idx, ac_idx) = ChanPair_Array_Coh(pfc_idx, ac_idx) + 1;
                    p_count = p_count + 1;
                end

            elseif strcmp(Statistic, 'Granger')
                % load in p-value
                load(datadir,'pvalue_PFC_AC','pvalue_AC_PFC');
                pvalue_PFC_AC = pvalue_PFC_AC.(Frequency_Band);
                pvalue_AC_PFC = pvalue_AC_PFC.(Frequency_Band);
    
                % increase passing p-value count at appropriate index if p-value < 0.05
                pfc_idx    = str2double(extractAfter(Send_Num, 'ch')) - 2;    % channels start at ch03
                ac_idx     = str2double(extractAfter(Rec_Num, 'ch')) - 2;     % adjust index to fit all pairs in 20x20 array
    
                if pvalue_PFC_AC <= 0.05
                    ChanPair_Array_PFC_AC(pfc_idx, ac_idx) = ChanPair_Array_PFC_AC(pfc_idx, ac_idx) + 1;
                    p_count = p_count + 1;
                end
                if pvalue_AC_PFC <= 0.05
                    ChanPair_Array_AC_PFC(ac_idx, pfc_idx) = ChanPair_Array_AC_PFC(ac_idx, pfc_idx) + 1;
                    p_count = p_count + 1;
                end
            end
        end
        if ~isempty(files)
            disp([num2str(p_count) ' passing p-values found.']);
        end
    end
end

%% generate grid and diplot for coherence values 

if strcmp(Statistic, 'Coherence')

        figure;
        
        % configure grid
        s = [20 20]; % [y x]
        xrange = [3 22]; % imagesc only needs the endpoints
        yrange = [3 22];
        dx = diff(xrange)/(s(2)-1);
        dy = diff(yrange)/(s(1)-1);
        xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
        yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);
        
        % plot heatmap
        imagesc(xrange,yrange,ChanPair_Array_Coh); hold on
        colorbar
        
        % adjust coloring
        hm = mesh(xg,yg,zeros(s+1));
        hm.FaceColor = 'none';
        hm.EdgeColor = 'k';
        
        title(['Coherence pvalues in ' Frequency_Band]);
        xlabel('AC');
        ylabel('PFC');
       
        % plot digraph (coherence)
        
        % PFC to AC
        
        % create adjacency matrix (directional array representing edges/arrows in the resulting digraph)
        adj_matrix = zeros(40);  
        adj_matrix(1:20, 21:40) = ChanPair_Array_Coh;    % fill adjacency matrix at proper indices (1:20 are PFC channels; 21:40 are AC channels)
        adj_matrix = adj_matrix + eye(40);   % add self-loops to make sure all nodes appear if they have no edges/arrows
        
        % find the source nodes, target nodes, and weights ('find' is a Matlab function)
        [s, t, weights] = find(adj_matrix);
        
        % create digraph
        G = digraph(s, t, weights, 'omitselfloops');    % 'omitselfloops' to have self-looping edges not appear
        
        % plot digraph
        figure;
        labels = [arrayfun(@(x) sprintf('PFC %d', x), 3:22, 'UniformOutput', false), ...    % labels nodes PFC 3, PFC 4,..., PFC 22
                  arrayfun(@(x) sprintf('AC %d', x), 3:22, 'UniformOutput', false)];        % labels nodes AC 3, AC 4,..., AC 22
        
        h = plot(G, 'NodeLabel', labels, 'Layout', 'layered');
        
        % adjust edge thickness/color based on weights
        h.LineWidth = G.Edges.Weight; % use edge weights as line width
        for e = 1:height(G.Edges)
            if G.Edges.Weight(e) > 1
                highlight(h, 'Edges', e, 'EdgeColor', 'r');     % 'highlight' function changes edge color
            end
        end
        
        % color nodes and edges
        h.NodeColor = 'k'; % set node color to black
        h.MarkerSize = 7; % set node size
        
        % set node positions
        x_senders = ones(1, 20) * 1; % x = 1 for senders
        y_senders = linspace(20, 1, 20); % evenly spaced y positions for senders
        x_receivers = ones(1, 20) * 2; % x = 2 for receivers
        y_receivers = linspace(20, 1, 20); % evenly spaced y positions for receivers
        x = [x_senders, x_receivers];
        y = [y_senders, y_receivers];
        h.XData = x;
        h.YData = y;
        
        title(['Coherence betweeen PFC and AC ' Frequency_Band, ' ', Animal] );

end

%% generate grid and diplot for granger values in both directions  


if strcmp(Statistic, 'Granger')

        % PFC to AC
        
        figure;
        
        % configure grid
        s = [20 20]; % [y x]
        xrange = [3 22]; % imagesc only needs the endpoints
        yrange = [3 22];
        dx = diff(xrange)/(s(2)-1);
        dy = diff(yrange)/(s(1)-1);
        xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
        yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);
        
        % plot heatmap
        imagesc(xrange,yrange,ChanPair_Array_PFC_AC); hold on
        colorbar
        
        % adjust coloring
        hm = mesh(xg,yg,zeros(s+1));
        hm.FaceColor = 'none';
        hm.EdgeColor = 'k';
        
        title(['PFC to AC Granger pvalues in ' Frequency_Band]);
        xlabel('AC');
        ylabel('PFC');
        
        % AC to PFC
        
        figure;
        
        % configure grid
        s = [20 20]; % [y x]
        xrange = [3 22]; % imagesc only needs the endpoints
        yrange = [3 22];
        a = rand(s);
        dx = diff(xrange)/(s(2)-1);
        dy = diff(yrange)/(s(1)-1);
        xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
        yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);
        
        % plot heatmap
        hi = imagesc(xrange,yrange,ChanPair_Array_AC_PFC); hold on
        colorbar
        
        % adjust coloring
        hm = mesh(xg,yg,zeros(s+1));
        hm.FaceColor = 'none';
        hm.EdgeColor = 'k';
        
        title(['AC to PFC Granger pvalues in ' Frequency_Band]);
        xlabel('PFC');
        ylabel('AC');
        
        % plot digraph (granger)
        
        % PFC to AC
        
        % create adjacency matrix (directional array representing edges/arrows in the resulting digraph)
        adj_matrix = zeros(40);  
        adj_matrix(1:20, 21:40) = ChanPair_Array_PFC_AC;    % fill adjacency matrix at proper indices (1:20 are PFC channels; 21:40 are AC channels)
        adj_matrix = adj_matrix + eye(40);   % add self-loops to make sure all nodes appear if they have no edges/arrows
        
        % find the source nodes, target nodes, and weights ('find' is a Matlab function)
        [s, t, weights] = find(adj_matrix);
        
        % create digraph
        G = digraph(s, t, weights, 'omitselfloops');    % 'omitselfloops' to have self-looping edges not appear
        
        % plot digraph
        figure;
        labels = [arrayfun(@(x) sprintf('PFC %d', x), 3:22, 'UniformOutput', false), ...    % labels nodes PFC 3, PFC 4,..., PFC 22
                  arrayfun(@(x) sprintf('AC %d', x), 3:22, 'UniformOutput', false)];        % labels nodes AC 3, AC 4,..., AC 22
        
        h = plot(G, 'NodeLabel', labels, 'Layout', 'layered');
        
        % adjust edge thickness/color based on weights
        h.LineWidth = G.Edges.Weight; % use edge weights as line width
        for e = 1:height(G.Edges)
            if G.Edges.Weight(e) > 1
                highlight(h, 'Edges', e, 'EdgeColor', 'r');     % 'highlight' function changes edge color
            end
        end
        
        % color nodes and edges
        h.NodeColor = 'k'; % set node color to black
        h.MarkerSize = 7; % set node size
        
        % set node positions
        x_senders = ones(1, 20) * 1; % x = 1 for senders
        y_senders = linspace(20, 1, 20); % evenly spaced y positions for senders
        x_receivers = ones(1, 20) * 2; % x = 2 for receivers
        y_receivers = linspace(20, 1, 20); % evenly spaced y positions for receivers
        x = [x_senders, x_receivers];
        y = [y_senders, y_receivers];
        h.XData = x;
        h.YData = y;
        
        title(['PFC to AC Granger in ' Frequency_Band]);
        
        
        % AC to PFC
        
        adj_matrix = zeros(40);  % create adjacency matrix
        adj_matrix(1:20, 21:40) = ChanPair_Array_AC_PFC;
        adj_matrix = adj_matrix + eye(40);   % add self-loops to make sure all nodes appear even without an edge
        
        % Find the source nodes, target nodes, and weights
        [s, t, weights] = find(adj_matrix);
        
        % Create the digraph
        G = digraph(s, t, weights, 'omitselfloops');
        
        % Plot the digraph
        figure;
        labels = [arrayfun(@(x) sprintf('AC %d', x), 3:22, 'UniformOutput', false), ...
                  arrayfun(@(x) sprintf('PFC %d', x), 3:22, 'UniformOutput', false)];
        
        h = plot(G, 'NodeLabel', labels, 'Layout', 'layered');
        
        % Adjust edge thickness/color based on weights
        h.LineWidth = G.Edges.Weight; % Use edge weights as line width
        for e = 1:height(G.Edges)
            if G.Edges.Weight(e) > 1
                highlight(h, 'Edges', e, 'EdgeColor', 'r');
            end
        end
        
        % Color nodes and edges
        h.NodeColor = 'k'; % Set node color to black
        % h.EdgeColor = 'b';
        h.MarkerSize = 7; % Set node size
        
        % Set the node positions
        x_senders = ones(1, 20) * 1; % x = 1 for senders
        y_senders = linspace(20, 1, 20); % evenly spaced y positions for senders
        x_receivers = ones(1, 20) * 2; % x = 2 for receivers
        y_receivers = linspace(20, 1, 20); % evenly spaced y positions for receivers
        x = [x_senders, x_receivers];
        y = [y_senders, y_receivers];
        h.XData = x;
        h.YData = y;
        
        title(['AC to PFC Granger in ' Frequency_Band]);

end
