
% For a specified frequency band, loops through all sessions and counts the 
% number of times the p value for coherence or granger passes (< 0.05) for 
% all channel pairs. Heatmap is generated to visualize the spatial
% distribution of the channel pairs that significantly modulated by context 

%% Define data

Frequency_Band  = 'theta';               % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'Coherence';             % 'Granger' or 'Coherence'
animals         = {'MrCassius'};         % 'MrCassius' and/or 'MrM'

rootdir  = 'E:\2024_09_27_TestTone_Correct';
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
        Epoch  = 'testToneOnset';   % Epoch  = 'testToneOnset'; 'preCueOnset'
     
        if exist(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct'), 'dir')
            if strcmp(Statistic, 'Coherence')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', ['*' Frequency_Band '*Coh_data.mat']));
            elseif strcmp(Statistic, 'Granger')
            files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', ['*' Frequency_Band '*Granger_data.mat']));
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

            tokens = regexp(file_name, '(.*?)_(\d+)_testToneOnset_(.*?)_ch(\d+)_AC_ch(\d+)', 'tokens');
            if isempty(tokens)
                continue;
            end

            Send_Cort = tokens{1}{3}; % Updated cortical region parsing
            Send_Num  = ['ch' tokens{1}{4}]; % Updated channel number parsing
            Rec_Cort  = 'AC'; % Fixed cortical region as 'AC' from the pattern
            Rec_Num   = ['ch' tokens{1}{5}];

            % pull all present data files for channel pair
            if strcmp(Statistic, 'Coherence')
            channel_files = dir(fullfile(rootdir, RecDate, Animal, 'testTone', 'Correct', [extractBefore(file_name, 'Coh') '*']));
            elseif strcmp(Statistic, 'Granger')
                channel_files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), [extractBefore(file_name, 'Granger') '*']));
            end

            % check if Maris statistical test has been run on the channel pair yet (pair can have Coh_data/Granger_data file without
            % % having statistical test results stored in that file since RunMaris_v3.m appends data from statistical test)
            % if length(channel_files) ~= 9      % number of files for a channel pair is 15 after Maris runs
            %     continue;
            % end

            % Load p value of frequency band
            datadir = fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), 'Correct', file_name);


            if strcmp(Statistic, 'Coherence')
                % load in p-value
                load(datadir,'pvalue_c');
                pvalue_c = pvalue_c.(Frequency_Band);
    
                % increase passing p-value count at appropriate index if p-value < 0.05
                pfc_idx = str2double(extractAfter(Send_Num, 'ch')) - 2; % PFC index
                ac_idx  = str2double(extractAfter(Rec_Num, 'ch')) - 2;  % AC index


    
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

        figure

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
        clim([0 8]);   % Set minimum color value to 0 and maximum to 10
        colorbar
        cb = colorbar;
        ylabel(cb, '# of Sig.Interactions', 'FontSize', 12)
        cb.Label.Rotation = 270;

        % adjust coloring
        hm = mesh(xg,yg,zeros(s+1));
        hm.FaceColor = 'none';
        hm.EdgeColor = 'k';

        title([Statistic, '-', Frequency_Band, '-', Animal],...
            'Fontsize', 12);
        xlabel('AC Channels', 'FontSize', 12);
        ylabel('PFC Channels', 'FontSize', 12);

        ax = gca; % Get current axis
        ax.FontSize = 12; % Set font size for tick marks

        % configure grid
% s = [20 20]; % [y x]
% xrange = [3 22]; % imagesc only needs the endpoints
% yrange = [3 22];
% dx = diff(xrange)/(s(2)-1);
% dy = diff(yrange)/(s(1)-1);
% xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
% yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);
% 
% % plot heatmap
% figure; % Create new figure
% imagesc(xrange, yrange, ChanPair_Array_Coh); hold on
% clim([0 8]);   % Set minimum color value to 0 and maximum to 10
% colorbar;
% cb = colorbar;
% ylabel(cb, '# of Sig.Interactions', 'FontSize', 24);
% cb.Label.Rotation = 270;
% 
% % adjust coloring
% hm = mesh(xg, yg, zeros(s+1));
% hm.FaceColor = 'none';
% hm.EdgeColor = 'k';
% 
% % Labels and axis settings
% xlabel('AC Channels', 'FontSize', 32);
% ylabel('PFC Channels', 'FontSize', 32);
% ax = gca; % Get current axis
% ax.FontSize = 20; % Set font size for tick marks
% 
% % Set figure size to 7.65 inches wide and 6 inches tall
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 7.65 6]); % [left, bottom, width, height]
% 
% % Define save directory and file name
% saveDir = 'C:\Users\Corey Roach\Desktop\2024_SFN_Poster'; % Replace with your desired path
% saveFileName = 'heatmap_highres.png'; % File name
% 
% % Full file path
% fullFilePath = fullfile(saveDir, saveFileName);
% 
% % Save the figure as a high-resolution image (e.g., 300 dpi)
% print(gcf, fullFilePath, '-dpng', '-r300'); % Save as PNG at 300 dpi


       
        %
%         % plot digraph (coherence)
% 
%         % PFC to AC
% 
%         % create adjacency matrix (directional array representing edges/arrows in the resulting digraph)
%         adj_matrix = zeros(40);  
%         adj_matrix(1:20, 21:40) = ChanPair_Array_Coh;    % fill adjacency matrix at proper indices (1:20 are PFC channels; 21:40 are AC channels)
%         adj_matrix = adj_matrix + eye(40);   % add self-loops to make sure all nodes appear if they have no edges/arrows
% 
%         % find the source nodes, target nodes, and weights ('find' is a Matlab function)
%         [s, t, weights] = find(adj_matrix);
% 
%         % create digraph
%         G = digraph(s, t, weights, 'omitselfloops');    % 'omitselfloops' to have self-looping edges not appear
% 
%         % plot digraph
%         figure;
%         labels = [arrayfun(@(x) sprintf('PFC %d', x), 3:22, 'UniformOutput', false), ...    % labels nodes PFC 3, PFC 4,..., PFC 22
%                   arrayfun(@(x) sprintf('AC %d', x), 3:22, 'UniformOutput', false)];        % labels nodes AC 3, AC 4,..., AC 22
% 
%         h = plot(G, 'NodeLabel', labels, 'Layout', 'layered','ArrowSize', 0);
% 
%         % adjust edge thickness/color based on weights
%         h.LineWidth = G.Edges.Weight; % use edge weights as line width
%         for e = 1:height(G.Edges)
%             if G.Edges.Weight(e) > 1
%                 highlight(h, 'Edges', e, 'EdgeColor', 'r');     % 'highlight' function changes edge color
%             end
%         end
% 
%         % color nodes and edges
%         h.NodeColor = 'k'; % set node color to black
%         h.MarkerSize = 7; % set node size
% 
%         % set node positions
%         x_senders = ones(1, 20) * 1; % x = 1 for senders
%         y_senders = linspace(20, 1, 20); % evenly spaced y positions for senders
%         x_receivers = ones(1, 20) * 2; % x = 2 for receivers
%         y_receivers = linspace(20, 1, 20); % evenly spaced y positions for receivers
%         x = [x_senders, x_receivers];
%         y = [y_senders, y_receivers];
%         h.XData = x;
%         h.YData = y;
% 
%         %title(['Distribution of Channel Pairs Significantly Modulated by Context', ':  ', Statistic, '-', Frequency_Band, '-', Animal]);
% % 

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
        
        %title(['Distribution of Channel Pairs Significantly Modulated by Context', ':  ', Statistic, '-', Frequency_Band, '-', Animal]);
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
        
        %title(['Distribution of Channel Pairs Significantly Modulated by Context', ':  ', Statistic, '-', Frequency_Band, '-', Animal]);
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
        
       %title(['Distribution of Channel Pairs Significantly Modulated by Context', ':  ', Statistic, '-', Frequency_Band, '-', Animal]);
    
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
        
        title(['Distribution of Channel Pairs Significantly Modulated by Context', ':  ', Statistic, '-', Frequency_Band, '-', Animal]);
end
