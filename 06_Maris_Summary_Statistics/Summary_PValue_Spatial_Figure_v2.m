
% For a specified frequency band, loops through all sessions and counts the 
% number of times the p value for coherence or granger passes (< 0.05) for 
% all channel pairs. Heatmap is generated to visualize the spatial
% distribution of the channel pairs that significantly modulated by context 

%% Define data

Frequency_Band  = 'alpha';         % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'Coherence';     % 'Granger' or 'Coherence'
animals         = {'MrM'};         % 'MrCassius' and/or 'MrM'

rootdir  = '\\Kilosort\d\Top_Down_Coherence_Project\00_DATA\zz_METADATA\2024_10_15_PreCue_Correct\00_Core';
sessions = dir(fullfile(rootdir, '19*'));

if strcmp(Statistic, 'Coherence')
    ChanPair_Array_Coh_sig    = zeros(20,20);   % used to store passing p-value count for coherence
    ChanPair_Array_Coh_total  = zeros(20, 20);

end

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        Epoch  = 'preCueOnset';   % Epoch  = 'testToneOnset'; 'preCueOnset'

        if exist(fullfile(rootdir, RecDate, Animal), 'dir') % if the session has the animal (one animal per session)
            if strcmp(Statistic, 'Coherence')
                files = dir(fullfile(rootdir, RecDate, Animal, extractBefore(Epoch, 'Onset'), ['*' Frequency_Band '*Coh_data.mat']));
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

                % load in p-value
                load(datadir,'pvalue_c');
                pvalue_c = pvalue_c.(Frequency_Band);
    
                % increase passing p-value count at appropriate index if p-value < 0.05
                pfc_idx    = str2double(extractAfter(Send_Num, 'ch')) - 2;    % channels start at ch03
                ac_idx     = str2double(extractAfter(Rec_Num, 'ch')) - 2;     % adjust index to fit all pairs in 20x20 array
    
                % Increment total appearances for this pair
                ChanPair_Array_Coh_total(pfc_idx, ac_idx) = ChanPair_Array_Coh_total(pfc_idx, ac_idx) + 1;

                if pvalue_c <= 0.05
                    ChanPair_Array_Coh_sig(pfc_idx, ac_idx) = ChanPair_Array_Coh_sig(pfc_idx, ac_idx) + 1;
                    p_count = p_count + 1;
                end            
        end

        if ~isempty(files)
            disp([num2str(p_count) ' passing p-values found.']);
        end

        
    end
end

%% generate grid and diplot for coherence values 

% Avoid divide-by-zero warnings

ChanPair_Array_Coh_prop = ChanPair_Array_Coh_sig ./ max(ChanPair_Array_Coh_total, 1);

figure (1)

% configure grid
s = [20 20]; % [y x]
xrange = [3 22]; % imagesc only needs the endpoints
yrange = [3 22];
dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);

% plot heatmap
imagesc(xrange,yrange,ChanPair_Array_Coh_prop); hold on
clim([0 1]);   % Set minimum color value to 0 and maximum to 10
colorbar
cb = colorbar;
ylabel(cb, 'Proportion Significant', 'FontSize', 12)
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

%% Plot Counts

figure (2)

% configure grid
s = [20 20]; % [y x]
xrange = [3 22]; % imagesc only needs the endpoints
yrange = [3 22];
dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);

% plot heatmap
imagesc(xrange,yrange, ChanPair_Array_Coh_sig); hold on
clim([0 8]);   % Set minimum color value to 0 and maximum to 10
colorbar
cb = colorbar;
ylabel(cb, 'Number of Significant Interactions', 'FontSize', 12)
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




