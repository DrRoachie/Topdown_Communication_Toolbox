
% For a specified frequency band, loops through all sessions and counts the 
% number of times the p value for coherence (< 0.05) for 
% all channel pairs. Heatmap is generated to visualize the spatial
% distribution of the channel pairs that significantly modulated by context 

%% Define data

Frequency_Band  = 'theta';         % 'theta', 'alpha', 'beta', 'gamma', 'highGamma'
Statistic       = 'Coherence';     % 'Granger' or 'Coherence'
animals         = {'MrCassius'};         % 'MrCassius' and/or 'MrM'

rootdir  = '\\Kilosort\u\Top_Down_Coherence_Project\00_DATA\zz_METADATA\2024_09_27_TestTone_Correct\00_data';
sessions = dir(fullfile(rootdir, '19*'));

if strcmp(Statistic, 'Coherence')
    ChanPair_Array_Coh_sig    = zeros(20,20);   % used to store passing p-value count for coherence
    ChanPair_Array_Coh_total  = zeros(20, 20);

end

for i = 1:length(sessions)

    RecDate = sessions(i).name;

    for j = 1:length(animals)

        Animal = animals{j};
        Epoch  = 'testToneOnset';   % Epoch  = 'testToneOnset'; 'preCueOnset'

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

% %% Proportion Plot [DO NOT DELETE]
% 
% % generate grid for proportion significant out of the total
% 
% % Avoid divide-by-zero warnings
% 
% figure(1)
% 
% ChanPair_Array_Coh_prop = ChanPair_Array_Coh_sig ./ max(ChanPair_Array_Coh_total, 1);
% 
% 
% % configure grid
% s = [20 20]; % [y x]
% xrange = [3 22]; % imagesc only needs the endpoints
% yrange = [3 22];
% dx = diff(xrange)/(s(2)-1);
% dy = diff(yrange)/(s(1)-1);
% xg = linspace(xrange(1)-dx/2,xrange(2)+dx/2,s(2)+1);
% yg = linspace(yrange(1)-dy/2,yrange(2)+dy/2,s(1)+1);
% 
% % plot heatmap
% imagesc(xrange,yrange,ChanPair_Array_Coh_prop); hold on
% clim([0 1]);   % Set minimum color value to 0 and maximum to 10
% colorbar
% cb = colorbar;
% ylabel(cb, 'Proportion Significant', 'FontSize', 12)
% cb.Label.Rotation = 270;
% 
% % adjust coloring
% hm = mesh(xg,yg,zeros(s+1));
% hm.FaceColor = 'none';
% hm.EdgeColor = 'k';
% 
% title([Statistic, '-', Frequency_Band, '-', Animal],...
%     'Fontsize', 12);
% xlabel('AC Channels', 'FontSize', 12);
% ylabel('PFC Channels', 'FontSize', 12);
% 
% ax = gca; % Get current axis
% ax.FontSize = 12; % Set font size for tick marks

%% Plot Counts with zero-valued cells in dark navy (aligned overlay)
figure; clf

% --- Global formatting (Arial, 12 pt, bold) ---
set(groot, 'defaultAxesFontName', 'Arial', ...
           'defaultTextFontName', 'Arial', ...
           'defaultAxesFontSize', 16, ...
           'defaultTextFontSize', 16);

% --- config ---
s = [20 20]; % [y x]
xrange = [3 22];
yrange = [3 22];

% Cell geometry
dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xc = linspace(xrange(1), xrange(2), s(2));
yc = linspace(yrange(1), yrange(2), s(1));
xedges = linspace(xrange(1)-dx/2, xrange(2)+dx/2, s(2)+1);
yedges = linspace(yrange(1)-dy/2, yrange(2)+dy/2, s(1)+1);

% --- data prep ---
plotData  = ChanPair_Array_Coh_sig;
zeroMask  = (plotData == 0) | isnan(plotData);

% Default auto-scaling
posVals = plotData(~zeroMask);
minPos  = iff(~isempty(posVals), min(posVals), 0);
maxVal  = max(plotData(:), [], 'omitnan');

% --- user-defined color range (set manually as needed) ---
% Example: fix all plots to 0–8
cRange = [0 6];      
% Or use auto: cRange = [minPos maxVal];

% --- base heatmap ---
hIm = imagesc(xrange, yrange, plotData); hold on;
set(gca,'XDir','normal','YDir','reverse');
clim(cRange);
hIm.AlphaData = double(~zeroMask); % transparent where 0/NaN

% --- overlay dark navy tiles ---
darkNavy = [0 0 0.35];
for r = 1:s(1)
    ylow = yedges(r);
    yhei = yedges(r+1) - yedges(r);
    for c = 1:s(2)
        if zeroMask(r,c)
            xlow = xedges(c);
            xwid = xedges(c+1) - xedges(c);
            rectangle('Position', [xlow, ylow, xwid, yhei], ...
                      'FaceColor', darkNavy, 'EdgeColor', 'none');
        end
    end
end

% --- grid overlay ---
hm = mesh(xedges, yedges, zeros(s+1));
hm.FaceColor = 'none';
hm.EdgeColor = 'k';
hm.LineWidth = 0.5;
hm.EdgeAlpha = 0.8;

% --- ticks ---
xticks(xc); xticklabels(1:s(2));
yticks(linspace(yrange(1), yrange(2), s(1)));
yticklabels(1:s(1));  % top row = 1

% --- colorbar & labels ---
cb = colorbar;
ylabel(cb, 'Number of Significant Interactions', ...
       'FontSize', 16, 'FontName','Arial', 'FontWeight','bold');
cb.Label.Rotation = 270;

xlabel('AC Channels', 'FontSize', 12, 'FontName','Arial', 'FontWeight','bold');
ylabel('PFC Channels', 'FontSize', 12, 'FontName','Arial', 'FontWeight','bold');

ax = gca;
ax.FontSize = 16; ax.FontName = 'Arial'; ax.FontWeight = 'bold';
box on; axis square;

% --- export ---
set(gcf,'Units','inches','Position',[0 0 11 8]);
exportgraphics(gcf, 'Counts_11x8in_600dpi_bold.tif', 'Resolution', 600);


% %% [DO NOT DELETE]
% 
% % plot digraph (coherence)
% 
% figure(3)
%         % PFC to AC
% 
%         % Mean of all elements
%         mu_all = mean(ChanPair_Array_Coh_sig(:));
%         std_all = std(ChanPair_Array_Coh_sig(:));
%         upper_2sd = mu_all + 2*std_all;      % mean + 2*SD
% 
%         % create adjacency matrix (directional array representing edges/arrows in the resulting digraph)
%         adj_matrix = zeros(40);  
%         adj_matrix(1:20, 21:40) = ChanPair_Array_Coh_sig;    % fill adjacency matrix at proper indices (1:20 are PFC channels; 21:40 are AC channels)
%         adj_matrix = adj_matrix + eye(40);   % add self-loops to make sure all nodes appear if they have no edges/arrows
% 
%         % find the source nodes, target nodes, and weights ('find' is a Matlab function)
%         [s, t, weights] = find(adj_matrix);
% 
%         % create digraph
%         G = digraph(s, t, weights, 'omitselfloops');    % 'omitselfloops' to have self-looping edges not appear
% 
%         % plot digraph
% 
%         labels = [arrayfun(@(x) sprintf('PFC %d', x), 1:20, 'UniformOutput', false), ...    % labels nodes PFC 3, PFC 4,..., PFC 22
%                   arrayfun(@(x) sprintf('AC %d', x), 1:20, 'UniformOutput', false)];        % labels nodes AC 3, AC 4,..., AC 22
% 
%         h = plot(G, 'NodeLabel', labels, 'Layout', 'layered','ArrowSize', 0);
% 
%         % set all edges to same baseline style (blue, thin)
%         h.LineWidth = 1;
% 
%         % now selectively make the strong edges bold + red
%         for e = 1:height(G.Edges)
%             if G.Edges.Weight(e) > upper_2sd
%                 highlight(h, 'Edges', e, ...
%                     'EdgeColor', 'r', ...
%                     'LineWidth', 4);   % fixed bolder width
%             end
%         end
% 
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


%% plots the number of total number of sig interactions per session

% figure(4)
% if strcmp(Animal, 'MrCassius') && strcmp(Frequency_Band, 'theta')
%     Cassius_Sig_per_Session = [29,5,4,1,13,12,55,7,6,9,4,68,71,156,2,15,13,82,16,21];
%     nSessions = length(Cassius_Sig_per_Session);
% 
%     plot(1:nSessions, Cassius_Sig_per_Session, 'LineWidth', 3);
%     xlabel('Session');
%     ylabel('Number of Modulated Channel Pairs');
%     ylim([0 200]);
% 
%     % Add ticks at each session
%     xticks(1:nSessions);
%     xlim([0 nSessions+1]);  % keep some spacing on both ends
% end
% 
% 
% if strcmp(Animal, 'MrM') && strcmp(Frequency_Band, 'theta')
%     Cassius_Sig_per_Session = [65,108,1,48,7,33,190,189,0,0,12,96,28,23,197,23];
%     nSessions = length(Cassius_Sig_per_Session);
% 
%     plot(1:nSessions, Cassius_Sig_per_Session, 'LineWidth', 3);
%     xlabel('Session');
%     ylabel('Number of Modulated Channel Pairs');
%     ylim([0 200]);
% 
%     % Add ticks at each session
%     xticks(1:nSessions);
%     xlim([0 nSessions+1]);  % keep some spacing on both ends
% end
% 
% %MrCassius_beta_testTone = [3,1,3,3,4,3,3,39,8,14]
% %MrM_beta_testTone = [48,17,23,4,40]


%% Binomial test probability heatmap

ChanPair_Array_Coh_binom = nan(20,20);

for r = 1:20
    for c = 1:20
        k = ChanPair_Array_Coh_sig(r,c);    % successes
        n = ChanPair_Array_Coh_total(r,c);  % total trials
        if n > 0
            % Probability under null (p=0.05) of >=k successes
            p_binom = 1 - binocdf(k-1, n, 0.05);
            % Store as -log10 for visualization (larger = stronger evidence)
            ChanPair_Array_Coh_binom(r,c) = -log10(p_binom);
        end
    end
end

% Precompute sizes
s = [20 20];

% Mask + plotting data
plotData = ChanPair_Array_Coh_binom;
L10_thresh = -log10(0.05);
epsL = 1e-12;
mask = (plotData < (L10_thresh - epsL)) | isnan(plotData);
plotData(mask) = NaN;

% Plot
figure(4); clf

% Global font formatting
set(groot, 'defaultAxesFontName','Arial', ...
           'defaultTextFontName','Arial', ...
           'defaultAxesFontSize',16, ...
           'defaultTextFontSize',16);

hImg = imagesc(xrange, yrange, plotData);
set(hImg,'Interpolation','nearest');  % crisp pixels
yc = fliplr(linspace(yrange(1), yrange(2), s(1)));  % flip to match image top-down
hold on;

clim([L10_thresh, max(ChanPair_Array_Coh_binom(:), [], 'omitnan')]);

% Cell geometry (axis coords)
dx = diff(xrange)/(s(2)-1);
dy = diff(yrange)/(s(1)-1);
xc = linspace(xrange(1), xrange(2), s(2));
yc = linspace(yrange(1), yrange(2), s(1));

% Overlay dark navy for subthreshold/NaN
darkNavy = [0 0 0.35];
for r = 1:s(1)
    for c = 1:s(2)
        if mask(r,c)
            rectangle('Position', [xc(c)-dx/2, yc(r)-dy/2, dx, dy], ...
                      'FaceColor', darkNavy, 'EdgeColor', 'none');
        end
    end
end

% Grid overlay
xg = linspace(xrange(1)-dx/2, xrange(2)+dx/2, s(2)+1);
yg = linspace(yrange(1)-dy/2, yrange(2)+dy/2, s(1)+1);
hm = mesh(xg, yg, zeros(s+1));
hm.FaceColor = 'none'; hm.EdgeColor = 'k'; hm.LineWidth = 0.5; hm.EdgeAlpha = 0.8;

% Axis ticks match count plot
xticks(linspace(xrange(1), xrange(2), s(2)));  xticklabels(1:s(2));
yticks(linspace(yrange(1), yrange(2), s(1)));  yticklabels(1:s(1));

% Colorbar & labels
cb = colorbar;
cb.FontName   = 'Arial';
cb.FontSize   = 16;
cb.FontWeight = 'bold';
ylabel(cb, '-log_{10}(p) Value', ...
       'FontName','Arial','FontSize',16,'FontWeight','bold');
cb.Label.Rotation = 270;

xlabel('AC Channels', 'FontName','Arial','FontSize',16,'FontWeight','bold');
ylabel('PFC Channels','FontName','Arial','FontSize',16,'FontWeight','bold');

ax = gca;
ax.FontName   = 'Arial';
ax.FontSize   = 16;
ax.FontWeight = 'bold';
box on; axis square;

hold off;

% --- export (sharp at 600 dpi, 11x8 in) ---
set(gcf,'Units','inches','Position',[0 0 11 8]); % final size
exportgraphics(gcf, 'Counts_11x8in_600dpi_bold.tif', 'Resolution', 600);


%% Theta Statistics [DO NOT DELETE]

% % for theta bottom and top rows 
% 
% % Apply threshold: keep only values above 1.3
% thresh = 1.3;
% 
% top_vals = ChanPair_Array_Coh_binom(1:10, :);
% bot_vals = ChanPair_Array_Coh_binom(11:20, :);
% 
% % Flatten and filter
% top_vals = top_vals(~isnan(top_vals) & top_vals >= thresh);
% bot_vals = bot_vals(~isnan(bot_vals) & bot_vals >= thresh);
% 
% % Wilcoxon rank-sum test
% [p,~,stats] = ranksum(top_vals, bot_vals);
% fprintf('Wilcoxon rank-sum (thresholded): z = %.2f, p = %.4f\n', stats.zval, p);
% 
% % Bootstrap mean difference
% nBoot = 10000;
% bootDiff = nan(nBoot,1);
% 
% for b = 1:nBoot
%     resamp_top = datasample(top_vals, numel(top_vals));
%     resamp_bot = datasample(bot_vals, numel(bot_vals));
%     bootDiff(b) = mean(resamp_top) - mean(resamp_bot);
% end
% 
% ci = prctile(bootDiff,[2.5 97.5]);
% meanDiff = mean(top_vals) - mean(bot_vals);
% 
% fprintf('Mean difference (top - bottom, thresholded) = %.3f\n', meanDiff);
% fprintf('95%% bootstrap CI = [%.3f, %.3f]\n', ci(1), ci(2));



% %% --- Save selected variables to a specified directory --- [DO NOT DELETE]
% 
% % Specify output directory
% saveDir = 'C:\Users\Corey Roach\Desktop\Maris_Data';  % <-- CHANGE THIS to your folder
% 
% % Make sure directory exists; if not, create it
% if ~exist(saveDir, 'dir')
%     mkdir(saveDir);
% end
% 
% % Specify output filename (without extension)
% outputFileName = 'MrM_beta';  % <-- CHANGE THIS to your desired name
% 
% % Combine directory and filename into full path
% fullFileName = fullfile(saveDir, [outputFileName, '.mat']);
% 
% % List of variables to save
% varsToSave = {'ChanPair_Array_Coh_sig', 'ChanPair_Array_Coh_total', 'ChanPair_Array_Coh_binom'};
% 
% % Save the variables
% save(fullFileName, varsToSave{:});
% 
% fprintf('Saved variables to file: %s\n', fullFileName);



%% --- helper (inline) ---

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end


