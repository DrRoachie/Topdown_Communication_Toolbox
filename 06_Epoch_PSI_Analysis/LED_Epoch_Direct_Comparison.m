% function [preCueOnset_PSI_grid, testToneOnset_PSI_grid, ...
%           preCue_averagedGrid, testTone_averagedGrid, ...
%           PFC_2_AC_PreCue, AC_2_PFC_PreCue, ...
%           PFC_2_AC_TestTone, AC_2_PFC_TestTone, ...
%           group_data] = build_epoch_PSI(Condition, animals, rootdir)

animals             = {'MrM'};
OnlyPrior_rootdir   = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_09_PSI_Epoch_Analysis_OnlyPrior_All_Congruency';
OnlyPretone_rootdir = 'C:\Users\Corey Roach\Documents\00_DATA\2025_09_12_PSI_Epoch_Analysis_OnlyPretone_All_Congruency';

% --- Build (OnlyPrior) ---
[~, ~, ~, ~, ...
 OnlyPrior_PFC_2_AC_PreCue, OnlyPrior_AC_2_PFC_PreCue, ...
 ~, ~, ~] = build_epoch_PSI('OnlyPrior', animals, OnlyPrior_rootdir);

% --- Build (OnlyPretone) ---  (FIXED: condition + rootdir)
[~, ~, ~, ~, ...
 OnlyPretone_PFC_2_AC_PreCue, OnlyPretone_AC_2_PFC_PreCue, ...
 ~, ~, ~] = build_epoch_PSI('OnlyPretone', animals, OnlyPretone_rootdir);

% --- Collapse across direction (PreCue): concatenate |PSI| from both directions ---
prior_pre_all   = [OnlyPrior_PFC_2_AC_PreCue(:);   OnlyPrior_AC_2_PFC_PreCue(:)];
pretone_pre_all = [OnlyPretone_PFC_2_AC_PreCue(:); OnlyPretone_AC_2_PFC_PreCue(:)];

% Clean up NaNs / Infs just in case
prior_pre_all   = prior_pre_all(isfinite(prior_pre_all));
pretone_pre_all = pretone_pre_all(isfinite(pretone_pre_all));

% --- Independent-samples Wilcoxon rank-sum (two-sided) ---
[p_val, ~, stats] = ranksum(prior_pre_all, pretone_pre_all, 'tail', 'both');

% Mann-Whitney U and Cliff's delta (δ) for effect size
n1 = numel(prior_pre_all); 
n2 = numel(pretone_pre_all);
R1 = stats.ranksum;                      % sum of ranks for group 1 (prior_pre_all)
U1 = R1 - n1*(n1+1)/2;                   % Mann-Whitney U for group 1
cliffs_delta = (2*U1)/(n1*n2) - 1;       % δ in [-1, 1]

% --- Quick descriptives (mean ± 95% CI on |PSI|) ---
m1   = mean(prior_pre_all, 'omitnan');   s1 = std(prior_pre_all, [], 'omitnan');   ci1 = 1.96*s1/sqrt(max(1,n1));
m2   = mean(pretone_pre_all, 'omitnan'); s2 = std(pretone_pre_all, [], 'omitnan'); ci2 = 1.96*s2/sqrt(max(1,n2));

fprintf('\nPreCue |PSI| (collapsed directions):\n');
fprintf('OnlyPrior   : n=%d, mean=%.6f, 95%%CI=±%.6f\n', n1, m1, ci1);
fprintf('OnlyPretone : n=%d, mean=%.6f, 95%%CI=±%.6f\n', n2, m2, ci2);
fprintf('Rank-sum (indep): p=%.6g, z=%.3f, Cliff''s δ=%.3f\n', p_val, stats.zval, cliffs_delta);


%%
% ===================== FINAL PLOTTING CODE (4"x4") =======================
% Inputs expected to already exist in workspace:
%   prior_pre_all   : vector of |PSI| values (Informative / OnlyPrior, PreCue; collapsed across directions)
%   pretone_pre_all : vector of |PSI| values (Neutral / OnlyPretone, PreCue; collapsed across directions)
%
% Provide your save directory here:
savedir = 'C:\Users\Corey Roach\Desktop\Figure_Epoch_PSI_Analysis_LED_Comparison';

% ---- Colors (Informative=OnlyPrior=turquoise, Neutral=OnlyPretone=yellow)
informativeColor = [0 0.6 0.6];  % turquoise
neutralColor     = [0.9 0.7 0];  % yellow

% ---- Clean data (safety)
prior_pre_all   = prior_pre_all(isfinite(prior_pre_all));
pretone_pre_all = pretone_pre_all(isfinite(pretone_pre_all));

% ---- Descriptives for mean ± 95% CI
n1  = numel(prior_pre_all);
n2  = numel(pretone_pre_all);
m1  = mean(prior_pre_all, 'omitnan');
m2  = mean(pretone_pre_all, 'omitnan');
s1  = std(prior_pre_all, [], 'omitnan');
s2  = std(pretone_pre_all, [], 'omitnan');
ci1 = 1.96 * s1 / sqrt(max(1, n1));
ci2 = 1.96 * s2 / sqrt(max(1, n2));

% ---- Figure (4x4 inches), publication styling
fig = figure('Color','w', 'Units','inches', 'Position',[1 1 4 4], 'Renderer','painters'); 
ax  = axes(fig); hold(ax, 'on');

% ---- Scatter all points with light jitter
jitter = 0.18;
x1 = 1 + (rand(n1,1)-0.5)*jitter;
x2 = 2 + (rand(n2,1)-0.5)*jitter;

s1p = scatter(ax, x1, prior_pre_all, 12, 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', informativeColor, 'MarkerEdgeColor', 'none');
s2p = scatter(ax, x2, pretone_pre_all, 12, 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', neutralColor, 'MarkerEdgeColor', 'none');

% ---- Mean ± 95% CI overlays + mean markers
% e1 = errorbar(ax, 1, m1, ci1, 'vertical', 'LineWidth', 1.5, 'CapSize', 10, 'Color', informativeColor);
% e2 = errorbar(ax, 2, m2, ci2, 'vertical', 'LineWidth', 1.5, 'CapSize', 10, 'Color', neutralColor);
% plot(ax, 1, m1, 'o', 'MarkerFaceColor', informativeColor, 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
% plot(ax, 2, m2, 'o', 'MarkerFaceColor', neutralColor,     'MarkerEdgeColor', 'k', 'MarkerSize', 5);

% ---- Axes & labels
xlim(ax, [0.5 2.5]);
xticks(ax, [1 2]);
xticklabels(ax, {'Informative','Neutral'});
ylabel(ax, '|PSI| During LED Epoch');
grid(ax, 'off'); box(ax, 'off');
set(ax, 'FontName','Times New Roman', 'FontSize', 11);


% ---- Tight layout inside 4"x4"
outerpos = ax.OuterPosition; 
ti = ax.TightInset; 
left   = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width  = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
set(ax, 'Position', [left bottom ax_width ax_height]);

% ---- Export (super high resolution PNG + vector PDF)
baseName = 'PreCue_theta_Collapsed_AllPoints_4x4in';
pngPath  = fullfile(savedir, [baseName '.png']);
pdfPath  = fullfile(savedir, [baseName '.pdf']);

exportgraphics(fig, pngPath, 'Resolution', 1200);            % super high-res
exportgraphics(fig, pdfPath, 'ContentType', 'vector');       % perfect for print

fprintf('Saved figure to:\n  %s\n  %s\n', pngPath, pdfPath);
% =================== END FINAL PLOTTING CODE (4"x4") =====================


