
% close all;

%% Define data
Animal      = 'MrCassius';
RecDate     = '190419';
Epoch       = 'testToneOnset';
Send_Chan   = "PFC_ch11";
Rec_Chan    = 'AC_ch03';
Freq_Band   = 'theta';
Behavior    = 'Correct';

datadir     = fullfile('D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_07_01_Analysis', RecDate, Animal, extractBefore(Epoch, 'Onset'));

%%
summary_fig = figure;

% super title
Send_Chan_lb   = strrep(Send_Chan,'_ch','');
Rec_Chan_lb    = strrep(Rec_Chan,'_ch','');

sgtitle(sprintf('%s %s %s-%s Connectivity Analysis in %s', Animal, RecDate, Send_Chan_lb, Rec_Chan_lb, Freq_Band));


% Prior Spectrogram
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_Prior_Spectrogram.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.png');

figure(summary_fig);
s = subplot(4,4,1);

img = imread('temp.png');
imshow(img);

delete('temp.png');


% Pretone Spectrogram
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_Pretone_Spectrogram.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.png');

figure(summary_fig);
subplot(4,4,2);

img = imread('temp.png');
imshow(img);
axis image

delete('temp.png');


% Coherence Spectrum
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_Spec.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
subplot(4,6,7);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');


% Coherence zthresh
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_zthresh.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
subplot(4,6,8);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');


% Coherence montecarlo
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_montecarlo.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
subplot(4,6,9);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');

% add p value
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Coh_data.mat', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
load(fullfile(datadir, file_name), 'pvalue_c');
p_value_text = sprintf('p = %.4f', pvalue_c.(Freq_Band));
annotation('textbox', [0.46, 0.59, 0.1, 0.1], 'String', p_value_text, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8, 'FitBoxToText', 'on');


% Granger Spectrum PFC to AC
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_Spec_PFC_to_AC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
s1 = subplot(4,6,13);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');


% Granger zthresh PFC to AC
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_zthresh_PFC_to_AC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
s1 = subplot(4,6,14);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');


% Granger montecarlo PFC to AC
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_montecarlo_PFC_to_AC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
s1 = subplot(4,6,15);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');

% add p value
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
load(fullfile(datadir, file_name), 'pvalue_PFC_AC');
p_value_text = sprintf('p = %.4f', pvalue_PFC_AC.(Freq_Band));
annotation('textbox', [0.46, 0.37, 0.1, 0.1], 'String', p_value_text, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8, 'FitBoxToText', 'on');


% Granger Spectrum AC to PFC
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_Spec_AC_to_PFC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.png');

figure(summary_fig);
subplot(4,6,19);

img = imread('temp.png');
imshow(img);
axis image

delete('temp.png');


% Granger zthresh AC to PFC
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_zthresh_AC_to_PFC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.png');

figure(summary_fig);
subplot(4,6,20);

img = imread('temp.png');
imshow(img);
axis image

delete('temp.png');


% Granger montecarlo AC to PFC
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_montecarlo_AC_to_PFC.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.png');

figure(summary_fig);
subplot(4,6,21);

img = imread('temp.png');
imshow(img);
axis image

delete('temp.png');


% add p value
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);
load(fullfile(datadir, file_name), 'pvalue_AC_PFC');
p_value_text = sprintf('p = %.4f', pvalue_AC_PFC.(Freq_Band));
annotation('textbox', [0.46, 0.16, 0.1, 0.1], 'String', p_value_text, 'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8, 'FitBoxToText', 'on');


% XCorr
file_name   = sprintf('%s_%s_%s_%s_%s_%s_%s_XCorr.fig', Animal, RecDate, Epoch, Send_Chan, Rec_Chan, Freq_Band, Behavior);

h = openfig(fullfile(datadir, file_name),'invisible');
saveas(h, 'temp.jpg');

figure(summary_fig);
subplot(2,2,2);

img = imread('temp.jpg');
imshow(img);

delete('temp.jpg');


