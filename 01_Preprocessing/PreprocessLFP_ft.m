%% Info

% 2021_2022: Banu to preprocess the MrCassius data that had been precut by Lalitta. 
% 2023_12_01: Yale used the same pipeline to preprocess the data from MrMiyagi.
% 2024_01_01: Corey uses the code Yale got running to further filter the
% residual line noise out of the Cassius data set 

% modify data structure for fieldtrip toolbox
% tutorial from: www.fieldtriptoolbox.org/tutorial/connectivity
% use fieldtrip functions for preprocessing (resampling, filtering, etc.)

%%
clear all

global DATA_DIR SAVE_DIR Animal 

% set path
DATA_DIR    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\01_LALITTA_Cut\MrCassius';
SAVE_DIR    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed\03_Preprocessed_Data_Corey';

%% Add Fieldtrip and Chronux 

% The fieldtrip and chronux toolboxes are needed to run this code. If you
% do not have these in your MATLAB environment download and use this chunk

% TOOLBOX_DIR = '/home/taku/Documents/MATLAB';
% addpath(genpath(fullfile(TOOLBOX_DIR,'chronux_2_12','chronux_2_12')));
% addpath /Users/yalecohen/Library/CloudStorage/Dropbox/yec/DataAnalysis/fieldtrip-20231215;
% addpath(genpath((fullfile(TOOLBOX_DIR,'fieldtrip','connectivity'))));


%% make sure to use Matlab built-in filtfilt function!
rmpath C:\Users\Corey Roach\Documents\fieldtrip-20230613\external\signal 

%% define data to analyze
Animal  =  'MrCassius';                 %Options: 'MrCassius','MrM'; '
Epoch   = {'moveOnset','preCueOnset','testToneOnset'};

% MrCassius Data Set 

% RecDate = {'190324',...
%             '190326','190328','190330','190331','190404',...
%             '190406','190413','190414','190416','190418',...
%             '190419','190421','190423','190429','190515',...
%             '190517','190520','190522','190526','190531',...
%             '190603','190605','190703','190705','190711',...
%             '190713','190715','190718','190720','190721',...
%             '190723','190725','190728','190917','190918',...
%             '190921','190922','190924','190925','190926',...
%             '190901','190903','190904','190905','190402'};


% MrM Data Set 
% RecDate = {'190322',...
%             '190325','190415','190417','190420','190422',...
%             '190424','190425','190427','190428','190502',...
%             '190507','190509','190514','190516','190518',...
%             '190521','190525','190527','190530','190601',...
%             '190604','190704','190709','190712','190714',...
%             '190717','190719','190722','190729','190815',...
%             '190817','190818','190820','190821','190823',...
%             '190825','190826','190829','190830','190831',...
%             '190901','190903','190904','190905'};


nSession = numel(RecDate);

for ii=1:nSession
    for jj=1:3
        FormatLFP_ft_v2(Animal,DATA_DIR, SAVE_DIR, RecDate{ii},Epoch{jj});
    end
end