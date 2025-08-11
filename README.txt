Communication_Toolbox

Last Updated: 06/15/25

This is the code that CR and SF developed to calculate communication between PFC and AC as a function of prior information (i.e, often designated as OnlyPrior and OnlyPretone trials). More specifically, it takes oscillation (LFP) data held in the 00_DATA folder and calculates coherence, granger, and xcorr between significantly modulated channel pairs. Each directory (both code and data) are amended with descriptions of what is held within. This document serves as a high level description of the code to put everything into context. 

Dependencies
MATLAB 2020 or later
FieldTrip (20230613 release or later) 

Introduction

The main engine underneath the hood of this code is the FieldTrip Toolbox. Because of this, our data is often put into a unique format with time, data, and trial structure information held in different fields within a .mat structure. The data can be found in the directory below

D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA

This code begins to process the data after "LALITTA_Cut" step, where Lalitta took the raw tdt data, threw out incomplete trials, and aligned the data to specific EPOC (namely testToneOnset, preCueOnset, and movementOnset).

See the \\00_DATA\01_Lalitta_Cut for more details, but basically all the code that cuts the raw like Lalitta can either be found in the PennBox that Yale established in late 2023 (see below) or zipped in the zz_Utilities directory

https://www.dropbox.com/scl/fo/xb2cfxn1u0jddoojyzm3a/h?rlkey=hm9h36esn7k9dql652nzst5vm&dl=0 


01_Preprocessing 

This step takes Lalitta's aligned data, and re-references (24 to 20 channels, see code documentation for details), and down samples the time series to 1000 Hz 	

02_Spectral Evaluation
Our fundamental unit of measurement is the LFP spectra. This step can visualize the LFP spectra in the frequency domain. In addition, this code band-wise and condition-wise identify channels that are significantly different from their phase scrambled counterpart. The statistical test uses a stock cluster-corrected montecarlo estimation from the fieldtrip tool box (see ScrambleTest.m). The resulting significant channels are used in the rest of the analysis. Important note, CR specifically chose to run the ScrambleTest on the entire data series (not just a discrete epoc). The logic being that we want channels that are overall unmodulated for a long period of time. We care if a channel suddenly goes from modulated to unmodulated in a discete epoc, so it is best to test the entire data frame. Moreover, the longer the data frame the better the spectral resolution, which is ideal for this type of early filtration step. 

03_Behavioral_Analysis 
This step takes behavioral data than Nathan gave us and calculates the chronometric and psychometric curves for the monkey's performance. It also carries out ranksum test between the two conditions at each SNR level. 

04_Maris_Connectivity_Test

The README file within this directory is particular non-descript because in order to understand our code, you need to basically go through this file line-by-line. It runs through all bands, and all significant channels (sometime shared channel pairs when appropriate), and calculates time x frequency matrix (aka spectrogram), xcorr, coherence, and granger, and accompanying null distributions. It generates values for every sig channel or sig channel pair resulting in several thousands of data individual data files. It also carries out the statistical test on both the coherence and granger spectra generating an array of p values. The statistics simply test if the two coherence/granger spectra are different as a function of prior condition. Statistical comparisons are generated using a computationally heavy monte carlo estimation (to account for differences in trials number between conditions), and takes several days to bootstrap all the null distributions

05_Maris_Summary_Statistics 

This folder is the working version of code that generates reportable summaries from the 04_Maris_Connectivity_Test. For Granger and Coherence, it plots the number of significant p-values per channel pair. When organized vertically this metric displays which PFC/AC channels differentially interact as a function prior. SF built code that identify the (first, second, and third) layers groupings that change their communication as a function of prior. Conceptually, this step takes the place of how we originally grouped channels together (superficial, upper-mid, low-mid, and deep). Now we do not have to guess, the current analysis is 'blind' to predetermined laminar demarcation. We identify that channel with  +/- 2 channels, exclude that range, and carry out that process 2 more times. If there are two channels tied for the max, their ranges take up the first and second slots. This code has also been outfitted to calculate cross correlation and phase slope index shared channel pairs

06_Epoch_PSI_Analysis

Contains code that compares PFC-AC connectivity between different task epochs (PreCueOnset and TestToneOnset). PSI and coherence are calculated following trial subsampling, and compared with non-parametric statistical testing. PSI and Coherence are visualized both spatially (i.e. by Channel-Pair via heat map) or non-spatially (i.e., by condition via histogram). 

07_Expectation_PSI_Analysis

Same logic as 06_Epoch_Analysis, hard coded data selection code for convenience. Compares PFC-AC connectivity between OnlyPrior (informative LED, followed by silence) and OnlyPretone (neutral LED, followed by pretones) trials.

08_Congruency_PSI_Analysis

Under construction, but also leverages the same logic as 06_Epoch_anslysis and 07_Expectation_Analysis. Compares PFC-AC connectivity between congruent and noncongruent trials.  

zz_Directed_Connectivity_Summary_Statistics
This directory holds a stripped down wrapper that is used in the Epoch, Expectation, and Congruency Analysis. This is useful if you want to session-wise do a new analysis on preprocessed data. 

zz_XCorr_Analysis
Calculates cross-correlations between PFC and AC time series as a function of expectation or condition. Takes the first, second, and third layer groupings, and seeing which condition(s) varies as a function of the prior condition. To do this, we calculate xcorr for all channels 5 x 5 layer grouping for each condition. This enables us to add a time and an (E/I) descriptor when describing the communication between PFC and AC. We identify the max peak above the null distribution, and develop a histogram with the peak xcorr values for each condition. We then run a chi-squared to test whether those two histograms are different. A pvalue less than .05 suggests that the xcorr max values are sig different between the two conditions. The lag value of max from the two conditions indicates how the two cortical areas are modulated by prior information. 

zz_Utilities 

Functions for carrying out stuff not related to graphs in the paper...
Calculating STRF for each recording location
Checking the status or progress of any scripts that are left to run over several days.  



