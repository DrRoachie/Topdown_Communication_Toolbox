02_Spectral_Evaluation 

Last Updated: 09/09/25

This directory contains a scripts that visualize and compare the raw LFP and evoked potential across different experimental  conditions. Requires you feed in un-cut preprocessed data.

ChanWise_LFP_Plot: This basically just the average LFP power across across all significant channels, and does basic stimulation. 


ChanWise_LFPeval_Epoch.m; ChanWise_LFPeval_Expectation: This code compares the LFP spectra between two conditions, and can also generate a spectrogram (optional). For each channel pair, it asks how often it passes the phase shuffle filter and/or evoked potential filter, pulls the power when channel is significant, and finds if the probability of significance is correlated with the observed avg power on that channel. 

Spectrogram_Analysis_ScratchPad_v2: Moves across all shared significant channel pairs between two conditions, and calculates the mean evoked response. 

Spectrogram_Analysis_ScratchPad_v3: Averages across all levels to get a single times series per condition, also plots the spectrograms for each channel collapsed across sessions, does some stats. 

Spectrogram_Analysis_ScratchPad_v4: Updated to collapse separate responses by high and low. 




