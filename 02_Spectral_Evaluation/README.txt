02_Spectral_Evaluation 

Last Updated: 06/16/25

After a meeting with Bejian Pesaran at the end of February, we accepted that the LFP over trials should be the basic unit that we build all other analyses. Bejian urged us to do a spectral evaluation on each channel before feeding our data into the coherence analysis to identify channels with a dynamic spectral composition (i.e., compare the channel or epoc to an internal baseline). 

FieldTrip has a built in function for comparing time-frequency series using the monte carlo method, which seemed suitable for our needs. However, like all monte carlo estimations the functions require the two conditions to be the same length in time. This was problematic because the only time period where nothing is happening (i.e., the monkey is either doing nothing or anticipating the next trial) was before the LED was illuminated, and was only 100 ms long with a maximum spectra resolution of 10 Hz. Yale came up with the idea  to phase scramble each channel (i.e. generate random spectra from the data) and compare that data with the non-scrambled data.


ChanWise_LFPeval.m: This code compares the LFP spectra between two conditions, and can also generate a spectrogram (optional). For each channel pair, it asks how often it passes the phase shuffle filter and/or evoked potential filter, pulls the power when channel is significant, and finds if the probability of significance is correlated with the observed avg power on that channel. 

LFP_Spectral_Analysis.m: [Outdated] This routine will plot the spectra in the frequency domain with one of our preselected layer grouping and epocs. For each area, the code will return a simple line plot of the avg power at each frequency surrounded by a 95% CI. It also separates data by right and wrong trials.

SpecBalanced_Eval.m: This script, developed by SB, session-wise, identifies what channels spectrally overlap (i.e., exhibit similiar power values), and outputs the channel labels by calling SpecBalancedRangedBased.

ChanWise_PhaseShuffle.m: For each session, this code channel-wise compares the spectral properties of the channel with its phase scrambled counterpart using a monte carlo estimation. It outputs the indices of channels that pass the test, which can then be called by the coherence code further down the pipe. This code takes a 20 minutes to 2 hours depending on which EPOC and FrequencyBand you are analyzing.



