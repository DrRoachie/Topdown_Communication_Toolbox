02_Phase_Shuffle_Filter
Updated: 09/09/25
CR

The first step of our pipeline is identifying channels that are stimulus modulated. To do this, for any given condition, we compare the data to a phase shuffled match. We run tests independently for each frequency, epoch, and behavior. The statistical comparison between the data and phase-shuffled match is carried out by a function called 'ScambleTest.m'.

ChanWise_SpecEval.m: The original script that is ran manually for each session on bipolar re-referenced data

ChanWise_SpecEval_fullArray.m: Developed in the July 2025 to repeat the analysis on the non-referenced original data set. This improves on the original variant in that it will automatically run through all sessions in the folder. 

*********************************************************************************************************************************

NOTES

After a meeting with Bejian Pesaran at the end of February 2024, we accepted that the LFP over trials should be the basic unit that we build all other analyses. Bejian urged us to do a spectral evaluation on each channel before feeding our data into the coherence analysis to identify channels with a dynamic spectral composition (i.e., compare the channel or epoc to an internal baseline). 

FieldTrip has a built in function for comparing time-frequency series using the monte carlo method, which seemed suitable for our needs. However, like all monte carlo estimations the functions require the two conditions to be the same length in time. This was problematic because the only time period where nothing is happening (i.e., the monkey is either doing nothing or anticipating the next trial) was before the LED was illuminated, and was only 100 ms long with a maximum spectra resolution of 10 Hz. Yale came up with the idea  to phase scramble each channel (i.e. generate random spectra from the data) and compare that data with the non-scrambled data.

