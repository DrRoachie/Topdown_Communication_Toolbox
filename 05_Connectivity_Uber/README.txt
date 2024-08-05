README

05_Connectivity_Uber 

Last Updated: 07/30/24

This is probably the most important directory in the Communication_ToolBox 

Pairwise_Connectivity_Test_Uber: This is our main uber file that calculates the spectrogram, coherence, granger causality, and xcorr for each condition (OnlyPretone and OnlyPrior) for every significant channel. For coherence and granger, it runs statistical test between the shared channel-pairs between the two conditions, and outputs a savable array of pvalues for each test. Aside from looping through channel pairs, a major innovation in this code is that it calculates the coherence and granger spectras off the same fourier/cross-spectral density step radically improving the efficacy of the code. 

Pairwise_Connectivity_Wrapper: This was provided to Noga (a visiting Israeli student in our lab) as a launch pad for developing new connectivity metrics that can be incorporated in our pipeline. Basically it just the wrapper for looping through the shared channel-pairs within a folder of significant channels, and calculating a connectivity metric. In this example, we leave in the example function 'TrialbyTrialXcorr' to demonstrate how the code works 