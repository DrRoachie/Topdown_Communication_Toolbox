README

05_Connectivity_Uber 

Last Updated: 04/21/25
CR

This is probably the most important directory in the Communication_ToolBox 

Pairwise_Connectivity_Test_Uber: This is our main uber file that calculates the spectrogram, coherence, granger causality, and xcorr for each condition (OnlyPretone and OnlyPrior) for every significant channel. For coherence and granger, it runs statistical test between the shared channel-pairs between the two conditions, and outputs a savable array of pvalues for each test. Aside from looping through channel pairs, a major innovation in this code is that it calculates the coherence and granger spectras off the same fourier/cross-spectral density step radically improving the efficacy of the code. 

Pairwise_Connectivity_Test_Uber_v2: A stripped down version of Pairwise_Connectivity_Test_Uber_v2 that cuts out generating and saving all the plots, and simply outputs the XCorr analysis and the Maris Test Results. 

Pairwise_Connectivity_Test_Uber_v3: this version of the connectivity test can run the Maris Procedure on the data filtered by SNR and congruency. It also outputs xcorr and phase slope index for each shared channel pair between the two conditions (depending on the exp., this could just be channels from one condition used for both). T