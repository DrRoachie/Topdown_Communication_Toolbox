
08_Summary_Statistics

Last updated: 09/09/2025

This folder contains scripts to summarize and visualize the output of the Maris test code.

Summary_Figure_v3.m

Loops through a directory of sessions and generates complete summary figures for all sessions, animals, frequency bands, and pairs. The idea is that this would generate a library of summary figures by going through all the data within a directory. 


Summary_PValue_Spatial_theta_Figure_v3.m
Summary_PValue_Spatial_beta_Figure_v3.m


For a given frequency band, list of animals, and statistic (coherence), this script counts the number of passing montecarlo p-values for that statistic for each channel pair into an array. A passing p-value means that for that channel pair and session, the statistic varied as a function of condition. The array can be plotted as a colored heatmap grid or as a network diagram (refer to chunk titles). The heatmap can be just as well plotted with the Matlab 'heatmap' function on the corresponding arrays, although color will be lost. 

In addition to the heatmap of the raw counts, it also calculates a the results of the binomial test on each channel pair, and returns a -log(10)p value indicating how likelihood of observing that number of counts considering chance. 

For convenience this code is forked for each frequency band because some parameters are hard coded and band-specific. 


