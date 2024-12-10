
08_Summary_Statistics

Last updated: 07/29/2024

This folder contains scripts by SF to plot the statistics on processed data in the 2024_07_24_Analysis folder.

Summary_Figure.m

For a given animal, session, epoch, channel pair, frequency band, and behavior, this code plots the spectrogram of PFC and AC in both conditions, coherence and bidirectional granger (PFC to AC and AC to PFC) with their respective Maris statistical test results, and the cross correlation.

Summary_Figure_v2.m

Functionally identical to Summary_Figure.m but plots everything as a jpeg rather than a Matlab subplot. Individual figures end up being much smaller and more difficult to read.

Summary_Figure_v3.m

Loops through a directory of sessions and generates complete summary figures for all sessions, animals, frequency bands, and pairs. The idea is that this would generate a library of summary figures by going through all the data within a directory. 

Summary_PValue_Spatial_Figure.m

For a given frequency band, list of animals, and statistic (coherence or granger), this script counts the number of passing montecarlo p-values for that statistic for each channel pair into an array. A passing p-value means that for that channel pair and session, the statistic varied as a function of condition. The array can be plotted as a colored heatmap grid or as a network diagram (refer to chunk titles). The heatmap can be just as well plotted with the Matlab 'heatmap' function on the corresponding arrays, although color will be lost. This script is best run chunk-by-chunk. If Coherence is selected it generates one key variable ChanPair_Array_Coh, if Granger is selected it generates two key varibles that specify directionality 'ChanPair_Array_PFC_AC' and 'ChanPair_Array_AC_PFC'

Summary_PValue_Blur_Cluster_Detection.m

To identify which PFC and AC channels showed the most activity, an array from Summary_PValue_Spatial_Figure.m is blurred 3 times. This script uses a custom Matlab function 'blurArray' that decreases the size of the array by 1 and averages each 2x2 window of values. This creates a blurring effect from which peak clusters can be identified. The highest 3 clusters (defined as the 5x5 area centered around a maximum value) are saved in an array which contains information on peak rank (highest, second highest, third highest), sender channel number, receiver channel number, logical id for which channel pair is the maximum, and the number of p-values in the original unblurred array. Overlapping clusters are handled by removing duplicate channel pairs so all channel pairs only appear once. These values are then fed into the Xcorr pipeline, so that you can analyze the timing and areal nature of the communication as a function of prior. 

[To-Do] Summary_Figure_v3 needs some developing. specifically, it needs to be tested for several channel pairs within a directory and also get a saving or output chunk written into it.
[To-Do][high priority] Summary_PValue_Blur_Cluster_Detection.m does not output and we have not fed those values into the xcorr pipeline yet. 
[To-Do?} Could maybe add a save chunk in Summary_PValue_Spatial_Figure.m so that it does not always need to be run before Summary_PValue_Blur_Cluster_Detection.m. Cluster detection script only needs variables ChanPair_Array_PFC_AC and ChanPair_Array_AC_PFC 
