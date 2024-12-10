03 Cross Corr

Last Updated: 07/30/2024

XCorr_Mean.m

This script loops through sessions in the 04_Epoc_Cut folder and calculates xcorr and null with confidence intervals. The average across sessions is plotted.

XCorr_Peak_Histogram.m

The first chunk populates an array for a list of animals, conditions, and frequency bands with the time bins of passing xcorr peaks (defined as non-overlapping 95 confidence intervals between the xcorr and the null at the maximum xcorr value). This is done for all sessions where the pool of channel pairs for each session comes from the significant channels. An structure is saved for each condition with fields containing the frequency bands. Each frequency band contains a full array of all significant channels and a filtered array of only passing xcorr peaks. A saved copy of the two PriorOnly and PretoneOnly structures can be found in the XCorr_Histogram_Data folder in 2024_07_24_Analysis.

The second chunk loads the structures in preparation for plotting the histograms.

The third chunk plots the histograms. The 'filter_for_sign' variable indicates whether or not to separate positive and negative xcorr peaks. 'filter_for_sign' is usually true.

The fourth chunk plots the histograms by layer grouping with predefined layer combinations. The code only handles PFC superficial to AC superficial, PFC lowermid to AC uppermid, PFC deep to AC deep, and PFC deep to AC uppermid. An error is thrown if an unsupported layer group pair is inputted.

The fifth chunk plots the histograms by animal. Each session corresponds to an animal, so the arrays to be plotted are filtered by session number.

XCorr_Peak_Histogram_v2.m

This script is identical to XCorr_Peak_Histogram.m in chunk function, but instead of calculating xcorr on all significant channels, xcorr is calculated on the channels of max pvalue count selected in the Summary_PValue_Blur_Cluster_Detection.m script. This takes a while to run and as of 07/30/2024 has not been run/saved the structures.

XCorr_Peak_Histogram_Plot_Max.m

This is an alternative way to see the xcorr of the max pvalue channels. Similar to plotting by layer grouping, this script filters the bins down to only the channel pairs included in the arrays from Summary_PValue_Blur_Cluster_Detection.m. Importantly, the xcorr is still calculated using all significant channels (like in XCorr_Peak_Histogram.m).

XCorr_Chi_Square_Test.m

Unfinished script meant to run the chi square test of independence on xcorr histograms. First two chunks load the histogram data and generate histogram counts for both conditions. The third chunk runs the chi square test.
