README

04_Coherence_and_Granger

If you came here looking to carry out the channel-pair wise version of the analysis, this is not the folder to you. This step carries out a non-parametric between conditions for coherence. It was incorporated into a more comprehensive wrapper that calculates more connectivity metrics (xcorr and granger causality). This code was our working version of our code before we made the switch to a channel-pair wise approach and rebuilt the code to be more efficient. For each condition (OnlyPrior and OnlyPretone), this code will band-wise and layer-grouping wise calculate coherence, then run a monte carlo estimation between the two conditions. It ulimateluy outputs an array of pvalues for each combination. 

It was retained within this folder architecture because it still collapses across all combinations of channel pairs within a layer grouping and those chunks might be of value to someone. 

Coherence_Trial_Shuffle_Analysis Directory: The basic version of the code. It basically will calculate all the coherence metrics for a single session or directory. It will only calculate one layer FrequencyBand and LayerGrouping per run. This code represents a development step. Can still be used for conviencence if you want the coherence for a single band. 

Maris_Statistical_Test_Uber: This code basically takes the Coherence_Trial_Shuffle_Analysis procedure and repeats it on each FrequencyBand and the LayerGrouping combo for a single session, and calculates the coherence for correct and wrong trials. I

max_cluster_plot: Plots a distribution of clustermax values as a function of session. Animals are designated as different shapes and significance is coded as green or black. 
 
[To-Do] These summary plots might be useful if you only plot the significant values. Any apparent effect is blurred by including the wrong trials on the same plot. 

max_cluster_visualization: Similiar to max_cluster_plot. This takes the max cluster value of each session and plots it on the y. At this point, we were interested in how correct and wrong trials differed, so behavior is on the x axis. The correct and wrong trials for each session are connected by a line, so you can see trends on how the cluster values change as a function of whether the monkey got it wrong or not. The data points are coded with green values representing sessions where OnlyPrior and OnlyPretone trials were significantly different.  In retrospect, this was an incorrect or precocious comparison to the data because it does not really tell us much. 

Examining the data by pvalue and max cluster value was not very useful. We ended up axing the logic, and moving towards other methods for the visualizing the data. 

[To-Do] These summary plots might be useful if you only plot the significant values. Any apparent effect is blurred by including the wrong trials on the same plot. 

[To-Do] Because the fourier step is the same for both coherence and granger, this code is easily outfitted for granger calculations. Simply change the cfg.method from 'coh' to 'granger' and rename any variable name referring to coherence with granger 