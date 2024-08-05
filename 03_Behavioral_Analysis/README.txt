
07_Behavioral_Analysis

Last updated: 07/29/2024

This is code developed by SF to generate and analyze psychometric and chronometric curves. Both functions run animal-wise and condition-wise, meaning for 2 animals and 2 conditions each function must run 4 times. More detailed instruction about running a full Wilcoxon Rank-Sum Test can be found in the Psycho_Chrono_Wilcoxon_Test.m script.


Psychometric_Chronometric_Analysis.m

In general, this is a combined script that primarily filters data and extracts the choice and reaction times for a range of signal-to-noise ratio (SNR) values. The code is divided into 3 main chunks. 

The first chunk is for defining data. Data for psychometric analysis comes from the 04_Epoc_Cut folder and data for chronometric analysis comes from the DDM csv table called '20210511_audiDeci_monkeyBeh_DT_13-Jun-2023.csv' in the 03_DDM_Decision_Times folder. This is also where 	'Animal', 'Epoch', and 'Condition' are specified before running. A note about 'Condition', in this experiment we partition into prior only and pretone only. Trials exist in which there are both priors and pretones which we generally ignore. Throughout the script there are many instances of comparing 'Condition'. This is because for the prior only condition, there are 3 possible states: high LED prior, low LED prior, and neutral LED prior. The pretone only condition only has 2 states: high pretone and low pretone. The purpose of the if statements containing 'Condition' is to make the code robust to the 'Condition' variable.

The second chunk is the psychometric analysis. For each session, a 2D array is generated with prior (H, L, N) as rows and SNR values as columns (9 SNR values total). The array is populated with the proportion of choice high for each corresponding prior and SNR value. The array is appended onto a 3D array that saves all sessions. The code then takes the mean and standard error of the 3D array, collapsing across session to produce two final 2D arrays. Each row of the 2D arrays corresponding to a prior state (H, L, N) is plotted with their standard error to generate the psychometric curve.

The third chunk is the chronometric analysis. Due to naming differences, defined data is adjusted to match the csv table. Then, similar to the psychometric chunk, a 2D array is populated for each session with the reaction times of correct trials and added to the 3D array. After collapsing across sessions for mean and standard error, the chronometric curve is plotted.

**Note: the variable 'sessions' is used dynamically in the psychometric and chronometric chunks due to the differing structures of the Epoc_Cut folder and the DDM csv. For psychometric filtering, sessions is a directory of folders 190330, 190404, etc. In the chronometric chunk, sessions is a list of session identifiers that correspond to a session date but contain values 1, 2, 3, etc. Information about matching session identifiers with dates can be found in the 'audiDeci_monkey_sess_date_13-Jun-2023.csv' file in 03_DDM_Decision_Times.


Psycho_Chrono_Wilcoxon_Test.m

This code performs the Wilcoxon Rank-Sum Test on the psychometric and chronometric data. Importantly, the script relies on the variables 'total_proportion_H' and 'total_RTs' calcualted in Psychometric_Chronometric_Analysis.m to run the test for psychometric and chronometric, respectively. The first two chunks run the rank-sum test on psychometric and chronometric data. In the prior only condition, the test is done on high vs. neutral data and low vs. neutral data. In the pretone only condition, the test is done on high vs. low data. The last chunk saves the 6 p-value results of the tests into a table 'T' which can be saved manually. One table is generated for each animal and the last chunk cannot save successfully until both conditions are run, requiring separate calls to Psychometric_Chronometric_Analysis.m and Psycho_Chrono_Wilcoxon_Test.m.

