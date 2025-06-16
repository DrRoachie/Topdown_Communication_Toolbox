From: Lalitta Suriya-Arunroj <lalitta.suriya.arunroj@gmail.com>
Date: February 28, 2025 at 11:30:30 AM EST
To: "Cohen, Yale" <ycohen@pennmedicine.upenn.edu>
Subject: Re: [External] Re: 2 questions

﻿
Hi Yale,

The information on the best frequency are now on this folder BestFreqInfo
There are data from pure tones and STRF (the file names indicate which on is which)
For the STRF the parameters are standard from the STRF code that I got from Jaejin/Taku; they are straightforward.
For the pure tone tuning there are too many parameters that you may need, you may want to focus on the columns:
 
audiRespTuning_PT_su_stimOn_tuningField_tuningCodeFreq  - The fitted best frequency from the tuning curve during stimulus onset
audiRespTuning_PT_su_stimOn_tuningField_tuningLatency   - The latency of the peak activity during stimulus onset
audiRespTuning_PT_su_stimOn_tuningField_tuningWidth     - The width of the tuning curve during stimulus onset
audiRespTuning_PT_su_stimOn_tuningSpk_tuningMaxFreq     - The maximum-response frequency (among 30 frequencies we played) during stimulus onset
audiRespTuning_PT_su_stimOn_tuningSpk_tuningP           - Significance of auditory response selectivity (non-parametric test: kruskal-wallis) during stimulus onset

audiRespTuning_PT_su_stimOff_tuningField_tuningCodeFreq  - The fitted best frequency from the tuning curve after stimulus offset
audiRespTuning_PT_su_stimOff_tuningField_tuningLatency   - The latency of the peak activity after stimulus offset
audiRespTuning_PT_su_stimOff_tuningField_tuningWidth     - The width of the tuning curve after stimulus offset
audiRespTuning_PT_su_stimOff_tuningSpk_tuningMaxFreq     - The maximum-response frequency (among 30 frequencies we played) after stimulus offset
audiRespTuning_PT_su_stimOff_tuningSpk_tuningP           - Significance of auditory response selectivity (non-parametric test: kruskal-wallis) after stimulus offset

Please find attached two example units: one has peak activity during the on-period and another one during the off-period (there are some that show both peaks as well).
Let me know if you need further clarifications. 

All the best,
Lalitta 