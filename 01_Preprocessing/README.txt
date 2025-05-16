01_Preprocessing 

Last Updated: 03/03/25

PreprocessLFP_ft: This is a wrapper or uber file developed by Taku that runs all preprocessing steps on selected sessions. This code takes Lalitta's cut epoch data, and runs 'FormatLFP_ft_v2.m'. That function "baseline corrects" (code: basecorrectLFP) to adjust for an offset in the time series inserted by the Synapse (TDT, acquisition software). The data is also referenced by passing the data through a spatial derivation (code: bipolarLFP). Finally, the data is reformatted for future use in field trip functions (code: downsample_bdLFP_v2) and downsampled to 1000 Hz. 

"Bipolar derivation is a recommended prestep prior to Granger causality and coherence analysis, as the presence of a common reference can lead to spurious results (73, 74). In addition, bipolar derivation enhances the spatial localization of LFP signals and removes the common reference and any common noise or volume conduction in the signal (75). Here, we computed the sample-by-sample bipolar differences by subtracting contacts that were at a distance of 400 μm: next-nearest neighbors for the laminar probe data spaced at 200 μm between contacts, and next-next-nearest neighbors for the probe data spaced at 100 μm between contacts.

We then estimated power, coherence, and Granger causality on these bipolar derivations..." ~ Bastos et al. (2020)

24-channel probe, the bipolar-referenced channel 1 is now computed as:
(1)=V(1)−V(5)

This means that instead of each channel representing the absolute voltage at a specific location, each channel now represents the difference in voltage between two spatially separated sites. This changes how you interpret spatial information in a few key ways:

1. Loss of Absolute Spatial Location
Before bipolar referencing, each channel corresponded to a precise location along the probe.
Now, each channel represents the difference between two locations spaced 4 channels apart.
This means your new bipolar-referenced channel 1 is not "at" channel 1 anymore—instead, it reflects the voltage difference between channels 1 and 5.

Implication:
You no longer have a direct map of neural activity at each electrode position.
Instead, you are looking at a gradient of activity, which can reveal local sources but makes spatial localization less straightforward.

2. Shift in the Reference Frame
Originally, LFP recordings are referenced to a single electrode (e.g., a distant site or a common ground). After bipolar referencing, each channel is referenced to a nearby site instead. This reduces the influence of global signals (e.g., widespread oscillations or movement artifacts) but means the activity you see is relative to the next electrode in the pair. 

Implication: You must now interpret signals in terms of local differences, not absolute voltage.
If there’s synchronous activity across all channels, it cancels out, emphasizing localized sources instead.

3. Reduction of Volume-Conducted Signals
Distant field potentials affect all electrodes similarly, so when you subtract two nearby electrodes, this common signal cancels out.
What remains is a signal dominated by locally generated activity.

Implication:You are now looking more at local current dipoles (neuronal sources and sinks) rather than large-scale LFP oscillations.
This is useful for detecting laminar processing, as activity originating from different cortical layers will stand out more clearly.

4. Reduced Spatial Resolution
Since each new channel is derived from two electrodes, the effective spatial resolution is now coarser.
Your probe originally had 24 unique spatial locations, but after bipolar referencing, you now only have 20.
This slightly blurs fine-grained spatial details, but the trade-off is reduced noise.

Implication: You may lose some precision in pinpointing very small-scale features.
However, you gain a clearer picture of directional flow and local interactions.

Practical Takeaways for Your Analysis
Laminar Positioning: If you were mapping LFP across cortical layers, bipolar referencing alters how you define depth. Instead of a direct position, each channel is now "between" two sites.

PSI Analysis: Since you're analyzing phase slope index (PSI), bipolar referencing may actually enhance true local phase differences by removing global oscillations that could confound results.

Cross-Frequency Coupling: If you’re investigating cross-frequency interactions (e.g., gamma riding on theta), bipolar referencing ensures that observed phase-amplitude coupling isn’t just a global artifact.

CSD Compatibility: Since CSD requires the second spatial derivative, bipolar-referenced signals are less ideal for CSD. If you plan to do CSD later, you might want to work with raw LFP instead.

Summary
Bipolar referencing changes the spatial domain by:

Losing absolute spatial location: Channels now represent voltage differences between two sites, not a single electrode.
Emphasizing local activity: Volume-conducted signals cancel out, highlighting nearby sources.

Altering the reference frame: Signals are now relative to another electrode instead of a common reference.

Reducing spatial resolution: You now have fewer effective channels but cleaner signals.

In practical terms, your 24-channel probe now behaves more like a 20-channel probe, but with better noise suppression and stronger localization of local dipoles. For phase-based analyses, this should help clarify true local phase interactions rather than global signal contamination.

ChoppingData.m: This takes the data from the 02_ft_Preprocessed folder, and cuts out the time window that we care about. The reason that we do this is that fieldtrip has a hard time iteratively calling declared time windows across different functions. It is easiest to cut the data ahead of passing through fieldtrip function, and just fieldtrip read the entire data array. 

AdjustSampleInfo.m: CR used this to change the SampleInfo field after data has been cut with the ChoppingData.m script. Could be consolidated into ChoppingData.m  in future releases 

MoveOnset_Cut.m & MoveOnset_Cut_v2: SF made these two scripts to cut moveOnset data from the 02_ft_Preprocessed folder. The first version mostly adjusts the SampleInfo field after EPOC cutting. 

[To-Do] All the files in this directory should be consolidated into a single function. 

 