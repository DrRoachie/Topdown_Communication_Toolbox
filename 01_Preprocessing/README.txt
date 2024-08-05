01_Preprocessing 

Last Updated: 7/23/24

The  PreprocessLFP_ft was developed by Taku and takes the Lalitta cut data, downsamples to 1000 Hz, adjusts the intrinsic lag coming out of TDT LFP stream, filters out the base 60 Hz noise plus harmonics, and organizes the data for future analyses centered around the fieldtrip tool box. For more information about preprocessing the data see D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed

ChoppingData.m: This takes the data from the 02_ft_Preprocessed folder, and cuts out the time window that we care about. The reason that we do this is that fieldtrip has a hard time iteratively calling declared time windows across different functions. It is easiest to cut the data ahead of passing through fieldtrip function, and just fieldtrip read the entire data array. 

AdjustSampleInfo.m: CR used this to change the SampleInfo field after data has been cut with the ChoppingData.m script. Could be consolidated into ChoppingData.m  in future releases 

MoveOnset_Cut.m & MoveOnset_Cut_v2: SF made these two scripts to cut moveOnset data from the 02_ft_Preprocessed folder. The first version mostly adjusts the SampleInfo field after EPOC cutting. 

[To-Do] All the files in this directory should be consolidated into a single function. 

 