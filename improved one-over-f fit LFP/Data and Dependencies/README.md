Dependencies are .m files that the code calls on.

If the user wants to use pre-made seeds (i.e., to have the same results across different sessions of MATLAB), they must have the three seeds*.mat files in the current folder. These pre-made seeds (.mat files) work for any of our datasets, but they are too large to upload on GitHub. We have uploaded them onto FigShare, so the user will need to download from FigShare and have the seed files in the current folder. The number of cells within each seed is determined by the number of heirarchical resamples needed for the bootstrapping process, which is determined by sample size.

Alternatively, if the user decides to make their own seeds, the runtime of makeseeds.m depends on the number of samples in the dataset, but should not exceed 100 sec. They could then save these seeds to have their own reproducible results.

Data are .mat files which are cell arrays. The number of cells is the nuber of subjects, and each cell contains a 2-d array (number of epochs x length of timeseries*sampling frequency). Each subject usually has 28-30 epochs. 
