Dependencies are .m files that the code calls on.

Data are .mat files which are cell arrays. The number of cells is the nuber of subjects, and each cell contains a 2-d array (number of epochs x length of timeseries*sampling frequency). Each subject usually has 20 epochs. Datasets from both V1 and S1 are included and can be run by simply adjusting the filename in run_improved_EEG_fit.m as needed.

