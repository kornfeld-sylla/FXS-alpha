Dependencies are .m files that the code calls on.

IF the user decides to make their own seeds, the runtime of makeseeds.m depends on the number of samples in the dataset, but should not exceed 100 sec.

Data are .mat files which are cell arrays. The number of cells is the nuber of subjects, and each cell contains a 2-d array (number of epochs x length of timeseries*sampling frequency). Each subject usually has 20 epochs. Datasets from both V1 is included. Other datafiles can be accessed on FigShare and added to the current folder, and can be analyzed simply by uncommenting the filenames of interest in run_improved_EEG_fit.m as needed.

