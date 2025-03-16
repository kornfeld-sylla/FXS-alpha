Code for analysis of the temporal dynamics of a periodic signal. A periodic signal can either be relatively continuous or bursty, and this code quantifies the number of bursts and the burst duration within the frequency range of the periodic signal.

For analysis of alpha-like oscillations in mice, the frequency range of interest is 2-10 Hz. 

This code requires MATLAB and the signal processing toolbox in MATLAB. It was developed on version R2018b. 

To use this code, you can simply execute run_burst_analysis.m. To execute, you need to have all the dependencies in the path, and the data in the current folder. 

Two user-defined parameters are listed at the top. They do not need to be adjusted for the code to execute properly, but are modifiable parameters. The first parameter indicates if the data you are analyzing is EEG or LFP (the default here is LFP). The second parameter allows you to adjust the threshold for defining a burst (default 0.9, or 90th percentile of the amplitude values). 

Dependencies: This code calls on the BURST_ANALYSIS.m to quantify burst count and duration based on this threshold value from the Hilbert-transformed, bandpassed amplitude signal. Code for band-passing the signal (eegfilt.m) was developed by Scott Makeig (copyright 1997) and is included in the dependencies folder for convenience. 
