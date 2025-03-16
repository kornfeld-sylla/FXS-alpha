Code for analysis of the 1/f fit of LFP signals from a mouse. The code uses the same math as previously published code from the Voytek Lab at UCSD (SpecParam, https://fooof-tools.github.io/fooof/), namely Lorentzian Functions, but SpecParam has difficulty fitting the signal well at the lowest frequencies and between 10-20 Hz and often estimates underlying 1/f power in these regions. This code tailors the fit specifically to these troublesome regions for an improved 1/f fit.  

This code requires MATLAB, and was developed in version R2018b. This code also calls on functions from the Chronux Toolbox for signal processing. Please download Chronux version 2.12 here: https://chronux.org/

To use this code, simply execute run_improved_LFP_fit.m. To execute, you need to have all the dependencies in the path, including the Chronux code, and the data in the current folder. Run time for this code is difficult to estimate as the user must respond to prompts. The code takes about 5-6 minutes to run, not including added user response time. 

There are several user-defined parameters at the top of rum_improved_LFP_fit.m. The code will execute properly with the parameters as they are, but they are modifiable. The first parameter allows the user to indicate the length of the time series segments you would like to analyze in the frequency domain. This value must be 5 (sec) for the example data in this Github. The next parameter allows the user to indicate the frequency range over which they want to fit the 1/f, we recommend [1.5 175] (default). The accuracy of our fit comes from fitting out to 175 Hz, as there is still a lot of periodic power until about 150 Hz. The following parameter allows the user to indicate the frequency resolution they would like to have for the multi-taper power spectrum, which must be an even number greater than or equal to 2 (recommended, default). Finally, the last parameter controls the confidence interval you want for bootstrapping. Accepted values are 95 and 99. We recommend 99 (default).

Dependencies: This code calls on improved_LFP_fit.m to calculate the power spectra from the given time-series. This code calls upon another script, st_fit_LFP.m, which finds the aperiodic fit for the average power spectrum per subject and also finds the center frequencies and powers of the periodic pk1 subpeaks (a and b). It will attempt to find these peaks through automated methods, but will ask you to confirm each peak by entering yes (1) or no (0). Instructions for this process will pop up in the command window as you run the code. 

Additional dependencies: For plotting power spectra and running bootstrap statistics, the code also calls on plot_spectrum_two_group.m and fg_bootstrap_two_groups.m. 







