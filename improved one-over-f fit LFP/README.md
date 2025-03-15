Code for analysis of the 1/f fit of LFP signals from a mouse. The code uses the same math as SpecParam, namely Lorentzian Functions, but SpecParam has difficulty fitting the signal well at the lowest frequencies and between 10-20 Hz and often estimates underlying 1/f power in these regions. This code tailors the fit specifically to these troublesome regions for an improved 1/f fit.  

This code requires MATLAB, and was developed in version R2018b. This code also calls on functions from the Chronux Toolbox for signal processing. Please download Chronux version 2.12 here: https://chronux.org/

To use this code, execute run_improved_LFP_fit.m with the correct user-defined parameters. Indicate the length of the time series segments you would like to analyze in the frequency domain, the default is 5s for the example data. Please also indicate the frequency range over which you want to fit the 1/f, we recommend [1.5 175] (default). The accuracy of our fit comes from fitting out to 175 Hz, as there is still a lot of periodic power until about 150 Hz. Please also indicate the frequency resolution you would like to have for the multi-taper power spectrum, which must be an even number greater than or equal to 2 (recommended). Finally, list the confidence interval you want for bootstrapping. Accepted values are 95 and 99. We recommend 99 (default).

Dependencies: This code calls on improved_LFP_fit.m to calculate the power spectra from the given time-series. This code calls upon another script, st_fit_LFP.m, which finds the aperiodic fit for the average power spectrum per subject and also finds the center frequencies and powers of the periodic pk1 subpeaks (a and b). It will attempt to find these peaks through automated methods, but will ask you to confirm each peak by entering yes (1) or no (0). Instructions for this process will pop up in the command window as you run the code. 

Additional dependencies: For plotting power spectra and running bootstrap statistics, the code also calls on plot_spectrum_two_group.m and bootstrap_two_groups.m. 







