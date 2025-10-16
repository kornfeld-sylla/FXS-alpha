% code for improved 1/f fit of mouse LFP data
% requires the following in the path: the Chronux package and the following 
% called upon functions: improved_LFP_fit.m,st_fit_LFP.m,
% plot_spectrum_two_groups.m, fg_bootstrap_two_groups.m 
% fg_bootstrap_two_groups_seed.m, and makeseed.m
% requires the following in the current folder: the timeseries data
% spreadsheets and the seeds*.mat files if using pre-made seeds

% user-defined parameters

%which timeseries would you like to analyze? Please write the filenames here. 
% For example, to reproduce Fig. 3a-g, call these sheets (.mat files):
sheet1 = 'Fig 3 - FX_ADULT_gray_epoched_V1_LFP.mat';
sheet2 = 'Fig 3 - WT_ADULT_gray_epoched_V1_LFP.mat';

%to match the graphs in the mansucript, please make sure sheet 1 is the FXS
%data. To run any of the other datasets in the following list, please use
%the files assigned to sheet1 and sheet2 as written to match exactly the
%published approach.

%For Fig. 3h-n, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 3 - FX_JUV_gray_epoched_V1_LFP.mat';
%sheet2 = 'Fig 3 - WT_JUV_gray_epoched_V1_LFP.mat';

%For Fig. 4, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 4 - FX_cko_gray_epoched_V1_LFP.mat';
%sheet2 = 'Fig 4 - WT_cko_gray_epoched_V1_LFP.mat';

%For Fig. 6a-j, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 6 - PV_pre_CNO_epoched_V1_LFP.mat';
%sheet2 = 'Fig 6 - PV_post_CNO_epoched_V1_LFP.mat';

%For Fig. 6k-t, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 6 - SOM_epoched_V1_LFP.mat';
%sheet2 = 'Fig 6 - Control_SOM_epoched_V1_LFP.mat';

%For Fig. 7a-j, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 7 - WT_Arbac_halfmpk_epoched_V1_LFP.mat';
%sheet2 = 'Fig 7 - WT_Saline_forhalfmpk_epoched_V1_LFP.mat';

%For Fig. 7k-t, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 7 - WT_Arbac_1mpk_epoched_V1_LFP.mat';
%sheet2 = 'Fig 7 - WT_Saline_for1mpk_epoched_V1_LFP.mat';

%For Fig. 8a-e, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 8 - FX_Arbac_halfmpk_epoched_V1_LFP.mat';
%sheet2 = 'Fig 8 - FX_Saline_forhalfmpk_epoched_V1_LFP.mat';

%For Fig. 8f-j, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Fig 8 - FX_Arbac_1mpk_epoched_V1_LFP.mat';
%sheet2 = 'Fig 8 - FX_Saline_for1mpk_epoched_V1_LFP.mat';

%For Supp Fig. 5, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Supp Fig 5 - FX_epoched_V1_LFP.mat';
%sheet2 = 'Supp Fig 5 - WT_epoched_V1_LFP.mat';

%For Supp Fig. 8, add these sheets (.mat files) from figshare to the current folder:
%To compare EMX1-FMR1-KO v. WT (black screen):
%sheet1 = 'Supp Fig 5 - FX_cko_epoched_V1_LFP.mat';
%sheet2 = 'Supp Fig 5 - WT_cko_epoched_V1_LFP.mat';
%Or, to compare Control (WT/WT) v. WT (black screen):
%sheet1 = 'Supp Fig 5 - WTcontrol_cko_epoched_V1_LFP.mat';
%sheet2 = 'Supp Fig 5 - WT_cko_epoched_V1_LFP.mat';

%For Supp Fig. 10, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Supp Fig 10 - OptoSOM_epoched_V1_LFP.mat';
%sheet2 = 'Supp Fig 10 - Control_OptoSOM_epoched_V1_LFP.mat';

%For Supp Fig. 12a-f, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Supp Fig 12 - WT_CTEP_2mpk_epoched_V1_LFP.mat';
%sheet2 = 'Supp Fig 12 - WT_Saline_for2mpk_epoched_V1_LFP.mat';

%For Supp Fig. 12g-l, add these sheets (.mat files) from figshare to the current folder: 
%sheet1 = 'Supp Fig 12 - FX_CTEP_2mpk_epoched_V1_LFP.mat';
%sheet2 = 'Supp Fig 12 - FX_Saline_for2mpk_epoched_V1_LFP.mat';

%continue user-definied parameters here

epoch_length = 5;
%length of segments of time series in the input data. Always 5 here.

fitting_range = [1.5 175];
%frequency range over which to fit the 1/f curve

frequency_resolution = 2;
%affects the number of tapers used in the multitaper spectrum. Must be a
%an even number greater than or equal to 2. 

confidence_interval = 99;
%for statistical testing. We use 99 for mice. Can also input 95 here.

%end of user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin code
FX = improved_LFP_fit(sheet1,epoch_length,fitting_range,frequency_resolution);
WT = improved_LFP_fit(sheet2,epoch_length,fitting_range,frequency_resolution);


B_L = 1000;  %number of bootstrap repetitions  
if confidence_interval == 95
CIL = 25;
CIU = 975;
elseif confidence_interval == 99
CIL = 5;
CIU = 995;    
end

%find sample size. Certain statistical tests only run if n >= 10
sample_sizes = [size(FX.st_spectrum,1),size(WT.st_spectrum,1)];

%ask user if samples should be treated as independent (different genotypes)
%or paired (same animal pre/post drug treatment)
prompt = "Are your two data sets independent (different genotypes) or paired (same animal pre/post drug treatment)? Paired = 1, Independent = 0.";
paired_test = input(prompt);

%set bootstrap mode to deterministic or random
prompt = "Would you like your bootstrap results to be random or deterministic? Deterministic = 1, Random = 0.";
det_test = input(prompt);

%if deterministic, load pre-existing seeds if they are in the current
%folder. If the user answers yes and the files are not in the current
%folder, the code will not work and flag an error
if det_test == 1
    promptly = "Would you like to use pre-existing seeds? You MUST have the three seeds*.m files in your current folder, or the code will not run. Yes = 1, No = 0.";
    seed_test = input(promptly);
    if seed_test == 1
        if isfile('seeds1.mat') == 0 || isfile('seeds2.mat') == 0 || isfile('seeds3.mat') == 0
            error('Error: you do not have the three required pre-existing seed .mat files in your current folder');
        else
            load('seeds1.mat'); load('seeds2.mat'); load('seeds3.mat');
        end
    end
end

%plot aboslute spectrum
plot_spectrum_two_groups(1,10*FX.st_spectrum,10*WT.st_spectrum,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%bootstrap test for significance for absolute spectrum and plot
if det_test == 0
[x_data,xxline] = fg_bootstrap_two_groups(2,B_L,FX.spectrum,WT.spectrum,FX.trial_counter,WT.trial_counter,1,FX.all_freqs,CIU,CIL,fitting_range);
else    
    if exist('seeds1')
    else
        if size(FX.spectrum,1) == size(WT.spectrum,1)
            seeds1 = makeseed(B_L,max(sample_sizes),2);
        else
            seeds1 = makeseed(B_L,max(sample_sizes),4);
        end
    end
    [x_data,xxline] = fg_bootstrap_two_groups_seed(2,seeds1,B_L,FX.spectrum,WT.spectrum,FX.trial_counter,WT.trial_counter,1,FX.all_freqs,CIU,CIL,fitting_range);
end
xlim([0 fitting_range(2)+fitting_range(1)])
val = -10;
scatter(x_data,xxline*val,1,'r.')
ylim([val-1, max(ylim)])
hold off

% plot aperiodic spectrum
plot_spectrum_two_groups(3,10*FX.st_aperiodic,10*WT.st_aperiodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%bootstrap test for significance for aperiodic spectrum and plot
if det_test == 0
[x_data,xxline] = fg_bootstrap_two_groups(4,B_L,FX.st_aperiodic,WT.st_aperiodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
else
    if exist('seeds2')
    else
        if sample_sizes(1) == sample_sizes(2)
        seeds2 = makeseed(B_L,max(sample_sizes),2);
        else
        seeds2 = makeseed(B_L,max(sample_sizes),4);
        end
    end    
[x_data,xxline] = fg_bootstrap_two_groups_seed(4,seeds2,B_L,FX.st_aperiodic,WT.st_aperiodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
end
xlim([0 fitting_range(2)+fitting_range(1)])
val = -0.25;
scatter(x_data,xxline*val,1,'r.')
ylim([val-0.05, max(ylim)])
hold off

% plot aperiodic parameters in box-plots
figure(5)
subplot(1,2,1)
boxplot([10*FX.st_aperiodic(:,1)],{'FX'});
hold on
basis = ones(length(FX.st_intercept),1);
scatter(basis,10*FX.st_aperiodic(:,1));
hold off
subplot(1,2,2)
boxplot([10*WT.st_aperiodic(:,1)],{'WT'});
hold on
basis = ones(length(WT.st_intercept),1);
scatter(basis,10*WT.st_aperiodic(:,1));
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_p,~,stats] = signrank(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1),Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else
[uncorrected_p,~,stats] = ranksum(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

figure(6)
subplot(1,2,1)
boxplot([FX.st_exponent],{'FX'});
hold on
basis = ones(length(FX.st_exponent),1);
scatter(basis,FX.st_exponent);
hold off
subplot(1,2,2)
boxplot([WT.st_exponent],{'WT'});
hold on
basis = ones(length(WT.st_exponent),1);
scatter(basis,WT.st_exponent);
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_pp,~,statss] = signrank(FX.st_exponent,WT.st_exponent)
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_exponent,WT.st_exponent,Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else
[uncorrected_pp,~,statss] = ranksum(FX.st_exponent,WT.st_exponent)
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_exponent,WT.st_exponent,Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

figure(7)
subplot(1,2,1)
boxplot([FX.st_knee_freq],{'FX'});
hold on
basis = ones(length(FX.st_knee_freq),1);
scatter(basis,FX.st_knee_freq);
hold off
subplot(1,2,2)
boxplot([WT.st_knee_freq],{'WT'});
hold on
basis = ones(length(WT.st_knee_freq),1);
scatter(basis,WT.st_knee_freq);
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_ppp,~,statsss] = signrank(FX.st_knee_freq,WT.st_knee_freq)
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_knee_freq,WT.st_knee_freq,Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else    
[uncorrected_ppp,~,statsss] = ranksum(FX.st_knee_freq,WT.st_knee_freq)
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_knee_freq,WT.st_knee_freq,Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

% correct for multiple comparisons across the aperiodic parameters
corrected_aperiodic_pvals = FDR_correct([uncorrected_p,uncorrected_pp,uncorrected_ppp])

% plot periodic spectrum
plot_spectrum_two_groups(8,10*FX.st_periodic,10*WT.st_periodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%boostrap test for significance for periodic spectrum and plot
if det_test == 0
[x_data,xxline] = fg_bootstrap_two_groups(9,B_L,FX.st_periodic,WT.st_periodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
else
    if exist('seeds3')
    else
        if sample_sizes(1) == sample_sizes(2)
        seeds3 = makeseed(B_L,max(sample_sizes),2);
        else
        seeds3 = makeseed(B_L,max(sample_sizes),4);
        end
    end
[x_data,xxline] = fg_bootstrap_two_groups_seed(9,seeds3,B_L,FX.st_periodic,WT.st_periodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
end
xlim([0 fitting_range(2)+fitting_range(1)])
val = -0.5;
scatter(x_data,xxline*val,1,'r.')
ylim([val-0.1, max(ylim)])
set(gcf,'renderer','painters')

%plot periodic peak with 1 Hz multi-taper frequency resolution
lim = find(FX.fitted_freqs <= 10.5);
plot_spectrum_two_groups(10,10*FX.st_periodic_zoom(:,lim),10*WT.st_periodic_zoom(:,lim),FX.fitted_freqs(lim),0)

%boostrap test for significance for periodic spectrum with higher multi-taper frequency resolution and plot
if det_test == 0
[x_data,xxline] = fg_bootstrap_two_groups(11,B_L,FX.st_periodic_zoom(:,lim),WT.st_periodic_zoom(:,lim),[],[],0,FX.fitted_freqs(lim),CIU,CIL,[fitting_range(1) 10.5]);
else
[x_data,xxline] = fg_bootstrap_two_groups_seed(11,seeds3,B_L,FX.st_periodic_zoom(:,lim),WT.st_periodic_zoom(:,lim),[],[],0,FX.fitted_freqs(lim),CIU,CIL,[fitting_range(1) 10.5]);    
end    
xlim([0 11])
val = -0.5;
scatter(x_data,xxline*val,1,'r.')
ylim([val-0.1, max(ylim)])
set(gcf,'renderer','painters')

%plot parameters for periodic peak(s) as box-plots
%peak 1a center frequency
figure(12)
subplot(1,2,1)
boxplot([FX.st_peak_params(:,1)],{'KO'});
hold on
basis = ones(length(FX.st_peak_params(:,1)),1);
scatter(basis,FX.st_peak_params(:,1));
hold off
subplot(1,2,2)
boxplot([WT.st_peak_params(:,1)],{'WT'});
hold on
basis = ones(length(WT.st_peak_params(:,1)),1);
scatter(basis,WT.st_peak_params(:,1));
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_p1,~,stats1] = signrank(FX.st_peak_params(:,1),WT.st_peak_params(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_peak_params(:,1),WT.st_peak_params(:,1),Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else    
[uncorrected_p1,~,stats1] = ranksum(FX.st_peak_params(:,1),WT.st_peak_params(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_peak_params(:,1),WT.st_peak_params(:,1),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

%peak 1b frequency
figure(13)
subplot(1,2,1)
boxplot([FX.st_peak_params(:,2)],{'KO'});
hold on
basis = ones(length(FX.st_peak_params(:,2)),1);
scatter(basis,FX.st_peak_params(:,2));
hold off
subplot(1,2,2)
boxplot([WT.st_peak_params(:,2)],{'WT'});
hold on
basis = ones(length(WT.st_peak_params(:,2)),1);
scatter(basis,WT.st_peak_params(:,2));
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_p2,~,stats2] = signrank(FX.st_peak_params(:,2),WT.st_peak_params(:,2))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_peak_params(:,2),WT.st_peak_params(:,2),Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else    
[uncorrected_p2,~,stats2] = ranksum(FX.st_peak_params(:,2),WT.st_peak_params(:,2))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_peak_params(:,2),WT.st_peak_params(:,2),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

%peak 1a maximum power
figure(14)
subplot(1,2,1)
boxplot([10*FX.st_peak_params(:,3)],{'KO'});
hold on
basis = ones(length(FX.st_peak_params(:,3)),1);
scatter(basis,10*FX.st_peak_params(:,3));
hold off
subplot(1,2,2)
boxplot([10*WT.st_peak_params(:,3)],{'WT'});
hold on
basis = ones(length(WT.st_peak_params(:,3)),1);
scatter(basis,10*WT.st_peak_params(:,3));
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_p3,~,stats3] = signrank(10*FX.st_peak_params(:,3),10*WT.st_peak_params(:,3))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(10*FX.st_peak_params(:,3),10*WT.st_peak_params(:,3),Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else    
[uncorrected_p3,~,stats3] = ranksum(10*FX.st_peak_params(:,3),10*WT.st_peak_params(:,3))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(10*FX.st_peak_params(:,3),10*WT.st_peak_params(:,3),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

%peak 1b maximum power
figure(15)
subplot(1,2,1)
boxplot([10*FX.st_peak_params(:,4)],{'KO'});
hold on
basis = ones(length(FX.st_peak_params(:,4)),1);
scatter(basis,10*FX.st_peak_params(:,4));
hold off
subplot(1,2,2)
boxplot([10*WT.st_peak_params(:,4)],{'WT'});
hold on
basis = ones(length(WT.st_peak_params(:,4)),1);
scatter(basis,10*WT.st_peak_params(:,4));
hold off
%give uncorrected p-value and cliff's delta effect size for either matched
%(signrank) or independent (ranksum) samples
if paired_test == 1
[uncorrected_p4,~,stats4] = signrank(10*FX.st_peak_params(:,4),10*WT.st_peak_params(:,4))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(10*FX.st_peak_params(:,4),10*WT.st_peak_params(:,4),Paired=true,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
else
[uncorrected_p4,~,stats4] = ranksum(10*FX.st_peak_params(:,4),10*WT.st_peak_params(:,4))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(10*FX.st_peak_params(:,4),10*WT.st_peak_params(:,4),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end
end

%correct for multiple comparisons across periodic parameters
corrected_periodic_pvals = FDR_correct([uncorrected_p1,uncorrected_p2,uncorrected_p3,uncorrected_p4])

% if the user generated their own seeds, give them the option to save to
% reproduce the results across MATLAB sessions
if det_test == 1 && seed_test == 0 
final_prompt = "Would you like to save the seeds you have generated? Yes = 1, No = 0.";
    save_test = input(final_prompt);
    if save_test == 1
        save('seeds1','seeds1','-v7.3');
        save('seeds2','seeds2','-v7.3');
        save('seeds3','seeds3','-v7.3');
    end
end
