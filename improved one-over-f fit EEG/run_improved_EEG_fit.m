% code for improved 1/f fit of EEG data (mouse or human)

% requires the following in the path: the Chronux package and the following 
% called upon functions: improved_EEG_fit.m, st_fit_EEG.m, 
% plot_spectrum_two_groups.m, fg_bootstrap_two_groups.m,
% fg_bootstrap_two_groups_seed.m, and makeseed.m

% requires the following in the current folder: the timeseries data
% spreadsheets and the seeds*.mat files if using pre-made seeds

% user-defined parameters

%which timeseries would you like to analyze? Please write the filenames
%here. 
% For example, to reproduce Fig. 2, call these sheets (.mat files):
 sheet1 = 'Fig 2 - FX_epoched_V1_EEG.mat';
 sheet2 = 'Fig 2 - WT_epoched_V1_EEG.mat';

%To reproduce Supp. Fig 4, add these .mat files from Figshare to the
%current folder, and comment out the previous files and uncomment these files instead:
% sheet1 = 'Supp Fig 4 - FX_epoched_S1_EEG.mat';
% sheet2 = 'Supp Fig 4 - WT_epoched_S1_EEG.mat';

%to match the graphs in the mansucript, please make sure sheet 1 is the FXS data

epoch_length = 5;
%either 2 (sec) for human or 5 (sec) for mice

fitting_range = [2 55];
%in the paper, we used [2 55] for mouse data, [3 55] for human data
%however, the code can handle fitting out to higher frequencies, even if a
%notch filter was used. For example, you could also try [2 90] with the
%sample data.

num_peaks = 1;
%1 peak for mouse data, 2 peaks for human data

frequency_resolution = 2;
%affects the number of tapers used in the multitaper spectrum. please
%select a value greater than or equal to 1. Use a value of 2 to match the manuscript. 

confidence_interval = 99;
%we use 99% for mice, 95% for humans

%end of user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%begin code
%create data structures with absolute, aperiodic, and periodic power
%spectra (and frequency bins) and aperiodic and periodic peak parameters
FX = improved_EEG_fit(sheet1,epoch_length,fitting_range,num_peaks,frequency_resolution);
WT = improved_EEG_fit(sheet2,epoch_length,fitting_range,num_peaks,frequency_resolution);

B_L = 1000;  %number of bootstrap repetitions  
%set bounds of confidence intervals
if confidence_interval == 95
CIL = 25;
CIU = 975;
elseif confidence_interval == 99
CIL = 5;
CIU = 995;    
end

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

%find sample size. Certain statistical tests only run if n >= 10
sample_sizes = [size(FX.st_spectrum,1),size(WT.st_spectrum,1)];

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
val = -5;
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
%give uncorrected p-value and cliff's delta effect size
[uncorrected_p,~,stats] = ranksum(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
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
%give uncorrected p-value and cliff's delta effect size
[uncorrected_pp,~,statss] = ranksum(FX.st_exponent,WT.st_exponent)
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_exponent,WT.st_exponent,Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end

% correct for multiple comparisons across the aperiodic parameters
corrected_aperiodic_pvals = FDR_correct([uncorrected_p,uncorrected_pp])

% plot periodic spectrum
plot_spectrum_two_groups(7,10*FX.st_periodic,10*WT.st_periodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%boostrap test for significance for periodic spectrum and plot
if det_test == 0
[x_data,xxline] = fg_bootstrap_two_groups(8,B_L,FX.st_periodic,WT.st_periodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
else
    if exist('seeds3')
    else
        if sample_sizes(1) == sample_sizes(2)
        seeds3 = makeseed(B_L,max(sample_sizes),2);
        else
        seeds3 = makeseed(B_L,max(sample_sizes),4);
        end
    end
[x_data,xxline] = fg_bootstrap_two_groups_seed(8,seeds3,B_L,FX.st_periodic,WT.st_periodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
end    
xlim([0 fitting_range(2)+fitting_range(1)])
val = -0.5;
scatter(x_data,xxline*val,1,'r.')
ylim([val-0.1, max(ylim)])
set(gcf,'renderer','painters')

%plot parameters for periodic peak(s) as box-plots
%peak 1 center frequency
figure(9)
subplot(1,2,1)
boxplot([FX.st_first_peak(:,1)],{'KO'});
hold on
basis = ones(length(FX.st_first_peak(:,1)),1);
scatter(basis,FX.st_first_peak(:,1));
hold off
subplot(1,2,2)
boxplot([WT.st_first_peak(:,1)],{'WT'});
hold on
basis = ones(length(WT.st_first_peak(:,1)),1);
scatter(basis,WT.st_first_peak(:,1));
hold off
%give uncorrected p-value and cliff's delta effect size
[uncorrected_p1,~,stats1] = ranksum(FX.st_first_peak(:,1),WT.st_first_peak(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_first_peak(:,1),WT.st_first_peak(:,1),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end

%peak 1 maximum power
figure(10)
subplot(1,2,1)
boxplot([10*FX.st_first_peak(:,2)],{'FX'});
hold on
basis = ones(length(FX.st_first_peak(:,2)),1);
scatter(basis,10*FX.st_first_peak(:,2));
hold off
subplot(1,2,2)
boxplot([10*WT.st_first_peak(:,2)],{'WT'});
hold on
basis = ones(length(WT.st_first_peak(:,2)),1);
scatter(basis,10*WT.st_first_peak(:,2));
%give uncorrected p-value and cliff's delta effect size
[uncorrected_p2,~,stats2] = ranksum(FX.st_first_peak(:,2),WT.st_first_peak(:,2))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_first_peak(:,2),WT.st_first_peak(:,2),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end

%plot peak 2 if there is a second peak
%for the data provided on figshare, there will not be any peak 2 plotting because we are only including data from mice
%please contact authors for requsts for human data
if num_peaks > 1
    %peak 2 center frequency
figure(11)
subplot(1,2,1)
boxplot([FX.st_second_peak(:,1)],{'FX'});
hold on
basis = ones(length(FX.st_second_peak(:,1)),1);
scatter(basis,FX.st_second_peak(:,1));
hold off
subplot(1,2,2)
boxplot([WT.st_second_peak(:,1)],{'WT'});
hold on
basis = ones(length(WT.st_second_peak(:,1)),1);
scatter(basis,WT.st_second_peak(:,1));
[uncorrected_p3,~,stats3] = ranksum(FX.st_second_peak(:,1),WT.st_second_peak(:,1))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_second_peak(:,1),WT.st_second_peak(:,1),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end

%peak 2 maximum power
figure(12)
subplot(1,2,1)
boxplot([10*FX.st_second_peak(:,2)],{'FX'});
hold on
basis = ones(length(FX.st_second_peak(:,2)),1);
scatter(basis,10*FX.st_second_peak(:,2));
hold off
subplot(1,2,2)
boxplot([10*WT.st_second_peak(:,2)],{'WT'});
hold on
basis = ones(length(WT.st_second_peak(:,2)),1);
scatter(basis,10*WT.st_second_peak(:,2));
hold off
[uncorrected_p4,~,stats4] = ranksum(FX.st_second_peak(:,2),WT.st_second_peak(:,2))
if sample_sizes(1) && sample_sizes(2) >= 10
effect_size = meanEffectSize(FX.st_second_peak(:,2),WT.st_second_peak(:,2),Paired=false,Effect="cliff",Alpha = (100-confidence_interval)/100) 
end

%correct for multiple comparisons across periodic parameters
corrected_periodic_pvals = FDR_correct([uncorrected_p1,uncorrected_p2,uncorrected_p3,uncorrected_p4])

%test for a correlation between the maximum power of the two peaks
figure(13)
scatter((FX.st_first_peak(:,2)),(FX.st_second_peak(:,2)),30,'r')
hold on
scatter((WT.st_first_peak(:,2)),(WT.st_second_peak(:,2)),30,'k')
xlim([0 1.5])
ylim([0 1]);
FX_data_pk1 = FX.st_first_peak(:,2);
FX_data_pk2 = FX.st_second_peak(:,2);
WT_data_pk1 = WT.st_first_peak(:,2);
WT_data_pk2 = WT.st_second_peak(:,2);

p = polyfit(FX_data_pk1,FX_data_pk2,1);
yfit = p(1)*FX_data_pk1 + p(2);
yresid = FX_data_pk2 - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(FX_data_pk2)-1)*var(FX_data_pk2);
rsq1 = 1-SSresid/SStotal
ef_size1 = sqrt(rsq1)
plot(FX_data_pk1,yfit,'r-.')

p = polyfit(WT_data_pk1,WT_data_pk2,1);
yfit = p(1)*WT_data_pk1 + p(2);
yresid = WT_data_pk2 - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(WT_data_pk2)-1)*var(WT_data_pk2);
rsq2 = 1-SSresid/SStotal
ef_size2 = sqrt(rsq2)
plot(WT_data_pk1,yfit,'k-.')
hold off
else

% if there's no peak 2, just correct for multiple comparisons for peak 1
corrected_periodic_pvals = FDR_correct([uncorrected_p1,uncorrected_p2])

end

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
