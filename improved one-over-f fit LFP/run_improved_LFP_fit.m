% code for improved 1/f fit of mouse LFP data
% requires the following in the path: the Chronux package and the called
% upon functions
% requires the following in the current folder: the timeseries data spreadsheet

% user-defined parameters
epoch_length = 5;
%length of segments of time series in the input data. Always 5 here.

fitting_range = [1.5 175];

frequency_resolution = 2;
%affects the number of tapers used in the multitaper spectrum. Must be a
%an even number greater than or equal to 2. 

confidence_interval = 99;
%for statistical testing. We use 99 for mice. Can also input 95 here.

%end of user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%begin code
FX = improved_LFP_fit('FX_sample_epoched_V1_LFP.mat',epoch_length,fitting_range,frequency_resolution);
WT = improved_LFP_fit('WT_sample_epoched_V1_LFP.mat',epoch_length,fitting_range,frequency_resolution);

B_L = 1000;  %number of bootstrap repetitions  
if confidence_interval == 95
CIL = 25;
CIU = 975;
elseif confidence_interval == 99
CIL = 5;
CIU = 995;    
end

%plot aboslute spectrum
plot_spectrum_two_groups(1,10*FX.st_spectrum,10*WT.st_spectrum,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%bootstrap test for significance
[x_data,xxline] = bootstrap_two_groups(2,[],B_L,1,FX.spectrum,WT.spectrum,FX.trial_counter,WT.trial_counter,1,FX.all_freqs,CIU,CIL,fitting_range);
xlim([0 fitting_range(2)+fitting_range(1)])
val = -10;
scatter(x_data,xxline*val,1,'r.')
hold off

% plot aperiodic spectrum
plot_spectrum_two_groups(3,10*FX.st_aperiodic,10*WT.st_aperiodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%bootstrap test for significance
[x_data,xxline] = bootstrap_two_groups(4,[],B_L,1,FX.st_aperiodic,WT.st_aperiodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
xlim([0 fitting_range(2)+fitting_range(1)])
val = -0.05;
scatter(x_data,xxline*val,1,'r.')
hold off

% aperiodic parameters
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
[p,~,stats] = ranksum(FX.st_aperiodic(:,1),WT.st_aperiodic(:,1), 'Alpha', 0.05)
effect_size = abs(stats.zval/sqrt(length(FX.st_aperiodic(:,1))+length(WT.st_aperiodic(:,1))))

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
[pp,~,statss] = ranksum(FX.st_exponent,WT.st_exponent, 'Alpha', 0.05)
effect_size = statss.zval/sqrt(length(FX.st_exponent)+length(WT.st_exponent))

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
[pp,~,statss] = ranksum(FX.st_knee_freq,WT.st_knee_freq, 'Alpha', 0.05)
effect_size = statss.zval/sqrt(length(FX.st_knee_freq)+length(WT.st_knee_freq))

% plot periodic spectrum
plot_spectrum_two_groups(8,10*FX.st_periodic,10*WT.st_periodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%boostrap test for significance
[x_data,xxline] = bootstrap_two_groups(9,[],B_L,1,FX.st_periodic,WT.st_periodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
xlim([0 fitting_range(2)+fitting_range(1)])
val = -0.5;
scatter(x_data,xxline*val,1,'r.')
set(gcf,'renderer','painters')

%plot periodic peak
lim = find(FX.fitted_freqs <= 10.5);
plot_spectrum_two_groups(10,10*FX.st_periodic_zoom(:,lim),10*WT.st_periodic_zoom(:,lim),FX.fitted_freqs(lim),0)

%bootstrap test for significance
[x_data,xxline] = bootstrap_two_groups(11,[],B_L,1,FX.st_periodic_zoom(:,lim),WT.st_periodic_zoom(:,lim),[],[],0,FX.fitted_freqs(lim),CIU,CIL,[fitting_range(1) 10.5]);
xlim([0 11])
val = -0.5;
scatter(x_data,xxline*val,1,'r.')
set(gcf,'renderer','painters')

%plot periodic peak(s)
%peak 1a frequency
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
[p,~,stats] = ranksum(FX.st_peak_params(:,1),WT.st_peak_params(:,1))
effect_size = abs(stats.zval/sqrt(length(FX.st_peak_params(:,1))+length(WT.st_peak_params(:,1))))

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
[p,~,stats] = ranksum(FX.st_peak_params(:,2),WT.st_peak_params(:,2))
effect_size = abs(stats.zval/sqrt(length(FX.st_peak_params(:,2))+length(WT.st_peak_params(:,2))))

%peak 1a power
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
[p,~,stats] = ranksum(10*FX.st_peak_params(:,3),10*WT.st_peak_params(:,3))
effect_size = abs(stats.zval/sqrt(length(FX.st_peak_params(:,3))+length(WT.st_peak_params(:,3))))

%peak 1b power
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
[p,~,stats] = ranksum(10*FX.st_peak_params(:,4),10*WT.st_peak_params(:,4))
effect_size = abs(stats.zval/sqrt(length(FX.st_peak_params(:,4))+length(WT.st_peak_params(:,4))))

