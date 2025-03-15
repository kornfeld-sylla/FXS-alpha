% code for improved 1/f fit of EEG data (mouse or human)
% requires the following in the path: the Chronux package and the following 
% called upon functions: improved_EEG_fit.m, st_fit_EEG.m, 
% plot_spectrum_two_groups.m, and bootstrap_two_groups.m
% requires the following in the current folder: the timeseries data
% spreadsheets

% user-defined parameters
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
%select a value greater than or equal to 1. 

confidence_interval = 99;
%we use 99% for mice, 95% for humans

%end of user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%begin code
%begin with KO data
FX = improved_EEG_fit('FX_epoched_V1_EEG.mat',epoch_length,fitting_range,num_peaks,frequency_resolution);
WT = improved_EEG_fit('WT_epoched_V1_EEG.mat',epoch_length,fitting_range,num_peaks,frequency_resolution);

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
[x_data,xxline] = fg_bootstrap_two_groups(2,B_L,FX.spectrum,WT.spectrum,FX.trial_counter,WT.trial_counter,1,FX.all_freqs,CIU,CIL,fitting_range);
xlim([0 fitting_range(2)+fitting_range(1)])
val = -10;
scatter(x_data,xxline*val,1,'r.')
hold off

% plot aperiodic spectrum
plot_spectrum_two_groups(3,10*FX.st_aperiodic,10*WT.st_aperiodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%bootstrap test for significance
[x_data,xxline] = fg_bootstrap_two_groups(4,B_L,FX.st_aperiodic,WT.st_aperiodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
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

% plot periodic spectrum
plot_spectrum_two_groups(7,10*FX.st_periodic,10*WT.st_periodic,FX.fitted_freqs,0)
xlim([0 fitting_range(2)+fitting_range(1)])

%boostrap test for significance
[x_data,xxline] = fg_bootstrap_two_groups(8,B_L,FX.st_periodic,WT.st_periodic,[],[],0,FX.fitted_freqs,CIU,CIL,fitting_range);
xlim([0 fitting_range(2)+fitting_range(1)])
val = -0.5;
scatter(x_data,xxline*val,1,'r.')
set(gcf,'renderer','painters')

%plot periodic peak(s)
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
[p,~,stats] = ranksum(FX.st_first_peak(:,1),WT.st_first_peak(:,1))
effect_size = abs(stats.zval/sqrt(length(FX.st_first_peak(:,1))+length(WT.st_first_peak(:,1))))

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
[pp,~,stats] = ranksum(FX.st_first_peak(:,2),WT.st_first_peak(:,2))
effect_size = abs(stats.zval/sqrt(length(FX.st_first_peak(:,2))+length(WT.st_first_peak(:,2))))

%plot peak 2 if there is a second peak
if num_peaks > 1
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
[p,~,stats] = ranksum(FX.st_second_peak(:,1),WT.st_second_peak(:,1))
effect_size = abs(stats.zval/sqrt(length(FX.st_second_peak(:,1))+length(WT.st_second_peak(:,1))))

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
[pp,~,statss] = ranksum(FX.st_second_peak(:,2),WT.st_second_peak(:,2))
effect_size = abs(statss.zval/sqrt(length(FX.st_second_peak(:,2))+length(WT.st_second_peak(:,2))))

%test for a correlation between the power of the two peaks
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
end