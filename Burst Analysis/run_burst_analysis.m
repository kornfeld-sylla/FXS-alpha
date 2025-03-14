% examine burst dynamics in mice

%user defined inputs
data_type = 'EEG';
%input either EEG or LFP here depending on data type

burst_threshold = 0.9;
%amplitude threshold (90th percentile) for being considered a burst

sampling_frequency = 1000;
% in Hz

%begin code
for GG = 1:2
    if GG == 1
        load('FX_continuous_S1_EEG.mat');
        ts_data = FX_burst_ts;
    else
        load('WT_continuous_S1_EEG.mat');
        ts_data = WT_burst_ts;
    end
    
%preallocate results arrays
mean_burst_power = zeros(size(ts_data,1),1);
num_bursts = zeros(size(ts_data,1),1);
mean_burst_length = zeros(size(ts_data,1),1);

fs = sampling_frequency;

%bandpass filter, requires eegfilt 
for a = 1:size(ts_data,1)
    if strcmp(data_type,'EEG') == 1
    freq_band = eegfilt(ts_data(a,:),fs,3,8);
    else
    freq_band = eegfilt(ts_data(a,:),fs,2,10);
    end
    freq_hil = hilbert(freq_band);
    freq_burst_anal = abs(freq_hil);
    
    [mean_burst_power(a),num_bursts(a),mean_burst_length(a)] = BURST_ANALYSIS(freq_burst_anal,burst_threshold);

end

if GG == 1
    FX_num_bursts = num_bursts;
    FX_mean_burst_length = mean_burst_length;
else
    WT_num_bursts = num_bursts;
    WT_mean_burst_length = mean_burst_length;
end

end

figure(1)
subplot(1,2,1)
boxplot([FX_num_bursts],{'FX'});
hold on
gosh = ones(length(FX_num_bursts),1);
scatter(gosh,FX_num_bursts);
hold off
subplot(1,2,2)
boxplot([WT_num_bursts],{'WT'});
hold on
gosh = ones(length(WT_num_bursts),1);
scatter(gosh,WT_num_bursts);
hold off
[pp,h,stats] = ranksum(FX_num_bursts,WT_num_bursts)
effect_size = abs(stats.zval/sqrt(length(FX_num_bursts)+length(WT_num_bursts)))

figure(2)
subplot(1,2,1)
boxplot([FX_mean_burst_length],{'FX'});
hold on
gosh = ones(length(FX_mean_burst_length),1);
scatter(gosh,FX_mean_burst_length);
hold off
subplot(1,2,2)
boxplot([WT_mean_burst_length],{'WT'});
hold on
gosh = ones(length(WT_mean_burst_length),1);
scatter(gosh,WT_mean_burst_length);
hold off
[pp,h,stats] = ranksum(FX_mean_burst_length,WT_mean_burst_length)
effect_size = abs(stats.zval/sqrt(length(FX_mean_burst_length)+length(WT_mean_burst_length)))
        