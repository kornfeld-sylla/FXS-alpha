% examine burst dynamics in mice
% user must select files they wish to analyze from the list below (lines
% 38-60). The code will automatically run with the sample data in the
% Github. If the user would like to analyze any of the other data, they
% must add it to the current folder from FigShare. Be sure to have the
% correct corresponding file from each genotype list selected (should match
% up in the file name, one FX and one WT).


burst_threshold = 0.9;
%amplitude threshold (90th percentile) for being considered a burst

a_prompt = "Would you like to analyze freely-moving or head-fixed mouse data? Freely-moving = 1, Head-fixed = 2." ;
resp = input(a_prompt);
if resp == 1
    data_type = 'EEG';
else
    data_type = 'LFP';
end

%ask user if they would like to implement zero padding or edge trimming
prompt = "Would you like to implement zero padding? Note that there are no major edge artifacts in the data for this manuscript. Yes = 1, No = 0.";
zero_pad = input(prompt);
if zero_pad ~= 1
prompts = "Would you like to implement edge trimming? Note that there are no major edge artifacts in the data for this manuscript. Yes = 1, No = 0.";
edge_trimming = input(prompts);
end

%ask use if they would like to see a spectogram
promptly = "Would you like to also see a time-frequency analysis of the data? (Requires Chronux toolbox) Yes = 1, No = 0.";
spect = input(promptly);

%begin file selection based on genotype - user must select which pairs of
%files from FigShare they would like to use. The sample files included in
%the Github are automatically used otherwise.
for GG = 1:2
    if GG == 1
        if strcmp(data_type,'EEG') == 1
        load('Supp Fig 7 - FX_continuous_V1_EEG.mat');
        %load('Supp Fig 7 - FX_continuous_V1_LFP_wEEG.mat');
        else
        load('Supp Fig 7 - FX_ADULT_gray_continuous_V1_LFP.mat');  
        %load('Supp Fig 7 - FX_JUV_gray_continuous_V1_LFP.mat');
        %load('Supp Fig 7 - FX_MONITOROFF_continuous_V1_LFP.mat');
        %load('Supp Fig 8 - FX_cko_continuous_V1_LFP.mat');
        %load('Supp Fig 8 - WTcontrol_cko_continuous_V1_LFP.mat'); %compare this to WT_cko_continuous_V1_LFP.mat as well (control vs. control for triple transgenic experiment)
        %load('Supp Fig 9 - PV_pre_CNO_continuous_V1_LFP.mat');
        end
        ts_data = FX_burst_ts;
    else
        if strcmp(data_type,'EEG') == 1
        load('Supp Fig 7 - WT_continuous_V1_EEG.mat');
        %load('Supp Fig 7 - WT_continuous_V1_LFP_wEEG.mat');
        else
        load('Supp Fig 7 - WT_ADULT_gray_continuous_V1_LFP.mat'); 
        %load('Supp Fig 7 - WT_JUV_gray_continuous_V1_LFP.mat');
        %load('Supp Fig 7 - WT_MONITOROFF_continuous_V1_LFP.mat');
        %load('Supp Fig 8 - WT_cko_continuous_V1_LFP.mat');
        %load('Supp Fig 9 - PV_post_CNO_continuous_V1_LFP.mat');
        end
        ts_data = WT_burst_ts;
    end
    
%preallocate results arrays
mean_burst_amplitude = zeros(size(ts_data,1),1);
num_bursts = zeros(size(ts_data,1),1);
mean_burst_length = zeros(size(ts_data,1),1);

fs = 1000; %sampling frequency

%bandpass filter, requires eegfilt 
for a = 1:size(ts_data,1)
    
    sample_ts = ts_data(a,:); 
    %make zero padding to next power of 2
    total_length = nextpow2(length(sample_ts));
    total_num_zeros = (2^total_length) - length(sample_ts);
    num_zeros = total_num_zeros/2;
    padded_sample_ts = padarray(sample_ts,[0 num_zeros],0,'both');

    if strcmp(data_type,'EEG') == 1
        band_edges = [3 8];
        %for freely-moving data, use narrower band edges of 3 and 8 Hz
        %because data has steeper aperioic slope over the range of periodic pk1
    else
        band_edges = [2 10];
        %for head-fixed data, use band edges that capture the range of
        %periodic pk1, 2 and 10 Hz
    end
    
    %band pass either the original or the zero-padded time series
    if zero_pad == 1
    freq_band_temp = eegfilt(padded_sample_ts,fs,band_edges(1),band_edges(2));
    else
    freq_band = eegfilt(sample_ts,fs,band_edges(1),band_edges(2));
    end

    %calculate hilbert transform 
    if zero_pad == 1
        freq_hil_temp = hilbert(freq_band_temp);
        freq_hil = freq_hil_temp(num_zeros+1:end-num_zeros);
    else
        %trim off 2 seconds on either end of continuous segment, if user
        %selected to do this
        if edge_trimming == 1
               freq_hil_temp = hilbert(freq_band);
               trim_length = 2*fs;
               freq_hil = freq_hil_temp(trim_length:end-trim_length-1);
               sample_ts = ts_data(a,trim_length:end-trim_length-1);
        else
        freq_hil = hilbert(freq_band);
        end
    end
    
    %find the amplitude envelope
    freq_burst_anal = abs(freq_hil);
    
    %run burst analysis, either with spectogram displayed or without
    if spect == 1
    [mean_burst_amplitude(a),num_bursts(a),mean_burst_length(a)] = BURST_ANALYSIS_W_SPECGM(freq_burst_anal,burst_threshold,sample_ts,band_edges,fs);
    else
    [mean_burst_amplitude(a),num_bursts(a),mean_burst_length(a)] = BURST_ANALYSIS(freq_burst_anal,burst_threshold);
    end

end

% save values based on group
if GG == 1
    FX_num_bursts = num_bursts;
    FX_mean_burst_length = mean_burst_length;
else
    WT_num_bursts = num_bursts;
    WT_mean_burst_length = mean_burst_length;
end

end

%plot the number of bursts
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
[pp,~,stats] = ranksum(FX_num_bursts,WT_num_bursts)
effect_size = meanEffectSize(FX_num_bursts,WT_num_bursts,Paired=false,Effect="cliff",Alpha=0.01) 


%plot the average burst duration
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
[pp,~,stats] = ranksum(FX_mean_burst_length,WT_mean_burst_length)
effect_size = meanEffectSize(FX_mean_burst_length,WT_mean_burst_length,Paired=false,Effect="cliff",Alpha=0.01) 

