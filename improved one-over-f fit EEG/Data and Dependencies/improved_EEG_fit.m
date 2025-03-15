function [fitted_data] =improved_EEG_fit(name_of_data_file,epoch_length,fitting_range,num_peaks,frequency_resolution)
%Find power spectra for each epoch and run 1/f fit for average spectra per animal
%   Output it a data structure with all spectra and fit information
%   Input is name of data file to analyze, length of epochs in s, rage over
%   which to do 1/f fit, number of periodic peaks to analyze, and frequency
%   resolution for multitaper spectra

load(name_of_data_file);

if strcmp(name_of_data_file,'FX_epoched_V1_EEG.mat')
    ts_data = FX_fit_ts;
else
    ts_data = WT_fit_ts;
end

%data will be a cell
%the length of the cell is the number of subjects
%dimension 1 of each cell = number of epochs 
%dimension 2 of each cell = length of epoch*sampling frequency

fs = size(ts_data{1},2)/epoch_length;
num_subj = length(ts_data);

%preallocate results arrays
le = (epoch_length*(fitting_range(2)-fitting_range(1)))+1;
le = floor(le);

st_spec = zeros(num_subj,le);
st_aperiodic = zeros(num_subj,le);
st_intercept = zeros(num_subj,1);
st_exponent = zeros(num_subj,1);
st_periodic = zeros(num_subj,le);
st_first_peak_params = zeros(num_subj,2);
if num_peaks > 1
st_second_peak_params = zeros(num_subj,2);
end

%if all subjects have the same number of epochs, organize into a 3d array
%for processing speed
trial_counter = zeros(num_subj,1);
truth = zeros(1,num_subj);
for aa = 1:num_subj
truth(aa) = size(ts_data{aa},1);
end
compare = ones(1,num_subj)*truth(1);
truth_test = truth == compare;
if truth_test == true(1,num_subj)
pre_spectrum = zeros(num_subj,truth(1),((epoch_length*fs)/2)+1);    
else
pre_spectrum = [];
end

%set multitaper parameters. You will need the Chronux package to run
tw_2 = (epoch_length*frequency_resolution)/2;
params2.tapers = [tw_2, (tw_2*2)-1];
params2.pad = -1;
params2.Fs = fs;

%calculate spectra for each epoch
for a = 1:num_subj
    num_epochs = size(ts_data{a},1);
    trial_counter(a) = num_epochs;
    temp_for_fit = [];
    for b = 1:num_epochs
        sts = ts_data{a}(b,:);
        [spec,all_freqs]=mtspectrumc(sts,params2);
        spec = spec';
        if truth_test == 1
            pre_spectrum(a,b,:) = spec;
        else
        pre_spectrum = [pre_spectrum; spec];
        end
        temp_for_fit = [temp_for_fit; spec];
    end
    fit_spec = mean(temp_for_fit,1);

%run code for fitting 1/f and periodic peaks    
if num_peaks > 1    
[st_intercept(a),st_exponent(a),st_spec(a,:),st_freqs,st_aperiodic(a,:),st_periodic(a,:),st_first_peak_params(a,:),st_second_peak_params(a,:)] = st_fit_EEG(all_freqs,fit_spec,fitting_range,num_peaks);
else
[st_intercept(a),st_exponent(a),st_spec(a,:),st_freqs,st_aperiodic(a,:),st_periodic(a,:),st_first_peak_params(a,:),~] = st_fit_EEG(all_freqs,fit_spec,fitting_range,num_peaks);
end    
end

%save to data structure
 fitted_data.spectrum = pre_spectrum;
 fitted_data.trial_counter = trial_counter;
 fitted_data.st_spectrum = st_spec; 
 fitted_data.st_aperiodic = st_aperiodic;
 fitted_data.st_intercept = st_intercept; 
 fitted_data.st_exponent = st_exponent;
 fitted_data.st_periodic = st_periodic;
 fitted_data.st_first_peak = st_first_peak_params;
 if num_peaks > 1
 fitted_data.st_second_peak = st_second_peak_params;
 end
 fitted_data.all_freqs = all_freqs;
 fitted_data.fitted_freqs = st_freqs;
end

