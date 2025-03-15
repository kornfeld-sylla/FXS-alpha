function [fitted_data] =improved_LFP_fit(name_of_data_file,epoch_length,fitting_range,frequency_resolution)
%Find power spectra for each epoch and run 1/f fit for average spectra per animal
%   Output it a data structure with all spectra and fit information
%   Input is name of data file to analyze, length of epochs in s, rage over
%   which to do 1/f fit, and frequency resolution for multitaper spectra

load(name_of_data_file);

if name_of_data_file(1:2) == 'FX'
    ts_data = FX_fit_ts_new;
else
    ts_data = WT_fit_ts_new;
end

%data will be a cell
% the length of the cell is the number of subjects
%dimension 1 of each cell = number of epochs 
%dimension 2 of each cell = length of epoch*sampling frequency

fs = size(ts_data{1},2)/epoch_length;
num_subj = length(ts_data);

%preallocate results arrays
pre_spectrum = [];
trial_counter = zeros(num_subj,1);

le = (epoch_length*(fitting_range(2)-fitting_range(1)))+1;
le = floor(le);

st_spec = zeros(num_subj,le);
st_aperiodic = zeros(num_subj,le);
st_intercept = zeros(num_subj,1);
st_exponent = zeros(num_subj,1);
st_knee_freq = zeros(num_subj,1);
st_periodic = zeros(num_subj,le);
st_periodic_zoom = zeros(num_subj,le);
st_peak_params = zeros(num_subj,4);    


%set multitaper parameters. You will need the Chronux package to run
tw_2 = (epoch_length*frequency_resolution)/2;
params2.tapers = [tw_2, (tw_2*2)-1];
params2.pad = -1;
params2.Fs = fs;

%set an additional parameter set with higher resolution for peak analysis
tw_3 = (epoch_length*(frequency_resolution/2))/2;
params3.tapers = [tw_3, (tw_3*2)-1];
params3.pad = -1;
params3.Fs = fs;

%find the multitaper spectrum for each subject. Save each epoched spectra
%for bootstrapping, and find the average spectra per subject for 1/f
%fitting.
for a = 1:num_subj
    num_epochs = size(ts_data{a},1);
    trial_counter(a) = num_epochs;
    temp_for_fit = [];
    temp_for_peak_fit = [];
    for b = 1:num_epochs
        sts = ts_data{a}(b,:);
        [spec,all_freqs]=mtspectrumc(sts,params2);
        spec = spec';
        pre_spectrum = [pre_spectrum; spec];
        temp_for_fit = [temp_for_fit; spec];
        [zoom_spec,~]=mtspectrumc(sts,params3);
        zoom_spec = zoom_spec';
        temp_for_peak_fit = [temp_for_peak_fit; zoom_spec];
    end
    fit_spec = mean(temp_for_fit,1);
    zoom_fit_spec = mean(temp_for_peak_fit,1);

%fit aperiodic spectra with lower frequency resolution (peak_analysis input set to 0)  
[st_intercept(a),st_exponent(a),st_knee_freq(a),st_spec(a,:),st_freqs,st_aperiodic(a,:),st_periodic(a,:),~] = st_fit_LFP(all_freqs, fit_spec, fitting_range, 0);
%find peak params with higher freqency resolution (peak_analysis input set to 1)
[~,~,~,~,~,~,st_periodic_zoom(a,:),st_peak_params(a,:)] = st_fit_LFP(all_freqs, zoom_fit_spec, fitting_range, 1);
   
end

%save to data structure
 fitted_data.spectrum = pre_spectrum;
 fitted_data.trial_counter = trial_counter;
 
 fitted_data.st_spectrum = st_spec; 
 fitted_data.st_aperiodic = st_aperiodic;
 fitted_data.st_intercept = st_intercept; 
 fitted_data.st_exponent = st_exponent;
 fitted_data.st_knee_freq = st_knee_freq;
 fitted_data.st_periodic = st_periodic;
 fitted_data.st_periodic_zoom = st_periodic_zoom;
 fitted_data.st_peak_params = st_peak_params;
 fitted_data.all_freqs = all_freqs;
 fitted_data.fitted_freqs = st_freqs;
end

