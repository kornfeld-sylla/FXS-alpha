function [offset,chai,log_spec,freqs,aperiodic,periodic,first_peak_params,second_peak_params] = st_fit_EEG(frequencies, spectrum, fitting_range, num_peaks)

%Function to fit 1/f to EEG data. 
%   Aperiodic component has 2 terms: offset and slope (chai)
%   output gives the absolute, aperiodic, and periodic spectra in log scale
%   with frequencies for x-axis
%   Output also gives center frequencies and powers of periodic peak 1,
%   and, if applicable, periodic peak 2. 

%   Input parameters are the mean absolute spectrum (NOT in log scale), the
%   corresponding frequencies, the range over which to fit the 1/f, and
%   how many periodic peaks you want to fit (1 or 2).


% first, isolate fitting range and convert to log-log scale for fitting
relevant_freq_inds = find(frequencies >= fitting_range(1) & frequencies <= fitting_range(2));
freqs = frequencies(relevant_freq_inds);
log_freqs = log10(freqs);
log_spec = log10(spectrum(relevant_freq_inds));

%next, find the point to fit to the bottom edge of the fitting range
endpoint_inds = find(freqs <= fitting_range(1)+0.6);
endpoint = [log_freqs(1),mean(log_spec(endpoint_inds))];

%test if there is a notch filter
line_point = find(freqs >= 60);
if isempty(line_point)
    notch_test = 0;
else
    notch_test = log_spec(line_point(1)) < log_spec(end);
end

%then, find the point to fit to the rest of the fitting range, which is
%affected by whether there is a notch filter. If there is no notch filter,
%scan the whole range in increments of 0.05 Log Frequency and find the
%minimum value within each section.
if notch_test == 0
initialize = 0.9;    
slope_upper_limit = log_freqs(end)-0.05;
num = round(10*(slope_upper_limit-initialize)*2);
firstpoint_temp = zeros(num,2);
slope_temp = zeros(num,1);
for iter = 1:num
    firstpoint_range = find(log_freqs >= initialize+(0.05*(iter-1)) & log_freqs <= initialize+0.05+(0.05*(iter-1)));
    [firstpoint_power_temp,firstpoint_ind] = min(log_spec(firstpoint_range));
    firstpoint_power = mean([log_spec(firstpoint_range(firstpoint_ind)-1),firstpoint_power_temp,log_spec(firstpoint_range(firstpoint_ind)+1)]);
    firstpoint = [log_freqs(firstpoint_range(firstpoint_ind)),firstpoint_power];
    slope_temp(iter) = (firstpoint(2)-endpoint(2))./(firstpoint(1)-endpoint(1));
    firstpoint_temp(iter,:) = firstpoint;
end 

% If there is a notch filter, do not inclue 55-65 Hz in the scan.
elseif notch_test == 1
initialize = 0.9; %search above 8 Hz     
slope_int_limit = 1.74-0.05; % about 55 Hz -0.05 Hz buffer 
slope_int_start = 1.813; %about 65 Hz
slope_upper_limit = log_freqs(end)-0.05;
num = round(10*(slope_int_limit-initialize)*2);
firstpoint_temp_1 = zeros(num,2);
slope_temp_1 = zeros(num,1);
for iter = 1:num
    firstpoint_range = find(log_freqs >= initialize+(0.05*(iter-1)) & log_freqs <=initialize+0.05+(0.05*(iter-1)));
    [firstpoint_power_temp,firstpoint_ind] = min(log_spec(firstpoint_range));
    firstpoint_power = mean([log_spec(firstpoint_range(firstpoint_ind)-1),firstpoint_power_temp,log_spec(firstpoint_range(firstpoint_ind)+1)]);
    firstpoint = [log_freqs(firstpoint_range(firstpoint_ind)),firstpoint_power];
    slope_temp_1(iter) = (firstpoint(2)-endpoint(2))./(firstpoint(1)-endpoint(1));
    firstpoint_temp_1(iter,:) = firstpoint;
end  
num = round(10*(slope_upper_limit-slope_int_start)*2);
firstpoint_temp_2 = zeros(num,2);
slope_temp_2 = zeros(num,1);
for iter = 1:num
    firstpoint_range = find(log_freqs >= slope_int_start+(0.05*(iter-1)) & log_freqs <= slope_int_start+0.05+(0.05*(iter-1)));
    [firstpoint_power_temp,firstpoint_ind] = min(log_spec(firstpoint_range));
    firstpoint_power = mean([log_spec(firstpoint_range(firstpoint_ind)-1),firstpoint_power_temp,log_spec(firstpoint_range(firstpoint_ind)+1)]);
    firstpoint = [log_freqs(firstpoint_range(firstpoint_ind)),firstpoint_power];
    slope_temp_2(iter) = (firstpoint(2)-endpoint(2))./(firstpoint(1)-endpoint(1));
    firstpoint_temp_2(iter,:) = firstpoint;
end  
slope_temp = [slope_temp_1; slope_temp_2];
firstpoint_temp = [firstpoint_temp_1; firstpoint_temp_2];
end

%of these possibilities, the best point is the one that maximizes the slope
[~,ind] = max(abs(slope_temp));
slope = slope_temp(ind);
firstpoint = firstpoint_temp(ind,:);

%use this point to find the aperiodic parameters
offset = firstpoint(2)-(slope*firstpoint(1));
chai = -1*slope;

%find aperiodic and periodic spectra
figure(623)
plot(log_freqs,log_spec)
hold on
line_output = (slope*log_freqs)+offset;
plot(log_freqs,line_output);
aperiodic = line_output;
plot(log_freqs,aperiodic);
hold off
periodic = log_spec-aperiodic;

%next, find the periodic parameters
figure(625)
% for mouse EEG data (requires smoothing)
if num_peaks < 2
 first_peak_inds = find(freqs > fitting_range(1) & freqs < 10);
 smoothed_theta = smoothdata(periodic(first_peak_inds),'lowess');
 plot(freqs(first_peak_inds),periodic(first_peak_inds))
 hold on
 plot(freqs(first_peak_inds),smoothed_theta)
[~,ind_temp] = max(smoothed_theta);
% isolate around the max identified from the filtered data
snip = length(first_peak_inds)-(ind_temp+5);
if snip >= 0
test_range = ind_temp-5:ind_temp+5;
else
test_range = ind_temp-5:ind_temp+snip+5;  
end
%find the actual peak within this range
[~,ind] = max(periodic(first_peak_inds(test_range)));
first_peak_freqs_temp = freqs(first_peak_inds(test_range(ind)));
first_peak_params = [first_peak_freqs_temp,periodic(first_peak_inds(test_range(ind)))];
scatter(first_peak_freqs_temp,periodic(first_peak_inds(test_range(ind))));

if first_peak_params(1) > 8
prompt = "It appears that the identified peak is greater than 8 Hz, which makes in an outlier. This could be because the fit is poor. Do you like this fit? Y = 1, No = 0. If not, please click on the graph where the peak is.";    
x = input(prompt);
      if x == 0 
          figure(625)
          [theta_peak_freq_temp,powt] = ginput(1);
          scatter(theta_peak_freq_temp,powt)
          first_peak_params = [theta_peak_freq_temp,powt];          
      end
end
hold off

second_peak_params = [0 0];

else
%find peak 1, which is the maximum value below 15 Hz   
figure(625)
 first_peak_inds = find(freqs > fitting_range(1) & freqs < 15);
 plot(freqs(first_peak_inds),periodic(first_peak_inds))
 hold on
 [pow,ind] = max(periodic(first_peak_inds(3:end-2)));
 first_peak_params = [freqs(first_peak_inds(ind)+2), pow];
 scatter(freqs(first_peak_inds(ind)+2),periodic(first_peak_inds(ind)+2))

 %it is unusal to find a peak about 12 Hz. The peak might be lower but the
 %spectrum could have an usual shape, causing a misidentification.
if first_peak_params(1) >= 12 
     prompt = "It appears that the identified peak is greater than 12 Hz, which makes in an outlier. This could be because the fit is poor. Do you like this fit? Y = 1, No = 0. If not, please click on the graph where the peak is.";
     x = input(prompt);
     if x == 0
      [theta_peak_freq_temp,powt] = ginput(1);
      scatter(theta_peak_freq_temp,powt)
      first_peak_params = [theta_peak_freq_temp,powt];
     end
end
hold off

%find periodic peak 2, the maximum power between 15 and 30 Hz
figure(626)
second_peak_inds = find(freqs > 15 & freqs < 30);
 plot(freqs(second_peak_inds),periodic(second_peak_inds))
 hold on
 [pow,ind] = max(periodic(second_peak_inds));
 second_peak_params = [freqs(second_peak_inds(ind)), pow];
 scatter(freqs(second_peak_inds(ind)),periodic(second_peak_inds(ind)))
 hold off    
end

