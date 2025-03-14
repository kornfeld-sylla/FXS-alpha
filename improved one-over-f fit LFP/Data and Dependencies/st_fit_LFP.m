function [offset,chai,knee_freq,log_spec,freqs,aperiodic,periodic,peak_params] = st_fit_LFP(frequencies, spectrum, fitting_range, peak_analysis)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

%first, isolate frequency range and convert to log-log scale for fitting
relevant_freq_inds = find(frequencies >= fitting_range(1) & frequencies <= fitting_range(2));
freqs = frequencies(relevant_freq_inds);
log_freqs = log10(freqs);
log_spec = log10(spectrum(relevant_freq_inds));

%next, find the point to fit to the upper edge of the fitting range
non_notched_freqs = find(freqs < 45 | freqs > 75);
[endpoint_pow_temp,endpoint_ind] = min(log_spec(non_notched_freqs));
if non_notched_freqs(endpoint_ind) < length(log_spec)
    endpoint_pow = mean([log_spec(non_notched_freqs(endpoint_ind)-1),endpoint_pow_temp,log_spec(non_notched_freqs(endpoint_ind)+1)]);
else
    endpoint_pow = mean([log_spec(non_notched_freqs(endpoint_ind)-2),endpoint_pow_temp,log_spec(non_notched_freqs(endpoint_ind)-1)]);
end
endpoint = [log_freqs(non_notched_freqs(endpoint_ind)), endpoint_pow];

%fit second point to problematic range, 9-25 Hz
%step through this range in units of 0.05 Log Frequency and identify the
%smallest point in the pin. Calculate a tentative slope from that point.
fit_range = [0.95 1.4];
numby = round(10*(fit_range(2)-fit_range(1))*2);
firstpoint_temp = zeros(numby,2);
slope_temp = zeros(numby,1);
for nash = 1:numby
    firstpoint_range = find(log_freqs >= fit_range(1)+(0.05*(nash-1)) & log_freqs <= fit_range(1)+0.05+(0.05*(nash-1)));
    [firstpoint_power_temp,firstpoint_ind] = min(log_spec(firstpoint_range));
    firstpoint_power = mean([log_spec(firstpoint_range(firstpoint_ind)-1),firstpoint_power_temp,log_spec(firstpoint_range(firstpoint_ind)+1)]);
    firstpoint = [log_freqs(firstpoint_range(firstpoint_ind)),firstpoint_power];
    slope_temp(nash) = (firstpoint(2)-endpoint(2))./(firstpoint(1)-endpoint(1));
    firstpoint_temp(nash,:) = firstpoint;
end  

%find the coordinate pair which minimizes the slope, as a less steep slope
%will not cut through the middle portions of the spectrum. Use this to find
%the slope and offset of the line
[~,ind] = min(abs(slope_temp));
slope = slope_temp(ind);
firstpoint = firstpoint_temp(ind,:);
offset = firstpoint(2)-(slope*firstpoint(1));
chai = -1*slope;

figure(623)
plot(log_freqs,log_spec)
hold on
line_output = (slope*log_freqs)+offset;
plot(log_freqs,line_output);
scatter(firstpoint(1),firstpoint(2))

% find the knee parameter based on local minima below 2.24 Hz
knee_inds = find(log_freqs < 0.35);
knee_data = -1*log_spec(knee_inds);
peak_locations = findpeaks(knee_data);
gramp = length(peak_locations.loc);
mean_knee_est = [];
for grump = 1:gramp
[~,peak_loc] = max(knee_data(peak_locations.loc));
flatpoint = [freqs(knee_inds(peak_locations.loc(peak_loc))),log_spec(knee_inds(peak_locations.loc(peak_loc)))];
if flatpoint(2) < line_output(knee_inds(peak_locations.loc(peak_loc))) %find point that is below y = mx+b line for knee
    scatter(log10(flatpoint(1)),flatpoint(2));
    knee_est_temp = (10^(offset-flatpoint(2)))-(flatpoint(1).^chai);
    mean_knee_est = [mean_knee_est, knee_est_temp];
    peak_locations.loc(peak_loc) = [];
else
    peak_locations.loc(peak_loc) = [];
end
end

%if there were no points that fell below the aperiodic line, then there is
%no knee. Otherwise, the knee value is the mean of the points. The knee
%frequency can be calulated from the knee value and the slope
if isempty(mean_knee_est)
    aperiodic = line_output;
    knee_freq = 0;
else
    knee_est = mean(mean_knee_est);
    knee_freq = knee_est^(1/chai);
    aperiodic = offset - log10(knee_est + (freqs.^chai));
end

plot(log_freqs,aperiodic);
hold off
periodic = log_spec-aperiodic;

% Next, analyze peak params if asked to do so
if peak_analysis > 0
% find frequencies below 10 Hz for peaks    
figure(625)
low_freq_inds = find(freqs > fitting_range(1) & freqs < 10);
smoothed_data = smoothdata(periodic(low_freq_inds),'lowess');
plot(freqs(low_freq_inds),periodic(low_freq_inds))
hold on
plot(freqs(low_freq_inds),smoothed_data)

%identify local minima
testing = findpeaks(-1*smoothed_data);
mid_min = find(freqs(low_freq_inds(testing.loc)) >= 3 & freqs(low_freq_inds(testing.loc)) <= 6);

%find the maximum frequency below 6 Hz
    a_inds = find(freqs < 6);
    smoothed_data_a = smoothdata(periodic(a_inds),'lowess');
    [~,ind] = max(smoothed_data_a);
    snip = a_inds(end)-(ind+5);
    if snip < 0
    test_range = ind-5:ind+snip+5;
    elseif ind <= 5
    test_range = 1:ind+5;
    else
    test_range = ind-5:ind+5;
    end 
      [~,peak_1] = max(periodic(test_range));
      peak_1_freqs_temp = freqs(test_range(peak_1));
      pow_1 = periodic(test_range(peak_1));
      scatter(peak_1_freqs_temp,pow_1)
      prompt = "Do you like this fit? Y = 1, N = 0. If not, you'll either be given a second option or be prompted to manually identify by clicking on the graph.";
      x = input(prompt);
      if x == 0 
          if isempty(mid_min)
          figure(625)
          [peak_1_freqs_temp,pow_1] = ginput(1);
          scatter(peak_1_freqs_temp,pow_1)
          else
               [~,peak_1_est] = max(smoothed_data(1:testing.loc(mid_min(1))));
                 if peak_1_est <= 5
                 test_range = 1:peak_1_est+5;
                 else
                test_range = peak_1_est-5:peak_1_est+5;
                 end 
                [~,peak_1] = max(periodic(low_freq_inds(test_range)));
                peak_1_freqs_temp = freqs(low_freq_inds(test_range(peak_1))); 
                pow_1 = periodic(low_freq_inds(test_range(peak_1)));
                scatter(peak_1_freqs_temp,pow_1);
                prompt = "Is this better? Y = 1, N = 0. If not, please click on the graph where the peak should be. ";
                x = input(prompt);
                    if x == 0 
                        figure(625)
                        [peak_1_freqs_temp,pow_1] = ginput(1);
                        scatter(peak_1_freqs_temp,pow_1)
                    end
          end
      end
    b_inds = find(freqs > 6 & freqs < 10);
    smoothed_data_b = smoothdata(periodic(b_inds),'lowess');
    [~,ind] = max(smoothed_data_b);
    snip = low_freq_inds(end)-(b_inds(ind)+5);
     if snip < 0
         test_range = b_inds(ind)-5:b_inds(ind)+snip+5;
     elseif ind == 1
         test_range = b_inds(ind):b_inds(ind)+5;
     else
     test_range = b_inds(ind)-5:b_inds(ind)+5;
     end
      [~,peak_2] = max(periodic(test_range));
      peak_2_temp = freqs(test_range(peak_2));
      pow_2 = periodic(test_range(peak_2));
      scatter(peak_2_temp,pow_2)
      prompt = "Do you like this fit? Y = 1, N = 0.If not, you'll either be given a second option or be prompted to manually identify by clicking on the graph.";
      x = input(prompt);
      if x == 0
          if isempty(mid_min)
             figure(625)
            [freq_peak_2,pow] = ginput(1);
            scatter(freq_peak_2,pow)
            peak_params = [peak_1_freqs_temp, freq_peak_2, pow_1, pow]; 
          else
              [~,lowest_min] = min(smoothed_data(testing.loc(mid_min)));
              [~,peak_2_est] = max(smoothed_data(testing.loc(mid_min(lowest_min)):end));
             snip = length(low_freq_inds)-(testing.loc(mid_min(lowest_min))+peak_2_est-1+5);
             if snip >= 0
                test_range = testing.loc(mid_min(lowest_min))+peak_2_est-1-5:testing.loc(mid_min(lowest_min))+peak_2_est-1+5;
             else
                test_range = testing.loc(mid_min(lowest_min))+peak_2_est-1-5:testing.loc(mid_min(lowest_min))+peak_2_est-1+snip+5;
             end
                 [~,peak_2] = max(periodic(low_freq_inds(test_range)));
                 peak_2_temp = freqs(low_freq_inds(test_range(peak_2))); 
                 pow_2 = periodic(low_freq_inds(test_range(peak_2)));
                scatter(peak_2_temp,pow_2);
                prompt = "Do you like this fit better? Y = 1, N = 0. If not, please click on the graph where the peak should be.";
                x = input(prompt);
                if x == 0
                    figure(625)
                    [freq_peak_2,pow] = ginput(1);
                    scatter(freq_peak_2,pow)
                    peak_params = [peak_1_freqs_temp, freq_peak_2, pow_1, pow]; 
                else
                    peak_params = [peak_1_freqs_temp, peak_2_temp, pow_1, pow_2];
                end
          end            
      else
      peak_params = [peak_1_freqs_temp, peak_2_temp, pow_1, pow_2];
      end
else
peak_params = [0 0 0 0];    
end
hold off


end

