function [mean_burst_amp,num_bursts,mean_burst_length] = BURST_ANALYSIS(bandpass_analysis,threshold_val)
% Identify burst count and duration, as well as mean burst amplitude 
%  Input is the hilbert transformed amplitude values, and the amplitude
%  threshold value for determining a burst
%  Output is the average amplitude of the bursts,the number of bursts, and the average burst length 

%plot hilbert transform amplitude values with threshold
figure(1000);
plot(bandpass_analysis,'b');
hold on
[cum_dist, power_vals] = ecdf(bandpass_analysis); 
toot = find(cum_dist > threshold_val);
bandpass_threshold = power_vals(toot(1));
yline(bandpass_threshold);

%determine burst count and burst duration
%find values that cross threshold (bursts)
bursts = find(bandpass_analysis > power_vals(toot(1)));
burst_tester = [bursts(1) bursts(1:end-1)];
jsn = bursts-burst_tester;
%find the beginning of each burst
lmp = find(jsn > 1); 
num_bursts = length(lmp) + 1;
if isempty(lmp)
    mean_burst_length = length(bursts);
else
%examine the first burst
burst_length_temp = length(bursts(1:lmp(1))-1);
%exclude burst if it is to short (less than 3 ms)
if burst_length_temp < 3
    bursts(1:lmp(1)) = zeros(1,burst_length_temp);
    num_bursts = num_bursts - 1;
    burst_length = [];
else
    burst_length = burst_length_temp;
end
for ken = 1:length(lmp)-1
    %examine rest of bursts
    rst = length(bursts(lmp(ken):(lmp(ken+1)-1)));
    %exclude bursts that are too short (less than 3ms)
    if rst < 3
        bursts(lmp(ken):lmp(ken+1)-1) = zeros(1,rst);
        num_bursts = num_bursts - 1;
    else
    burst_length = [burst_length, rst];  
    end
end
%exclude last burst if it is too short (less than 3ms)
if length(bursts(lmp(end):end)) < 3
    bursts(lmp(end):end) = zeros(1,length(bursts(lmp(end):end)));
    num_bursts = num_bursts - 1;
else
burst_length = [burst_length, length(bursts(lmp(end):end))];
end
mean_burst_length = mean(burst_length);
end

%remove burst amplitude values that were set to zero because the bursts
%were excluded (less than 3ms)
chop = find(bursts < 1);
bursts(chop) = [];

%plot burst values surperimposed on hilbert transform values
burst_vals = bandpass_analysis(bursts);
mean_burst_amp = mean(burst_vals);
scatter(bursts,burst_vals,1,'m*'); 

hold off

end

