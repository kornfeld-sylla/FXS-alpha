function plot_spectrum_two_groups(fig_num,group_1_2d_data,group_2_2d_data,frequencies,logical_log)
%Plot power spectra with mean and SEM. Group one is plotted in red, group
%two in black.

mean_group_1 = mean(group_1_2d_data,1);
SE_group_1 = zeros(1,length(mean_group_1));
for jj = 1:length(mean_group_1)
SE_group_1(jj) = std((group_1_2d_data(:,jj)))./sqrt(length(group_1_2d_data(:,jj)));
end
mean_group_2 = mean(group_2_2d_data,1);
SE_group_2 = zeros(1,length(mean_group_2));
for jj = 1:length(mean_group_2)
SE_group_2(jj) = std((group_2_2d_data(:,jj)))./sqrt(length(group_2_2d_data(:,jj)));
end

figure(fig_num)
color_1 = [1 0.8 0];
color_2 = [0.6 0.6 0.6];
x_axis = frequencies;
if logical_log == 1
x = 10*log10(mean_group_1+SE_group_1);
y = 10*log10(mean_group_1-SE_group_1);
else
x = mean_group_1+SE_group_1;
y = mean_group_1-SE_group_1;
end
h = patch([x_axis, fliplr(x_axis)], [x, fliplr(y)], color_1);
hold on
plot(frequencies,x,'r--');
plot(frequencies,y,'r--');
set(h, 'facealpha', 0.15);
set(h, 'edgecolor', color_1);
set(h, 'linestyle', '--');
if logical_log == 1
plot(frequencies,10*log10(mean_group_1), 'r', 'LineWidth', 0.8);
xx = 10*log10(mean_group_2+SE_group_2);
yy = 10*log10(mean_group_2-SE_group_2);
else
plot(frequencies,mean_group_1, 'r', 'LineWidth', 0.8);
xx = mean_group_2+SE_group_2;
yy = mean_group_2-SE_group_2;
end
h = patch([x_axis, fliplr(x_axis)], [xx, fliplr(yy)], color_2);
plot(frequencies,xx,'k--');
plot(frequencies,yy,'k--');
set(h, 'facealpha', 0.15);
set(h, 'edgecolor', color_2);
set(h, 'linestyle', '--');
if logical_log == 1
plot(frequencies,10*log10(mean_group_2), 'k', 'LineWidth', 0.8);
else
plot(frequencies,mean_group_2, 'k', 'LineWidth', 0.8);
hold off
xlabel('Frequency (Hz)')
ylabel('Spectrum Power (dB)')
end

