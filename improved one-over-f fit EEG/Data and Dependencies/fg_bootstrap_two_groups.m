function [x_data,xxline] = fg_bootstrap_two_groups(fig_num,B_L,group_1_data,group_2_data,group_1_counter,group_2_counter,log_scale,frequencies,CIU,CIL,freq_lims)
%Use bootstrapping for statistic of power spectra
%output gives the coordinates for points of significant difference
%input is figure number, the number of bootstrap repetitions (i.e., 1000),
%the two power spectra samples, any counters for number of
%trials per animal (if applicable), whether the output should be plotted in
%log scale (1 - Yes, 0 - No), the frequencies for the spectra, the
%confidence interval low and upper limits,
%and the frequency range over which to plot the resulting bootstrapped spectra

span = ndims(group_1_data);
%preallocate
bootstrapped_diff = zeros(B_L,size(group_1_data,span));
%determine if samples are the same size, which decides if the trials are
%resampled from the same random distribution (i.e., matched) or different
%distributions (i.e., independent)
if size(group_1_data,1) == size(group_2_data,1)
for boot = 1:B_L
    %select random subset of subjects
    boop = randi(size(group_1_data,1),size(group_1_data,1),1);
    %preallocate
    if span < 3
    group1_temp = zeros(size(group_1_data,1),size(group_1_data,span));
    group2_temp = zeros(size(group_2_data,1),size(group_2_data,span));
    elseif span == 3
    group1_temp = zeros(size(group_1_data,1)*size(group_1_data,2),size(group_1_data,span));
    group2_temp = zeros(size(group_2_data,1)*size(group_2_data,2),size(group_2_data,span));
    end
    for loop = 1:length(boop)
        %if 2-d data matrix, iterate through the random subset of subjects in boop
        if span < 3
        group1_temp(loop,:) = group_1_data(boop(loop),:);
        group2_temp(loop,:) = group_2_data(boop(loop),:);
        %if 3-d data matix, do another round of heirarchical bootstrapping
        %of trials per subject
        elseif span == 3
        %select a random subset of trials for each randomly selected subject in the boop list
        beep = randi(size(group_1_data,2),size(group_1_data,2),1);
            %loop through these trials and reshape into a 2-d data matrix,
            %as all beep trials from all boop subjects are averaged together
            %for each bootstrapped sample
            for leep = 1:length(beep)
            group1_temp(leep+(size(group_1_data,2)*(loop-1)),:) = reshape(group_1_data(boop(loop),beep(leep),:),[1,size(group_1_data,span)]);
            group2_temp(leep+(size(group_2_data,2)*(loop-1)),:) = reshape(group_2_data(boop(loop),beep(leep),:),[1,size(group_2_data,span)]);
            end
        end
    end
    %find the difference between group means for each boostrapped sample
    if log_scale == 1
    bootstrapped_diff(boot,:) = mean(group1_temp,1)./mean(group2_temp,1);
    else
    bootstrapped_diff(boot,:) = mean(group1_temp,1) - mean(group2_temp,1);
    end
end

%if they are not the same size, determine if there is a counter to
%breakdown by subject
else
for boot = 1:B_L
    %if there is not a counter, it means each row corresponds to a
    %different subjects and there are a different number of subjects per
    %group, meaning the datasets must be independent
    if isempty(group_1_counter) && isempty(group_2_counter)
        %generate two different random samples of animals for the
        %independent datasets
        boop_1 = randi(size(group_1_data,1),size(group_1_data,1),1);
        boop_2 = randi(size(group_2_data,1),size(group_2_data,1),1);
    %preallocate
    if span < 3
    group1_temp = zeros(size(group_1_data,1),size(group_1_data,span));
    group2_temp = zeros(size(group_2_data,1),size(group_2_data,span));
    elseif span == 3
    group1_temp = zeros(size(group_1_data,1)*size(group_1_data,2),size(group_1_data,span));
    group2_temp = zeros(size(group_2_data,1)*size(group_2_data,2),size(group_2_data,span));
    end
    %start with the first random subject list
    for loop_1 = 1:length(boop_1)
        if span < 3
            %as before, if it is a 2-d matrix, select out the random subset
            %of subjects
        group1_temp(loop_1,:) = group_1_data(boop_1(loop_1),:);
        elseif span == 3
            %as before, if it is a 3-d matrix, perform heirarchical
            %bootstrapping and select out a random subset of trials for
            %each randomly selected subject
        beep_1 = randi(size(group_1_data,2),size(group_1_data,2),1);
            for leep_1 = 1:length(beep_1)
            group1_temp(leep_1+(size(group_1_data,2)*(loop_1-1)),:) = reshape(group_1_data(boop_1(loop_1),beep_1(leep_1),:),[1,size(group_1_data,span)]);
            end
        end
    end
    %repeat for the second random subject list
    for loop_2 = 1:length(boop_2)
        if span < 3
        group2_temp(loop_2,:) = group_2_data(boop_2(loop_2),:);
        elseif span == 3
        beep_2 = randi(size(group_2_data,2),size(group_2_data,2),1);
            for leep_2 = 1:length(beep_2)
            group2_temp(leep_2+(size(group_2_data,2)*(loop_2-1)),:) = reshape(group_2_data(boop_2(loop_2),beep_2(leep_2),:),[1,size(group_2_data,span)]);
            end
        end
    end

    %if the counters are not empty, this means that the group data
    %matricies are organized by trials, with the number of trials per
    %subject corresponding to the value in the counter array per subject
    else
    %generate random sample of animals
    boop_1 = randi(length(group_1_counter),length(group_1_counter),1);
    boop_2 = randi(length(group_2_counter),length(group_2_counter),1);
    %preallocate
    if size(group_1_data,1) > length(group_1_counter) && size(group_2_data,1) > length(group_2_counter)
    group1_temp = [];
    group2_temp = [];
    else
        group1_temp = zeros(length(group_1_counter),size(group_1_data,span));
        group2_temp = zeros(length(group_2_counter),size(group_2_data,span));
    end
    for loop_1 = 1:length(boop_1)
        if size(group_1_data,1) > length(group_1_counter) 
        %find the randomly selected subjects from the counter list 
        num_trials_animal = group_1_counter(boop_1(loop_1));
        previous_trials = sum(group_1_counter(1:boop_1(loop_1)-1));
        data_list = group_1_data(previous_trials+1:previous_trials+num_trials_animal,:);
        %because the data is heirarchical by nature, randomly select a
        %subset of trials per randomly selected subject
        beep_1 = randi(num_trials_animal,num_trials_animal,1);
        %concatenate the selected trials for each subject
        group1_temp = [group1_temp; data_list(beep_1,:)];
        else
        group1_temp(loop_1,:) = group_1_data(boop_1(loop_1),:);
        end
    end
    %repeat for the second, independent set of subjects
    for loop_2 = 1:length(boop_2)
        if size(group_2_data,1) > length(group_2_counter) 
        num_trials_animal = group_2_counter(boop_2(loop_2));
        previous_trials = sum(group_2_counter(1:boop_2(loop_2)-1));
        data_list = group_2_data(previous_trials+1:previous_trials+num_trials_animal,:);
        beep_2 = randi(num_trials_animal,num_trials_animal,1);
        group2_temp = [group2_temp; data_list(beep_2,:)];
        else
        group2_temp(loop_2,:) = group_2_data(boop_2(loop_2),:);
        end
    end
    end
    %find the difference between group means for each boostrapped sample
    if log_scale == 1
    bootstrapped_diff(boot,:) = mean(group1_temp,1)./mean(group2_temp,1);
    else
    bootstrapped_diff(boot,:) = mean(group1_temp,1) - mean(group2_temp,1);
    end
end    
end

%plot the boostrapped difference, with the median difference represented as
%a line and +/- confidence intervals shaded on each side
sorted_boot = sort(bootstrapped_diff,1,'ascend');
figure(fig_num)
color = [0.6 0.6 0.6];
x_axis = frequencies;
%plot median
if log_scale == 1
x = 10*log10(sorted_boot(CIU,:));
y = 10*log10(sorted_boot(CIL,:));
plot(frequencies,10*log10(sorted_boot(500,:)), 'k', 'LineWidth', 0.8);
else
x = (sorted_boot(CIU,:));
y = (sorted_boot(CIL,:));
if size(sorted_boot,2) == 18
plot(frequencies,(sorted_boot(500,:)), 'k.-', 'LineWidth', 0.8);
else
plot(frequencies,(sorted_boot(500,:)), 'k', 'LineWidth', 0.8);
end
end
%plot confidence intervals
c = color;
h = patch([x_axis, fliplr(x_axis)], [x, fliplr(y)], c);
hold on
plot(frequencies,x,'k--');
plot(frequencies,y,'k--');
set(h, 'facealpha', 0.15);
set(h, 'edgecolor', c);
set(h, 'linestyle', '--');
ylabel('Bootstrapped Difference')
xlabel('Frequency')
yline(0);
relevant_indicies = find(frequencies >= freq_lims(1) & frequencies <= freq_lims(2));
if log_scale == 1
pin = find(10*log10(sorted_boot(CIL,relevant_indicies)) > 0);
nip = find(10*log10(sorted_boot(CIU,relevant_indicies)) < 0);
else
pin = find((sorted_boot(CIL,relevant_indicies)) > 0);
nip = find((sorted_boot(CIU,relevant_indicies)) < 0);   
end
new_ax = frequencies(relevant_indicies);
x_data = new_ax([pin,nip]);
xxline = ones(1,length([pin,nip]));

end

