function [x_data,xxline] = bootstrap_two_groups(fig_num,subplot_info, B_L,diff,group_1_data,group_2_data,group_1_counter,group_2_counter,log_scale,frequencies,CIU,CIL,freq_lims)
%UNTITLED Summary of this function goes here
%First, sample across same number of frequency bins and animals for power
%spectrum
%for bootstrapping mod index, include second version of each graph that
%just highlights portions that are significantly more than 0
spank = ndims(group_1_data);
if diff ~= 1
    bootstrapped_group_1 = zeros(B_L,1);
    bootstrapped_group_2 = zeros(B_L,1);
    diff_test = zeros(B_L,1);
if size(group_1_data,1) == size(group_2_data,1)    
for boot = 1:B_L
    boop = randi(size(group_1_data,1),size(group_1_data,1),1);    
    group1_temp = zeros(size(group_1_data,1)*size(group_1_data,spank),1);
    group2_temp = zeros(size(group_2_data,1)*size(group_2_data,spank),1);
    for loop = 1:length(boop)
        beep = randi(size(group_1_data,2),size(group_1_data,2),1);
        for leep = 1:length(beep)
        group1_temp(leep+(size(group_1_data,2)*(loop-1))) = group_1_data(boop(loop),beep(leep));
        group2_temp(leep+(size(group_2_data,2)*(loop-1))) = group_2_data(boop(loop),beep(leep));
        end
    end
bootstrapped_group_1(boot,:) = mean(group1_temp,1);
bootstrapped_group_2(boot,:) = mean(group2_temp,1);
diff_test(boot,:) = mean(group1_temp,1) - mean(group2_temp,1);
end
else
 for boot = 1:B_L
     if isempty(group_1_counter) && isempty(group_2_counter)
     boop_1 = randi(size(group_1_data,1),size(group_1_data,1),1); 
     boop_2 = randi(size(group_2_data,1),size(group_2_data,1),1); 
     group1_temp = zeros(size(group_1_data,1)*size(group_1_data,spank),1);
     group2_temp = zeros(size(group_2_data,1)*size(group_2_data,spank),1);
     for loop_1 = 1:length(boop_1)
        beep_1 = randi(size(group_1_data,2),size(group_1_data,2),1);
        for leep_1 = 1:length(beep_1)
        group1_temp(leep_1+(size(group_1_data,2)*(loop_1-1))) = group_1_data(boop_1(loop_1),beep_1(leep_1));
        end
     end
     for loop_2 = 1:length(boop_2)
         beep_2 = randi(size(group_2_data,2),size(group_2_data,2),1);
         for leep_2 = 1:length(beep_2)
        group2_temp(leep_2+(size(group_2_data,2)*(loop_2-1))) = group_2_data(boop_2(loop_2),beep_2(leep_2));
         end
     end
     else
     boop_1 = randi(length(group_1_counter),length(group_1_counter),1);
     boop_2 = randi(length(group_2_counter),length(group_2_counter),1);
     group1_temp = [];
     group2_temp = [];
     for loop_1 = 1:length(boop_1)
     num_trials_animal = group_1_counter(boop_1(loop_1));
        previous_trials = sum(group_1_counter(1:boop_1(loop_1)-1));
        data_list = group_1_data(previous_trials+1:previous_trials+num_trials_animal,:);
        beep_1 = randi(num_trials_animal,num_trials_animal,1);
            for leep_1 = 1:length(beep_1)  
                group1_temp = [group1_temp; data_list(beep_1(leep_1),:)];
            end
     end
     for loop_2 = 1:length(boop_2)
        num_trials_animal = group_2_counter(boop_2(loop_2));
        previous_trials = sum(group_2_counter(1:boop_2(loop_2)-1));
        data_list = group_2_data(previous_trials+1:previous_trials+num_trials_animal,:);
        beep_2 = randi(num_trials_animal,num_trials_animal,1);
            for leep_2 = 1:length(beep_2)  
                group2_temp = [group2_temp; data_list(beep_2(leep_2),:)];
            end
     end
     end
bootstrapped_group_1(boot,:) = mean(group1_temp,1);
bootstrapped_group_2(boot,:) = mean(group2_temp,1); 
diff_test(boot,:) = mean(group1_temp,1) - mean(group2_temp,1);
end   
end

sorted_boot_1 = sort(bootstrapped_group_1,1,'ascend');
sorted_boot_2 = sort(bootstrapped_group_2,1,'ascend');
sorted_diff_test = sort(diff_test,1,'ascend');
if isempty(subplot_info) == false
figure(fig_num)
subplot(1,subplot_info(1),subplot_info(2))
h = boxplot(sorted_boot_1,'Whisker',0);
 set(h(7,:), 'Visible', 'off')
 set(h(3,:), 'Visible', 'off')
 set(h(4,:), 'Visible', 'off')
 set(h(5,1), 'YData',[sorted_boot_1(CIL) sorted_boot_1(CIU) sorted_boot_1(CIU) sorted_boot_1(CIL) sorted_boot_1(CIL)]); %changes the box itself
 set(h(5,1), 'Color', [1 0 0]); set(h(6,1), 'Color', [1 0 0]);
 hold on
h1 = boxplot(sorted_boot_2,'Whisker',0); 
set(h1(7,:), 'Visible', 'off')
 set(h1(3,:), 'Visible', 'off')
 set(h1(4,:), 'Visible', 'off')
  set(h1(5,1), 'YData',[sorted_boot_2(CIL) sorted_boot_2(CIU) sorted_boot_2(CIU) sorted_boot_2(CIL) sorted_boot_2(CIL)]); %changes the box itself
 set(h1(5,1), 'Color', [0 0 0]); set(h1(6,1), 'Color', [0 0 0]);
 
figure(fig_num+100);
subplot(1,subplot_info(1),subplot_info(2)) 
h = boxplot(sorted_diff_test,'Whisker',0);
 set(h(7,:), 'Visible', 'off')
 set(h(3,:), 'Visible', 'off')
 set(h(4,:), 'Visible', 'off')
 set(h(5,1), 'YData',[sorted_diff_test(CIL) sorted_diff_test(CIU) sorted_diff_test(CIU) sorted_diff_test(CIL) sorted_diff_test(CIL)]); %changes the box itself
 set(h(5,1), 'Color', [0.7 0.7 0.7]); set(h(6,1), 'Color', [0.7 0.7 0.7]);
else
figure(fig_num)    
h = boxplot([sorted_boot_1, sorted_boot_2],{'FX' 'WT'},'Whisker',0);
 set(h(7,:), 'Visible', 'off')
 set(h(3,:), 'Visible', 'off')
 set(h(4,:), 'Visible', 'off')
 set(h(5,1), 'YData',[sorted_boot_1(CIL) sorted_boot_1(CIU) sorted_boot_1(CIU) sorted_boot_1(CIL) sorted_boot_1(CIL)]); %changes the box itself
 set(h(5,1), 'Color', [1 0 0]); set(h(6,1), 'Color', [1 0 0]);
 set(h(5,2), 'YData',[sorted_boot_2(CIL) sorted_boot_2(CIU) sorted_boot_2(CIU) sorted_boot_2(CIL) sorted_boot_2(CIL)]); %changes the box itself
 set(h(5,2), 'Color', [0 0 0]); set(h(6,2), 'Color', [0 0 0]);
 
 figure(fig_num+100)
 h = boxplot(sorted_diff_test,'Whisker',0);
 set(h(7,:), 'Visible', 'off')
 set(h(3,:), 'Visible', 'off')
 set(h(4,:), 'Visible', 'off')
 set(h(5,1), 'YData',[sorted_diff_test(CIL) sorted_diff_test(CIU) sorted_diff_test(CIU) sorted_diff_test(CIL) sorted_diff_test(CIL)]); %changes the box itself
 set(h(5,1), 'Color', [0.7 0.7 0.7]); set(h(6,1), 'Color', [0.7 0.7 0.7]);
end

else 
bootstrapped_diff = zeros(B_L,size(group_1_data,spank));
if size(group_1_data,1) == size(group_2_data,1)
for boot = 1:B_L
    boop = randi(size(group_1_data,1),size(group_1_data,1),1);
    if spank < 3
    group1_temp = zeros(size(group_1_data,1),size(group_1_data,spank));
    group2_temp = zeros(size(group_2_data,1),size(group_2_data,spank));
    elseif spank == 3
    group1_temp = zeros(size(group_1_data,1)*size(group_1_data,2),size(group_1_data,spank));
    group2_temp = zeros(size(group_2_data,1)*size(group_2_data,2),size(group_2_data,spank));
    end
    for loop = 1:length(boop)
        if spank < 3
        group1_temp(loop,:) = group_1_data(boop(loop),:);
        group2_temp(loop,:) = group_2_data(boop(loop),:);
        elseif spank == 3
        beep = randi(size(group_1_data,2),size(group_1_data,2),1);
            for leep = 1:length(beep)
            group1_temp(leep+(size(group_1_data,2)*(loop-1)),:) = reshape(group_1_data(boop(loop),beep(leep),:),[1,size(group_1_data,spank)]);
            group2_temp(leep+(size(group_2_data,2)*(loop-1)),:) = reshape(group_2_data(boop(loop),beep(leep),:),[1,size(group_2_data,spank)]);
            end
        end
    end
    if log_scale == 1
    bootstrapped_diff(boot,:) = mean(group1_temp,1)./mean(group2_temp,1);
    else
    bootstrapped_diff(boot,:) = mean(group1_temp,1) - mean(group2_temp,1);
    end
end

else
for boot = 1:B_L
    if isempty(group_1_counter) && isempty(group_2_counter)
        boop_1 = randi(size(group_1_data,1),size(group_1_data,1),1);
        boop_2 = randi(size(group_2_data,1),size(group_2_data,1),1);
    if spank < 3
    group1_temp = zeros(size(group_1_data,1),size(group_1_data,spank));
    group2_temp = zeros(size(group_2_data,1),size(group_2_data,spank));
    elseif spank == 3
    group1_temp = zeros(size(group_1_data,1)*size(group_1_data,2),size(group_1_data,spank));
    group2_temp = zeros(size(group_2_data,1)*size(group_2_data,2),size(group_2_data,spank));
    end
    for loop_1 = 1:length(boop_1)
        if spank < 3
        group1_temp(loop_1,:) = group_1_data(boop_1(loop_1),:);
        elseif spank == 3
        beep_1 = randi(size(group_1_data,2),size(group_1_data,2),1);
            for leep_1 = 1:length(beep_1)
            group1_temp(leep_1+(size(group_1_data,2)*(loop_1-1)),:) = reshape(group_1_data(boop_1(loop_1),beep_1(leep_1),:),[1,size(group_1_data,spank)]);
            end
        end
    end
    for loop_2 = 1:length(boop_2)
        if spank < 3
        group2_temp(loop_2,:) = group_2_data(boop_2(loop_2),:);
        elseif spank == 3
        beep_2 = randi(size(group_2_data,2),size(group_2_data,2),1);
            for leep_2 = 1:length(beep_2)
            group2_temp(leep_2+(size(group_2_data,2)*(loop_2-1)),:) = reshape(group_2_data(boop_2(loop_2),beep_2(leep_2),:),[1,size(group_2_data,spank)]);
            end
        end
    end
    else
    boop_1 = randi(length(group_1_counter),length(group_1_counter),1);
    boop_2 = randi(length(group_2_counter),length(group_2_counter),1);
    if size(group_1_data,1) > length(group_1_counter) && size(group_2_data,1) > length(group_2_counter)
    group1_temp = [];
    group2_temp = [];
    else
        group1_temp = zeros(length(group_1_counter),size(group_1_data,spank));
        group2_temp = zeros(length(group_2_counter),size(group_2_data,spank));
    end
    for loop_1 = 1:length(boop_1)
        if size(group_1_data,1) > length(group_1_counter) 
        num_trials_animal = group_1_counter(boop_1(loop_1));
        previous_trials = sum(group_1_counter(1:boop_1(loop_1)-1));
        data_list = group_1_data(previous_trials+1:previous_trials+num_trials_animal,:);
        beep_1 = randi(num_trials_animal,num_trials_animal,1);
%            for leep_1 = 1:length(beep_1)  
                %group1_temp = [group1_temp; data_list(beep_1(leep_1),:)];
                group1_temp = [group1_temp; data_list(beep_1,:)];
 %           end
        else
        group1_temp(loop_1,:) = group_1_data(boop_1(loop_1),:);
        end
    end
    for loop_2 = 1:length(boop_2)
        if size(group_2_data,1) > length(group_2_counter) 
        num_trials_animal = group_2_counter(boop_2(loop_2));
        previous_trials = sum(group_2_counter(1:boop_2(loop_2)-1));
        data_list = group_2_data(previous_trials+1:previous_trials+num_trials_animal,:);
        beep_2 = randi(num_trials_animal,num_trials_animal,1);
           % for leep_2 = 1:length(beep_2)  
                %group2_temp = [group2_temp; data_list(beep_2(leep_2),:)];
                group2_temp = [group2_temp; data_list(beep_2,:)];
           % end
        else
        group2_temp(loop_2,:) = group_2_data(boop_2(loop_2),:);
        end
    end
    end
    if log_scale == 1
    bootstrapped_diff(boot,:) = mean(group1_temp,1)./mean(group2_temp,1);
    else
    bootstrapped_diff(boot,:) = mean(group1_temp,1) - mean(group2_temp,1);
    end
end    
end


sorted_boot = sort(bootstrapped_diff,1,'ascend');
figure(fig_num)
if isempty(subplot_info) == false
subplot(1,subplot_info(1),subplot_info(2))
else
end
color = [0.6 0.6 0.6];
x_axis = frequencies;
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
end

