function plot_PSTH(data,time,stimulus_on,stimulus_off,odor,Cell,cycle_odors,Colors,date,PSTH_smooth)
data = permute(data,[2 3 1]);
data_mean = mean(data,1);
spline_interpolates = linspace(0,size(data,2),size(data,2)*PSTH_smooth);    % calculate number of spline interpolates
spline_mean = spline(1:size(mean(data,1),2),mean(data,1),spline_interpolates);
spline_sem = spline(1:size(data,2),std(data)/sqrt(size(data,1)),spline_interpolates);

% mean_plot = plot(1:size(data,2),data_mean, 'LineWidth', 1.5, 'Color', Colors(cycle_odors,:));
hold on;
PSTH_smooth_handle(cycle_odors) = plot(1:numel(spline_mean),spline_mean, 'LineWidth', 1.5, 'Color', Colors(cycle_odors,:));
sem_patch = patch([1:size(data,2)*PSTH_smooth, fliplr(1:size(data,2)*PSTH_smooth)], [(spline_mean - spline_sem), fliplr(spline_mean + spline_sem)], Colors(cycle_odors,:),'FaceAlpha',0.5,'LineStyle','none');   

samples_per_sec = (size(data,2)*PSTH_smooth)/(time(2)-time(1));

sgtitle(date);
title(gca,[odor ': Cell ' num2str(Cell)]);
stim_x = [-time(1)*samples_per_sec -time(1)*samples_per_sec ((stimulus_off-stimulus_on)-time(1))*samples_per_sec ((stimulus_off-stimulus_on)-time(1))*samples_per_sec];
stim_y = [0 max(spline_mean+spline_sem,[],'all') max(spline_mean+spline_sem,[],'all') 0]; 
stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);

xticks(0:samples_per_sec:samples_per_sec*(time(2)-time(1)));
xticklabels(0:(time(2)-time(1)));

xlim([0 samples_per_sec*(time(2)-time(1))]);

% xlim([0 size(data,2)]);

ylim([0 max(spline_mean+spline_sem,[],'all')]);
hold off;
return