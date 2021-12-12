function plot_Raster(data,time,stim_on,stim_off,odor,Cell,cycle_odors,Colors,date)

for cycle_trials = 1:size(data,2)        % bins spike times and counts number of spikes per bin
    data_temp = repmat(data(:,cycle_trials)',2,1);
    trial = [cycle_trials-1; cycle_trials]; %Y offset for raster plot
    plot(data_temp, trial, 'Color', Colors(cycle_odors,:))  
end

sgtitle(date);
title(gca,[odor ': Cell ' num2str(Cell)]);

xlim([0 (time(2)-time(1))]);
xticks(0:(time(2)-time(1)));
xticklabels(0:(time(2)-time(1)));
ylim([0 size(data,2)]);
 
stim_x = [-time(1) -time(1) (stim_off-stim_on)-time(1) (stim_off-stim_on)-time(1)];
stim_y = [0 size(data,1) size(data,1) 0]; 
stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);
   
return