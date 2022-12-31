function SS_Raster(spikes,stim_on,stim_off,time,odor,Cell)
hold on;
for cycle_trials = 1:size(spikes,2)        % bins spike times and counts number of spikes per bin
    spikes_trial = repmat(nonzeros(spikes(:,cycle_trials))',2,1);
    trial = [cycle_trials-1; cycle_trials]; %Y offset for raster plot
    plot(spikes_trial, trial, 'k')  
end

title(gca,[odor ': Cell ' num2str(Cell)]);

xlim([0 (time(2)-time(1))]);
xticks(0:(time(2)-time(1)));
xticklabels(0:(time(2)-time(1)));
ylim([0 size(spikes,2)]);
yticks(0.5:1:size(spikes,2)-.5);
yticklabels(1:1:size(spikes,2));
 
stim_x = [stim_on stim_on stim_off stim_off];
stim_y = [0 size(spikes,2) size(spikes,2) 0]; 
stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);
hold off;
return