clear; close all; clc; tic;
readpath = '/Users/Xander/Documents/MATLAB/Honeybee_LungCancer/Master_files';
date = '05_24_2022';
position = 1;
tetrode = 1;
channel = 1;

rms_window = 500;
smooth_window = 500;
bin_size = 10; %bin size (msecs)
time = [0 4];

filename = ['R' num2str(rms_window) '_S' num2str(smooth_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
load([readpath '/' filename '.mat']); 
files = who([readpath '/' filename '.mat'], '*_RMS');

date_temp = ismember(experiment,[date '/Position_' num2str(position) ' Tetrode_' num2str(tetrode)],'rows');

for cycle_stimuli = 8 %1:numel(files)    
    data = eval(files{cycle_stimuli});
    data = permute(data(date_temp,:,:),[3 2 1]);
        
    ylim_min = min(data,[],'all');
    ylim_max = max(data,[],'all');
    
    figure('Units','normalized','Position',[0 0 0.5 1]);  % generate single figure for all Odorant PCA trajectories
    for cycle_trials = 1:size(data,2)
        subplot(size(data,2),1,cycle_trials);
        plot(1:(time(2)-time(1))*1000/bin_size,data(:,cycle_trials));
    
        if cycle_trials ~= size(data,2)
            set(gca,'XTick',[])
        else
            xticks([0:1000/bin_size:(time(2)-time(1))*1000/bin_size])
            xticklabels(time(1):time(2));
        end

        ylim([ylim_min ylim_max]);
        stim_x = [abs(time(1))*1000/bin_size abs(time(1))*1000/bin_size (stim_off-stim_on+abs(time(1)))*1000/bin_size (stim_off-stim_on+abs(time(1)))*1000/bin_size];
        stim_y = [ylim_min ylim_max ylim_max ylim_min];
        patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);    
    end 
    sgtitle(files{cycle_stimuli})
end