clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Honeybee_LungCancer/MAT_files/';
date = '05_24_2022';
position = 1;
tetrode = 1;
channel = 2;
time = [-1 6];

filepath = ['Position_' num2str(position)];
files = dir([readpath '/' date '/' filepath]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));

Colors = [242,147,147 ; 189,67,67 ; 237,70,47 ; 165,245,167 ; 57,123,39 ;...
    69,148,39 ; 79,201,251 ; 49,74,251 ; 147,147,147 ; 38,38,38]./255;
    
for cycle_stimuli = [7:8] %1:numel(files)
    load([readpath '/' date '/' filepath '/' files(cycle_stimuli).name]);
    data = eval([files(cycle_stimuli).name(1:end-4) '_data_filt(tetrode,channel,:,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)']);
    PSVT_fig = figure('Units','normalized','Position',[0 0 0.75 1]);
    for cycle_trials = 1:size(data,3)
        data_temp = squeeze(data(:,:,cycle_trials,:));
        subplot(size(data,3),1,cycle_trials);
        plot(1:(time(2)-time(1))*sample_rate,data_temp,'Color',Colors(cycle_stimuli,:)); 
        if cycle_trials ~= size(data,3)
            set(gca,'XTick',[])
        else
            xticks([0:sample_rate:(time(2)-time(1))*sample_rate])
            xticklabels(time(1):time(2));
        end
        xlim([0 (time(2)-time(1))*sample_rate])
        ylim([min(data,[],'all') max(data,[],'all')]);
        
        stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
        stim_y = [min(data,[],'all') max(data,[],'all') max(data,[],'all') min(data,[],'all')];
        stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);
    end
    sgtitle(files(cycle_stimuli).name(1:end-4))
    PSVT_axes = axes(PSVT_fig,'visible','off');
    PSVT_axes.YLabel.Visible='on';
    ylabel(PSVT_axes, 'Voltage (µV)', 'FontSize', 16);
    PSVT_axes.XLabel.Visible='on';
    xlabel(PSVT_axes, 'Time (s)', 'FontSize', 16);
end