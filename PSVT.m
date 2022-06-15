clear; clc; close all;

date = '08_05_2021';
position = '1';
tetrode = 1;
channel = 3;
trial = 2;
time = [-1 6];

readpath = '/Users/Xander/Documents/MATLAB/Day_files';
filepath = [readpath '/' date '/Position_' position];
files = dir(filepath);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));
sample_rate = 20000;
for cycle_stimuli = 1:numel(files)
    file = files(cycle_stimuli).name(1:end-4);
    load([filepath '/' files(cycle_stimuli).name]);
    data = eval([file '_data_filt']);
    data = squeeze(data(tetrode,channel,trial,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate));
    figure;
    plot(1:(time(2)-time(1))*sample_rate,data)
    
    stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
    stim_y = [min(data) max(data) max(data) min(data)];
    stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);
        
    xlim([0 (time(2)-time(1))*sample_rate]);
    xticks([0:sample_rate:(time(2)-time(1))*sample_rate])
    xticklabels(time(1):time(2));
    xlabel('Time (s)', 'FontSize', 16);
    
    ylim([min(data) max(data)])
    ylabel('Voltage (µV)', 'FontSize', 16);
    
    title(file)
end