clear; close all; clc; tic;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Locust/Endometriosis/MAT_files';
date = '10_20_2022';
position = 1;
tetrode = 1;
channel = 3;

time = [-2 6];
sample_rate = 20000;

files = dir([readpath '/' date '/Position_' num2str(position)]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));
for cycle_stimuli = 1:numel(files)
    load([files(cycle_stimuli).folder '/' files(cycle_stimuli).name]);
    data_filt = eval([[files(cycle_stimuli).name(1:end-4)] '_data_filt']);

    for cycle_trials = 1:size(data_filt,3)
        subplot(size(data_filt,3),1,cycle_trials);
        plot(1:size(data_filt,4),squeeze(data_filt(tetrode,channel,cycle_trials,:)),'Color',[0 0 0 0.5]);
%         plot((stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate,squeeze(data_filt(tetrode,channel,cycle_trials,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)),'Color',[0 0 0 0.5]);
        xticks(0:sample_rate:total_time*sample_rate)
        xticklabels(0:total_time)
    end
end