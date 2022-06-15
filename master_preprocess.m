clear; close all; clc; tic;
filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
readpath = 'RHD_files';
week = 'Week11';
date = '08_12_2021';
trial_params = 4;
% trial_params = 1: stim_on = 5, stim_off = 9, trimmed_time = 19, num_samples = 400000 
% trial_params = 4: stim_on = 10, stim_off = 14, trimmed_time = 33, num_samples = 680000 
% trial_params = 4: stim_on = 10, stim_off = 14, trimmed_time = 36, num_samples = 740000 
% trial_params = 4: stim_on = 15, stim_off = 19, trimmed_time = 40, num_samples = 820000 
rms_window = 500; % number of samples
smooth_window = 500;
bin_size = 10; %bin size in msecs
time_window = [0 4];

positions = dir([filepath '/' readpath '/' week '/' date]);
positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));

for cycle_positions = 1:numel(positions)
    exp_path = [week '/' date '/' positions(cycle_positions).name];
    stimuli = dir([filepath '/' readpath '/' week '/' date '/' positions(cycle_positions).name]);
    stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
    for cycle_stimuli = 1:length(stimuli)
        stimulus = stimuli(cycle_stimuli).name;
        [data_temp, exp_pars] = read_intan(stimulus, filepath, readpath, exp_path, trial_params);
        rms_construct(stimulus, data_temp, exp_pars, rms_window, smooth_window, bin_size, time_window, filepath, exp_path)
    end
end
toc;

