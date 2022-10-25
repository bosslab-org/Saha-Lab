clear; close all; clc; tic;
filepath = '/Users/alexanderfarnum/Documents/MATLAB/Locust/Odorants';
readpath = 'RHD_files';
% week = '03';
% day = '04_15_2021';

trial_params = 2;
% trial_params = 1: stim_on = 5, stim_off = 9, trimmed_time = 19, num_samples = 400000 
% trial_params = 2: stim_on = 10, stim_off = 14, trimmed_time = 33, num_samples = 680000 
% trial_params = 3: stim_on = 10, stim_off = 14, trimmed_time = 36, num_samples = 740000 
% trial_params = 4: stim_on = 15, stim_off = 19, trimmed_time = 40, num_samples = 820000
% trial_params = 5: stim_on = 20, stim_off = 24, trimmed_time = 45, num_samples = 920000 

rms_window = 500;
smooth_window = 500;
bin_size = 10; %bin size (msecs)
time = [0 1];

if exist('week','var') && exist('day','var')
    positions = dir([filepath '/' readpath '/Week' week '/' day]);
    positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
    for cycle_positions = 1:numel(positions)
        exp_path = ['Week' week '/' day '/' positions(cycle_positions).name];
        stimuli = dir([filepath '/' readpath '/' exp_path]);
        stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
        for cycle_stimuli = 1:numel(stimuli)
            stimulus = stimuli(cycle_stimuli).name;
            [data_temp, exp_pars] = read_intan_2021(stimulus, filepath, readpath, exp_path, trial_params);
            rms_construct_2021(stimulus, data_temp, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)
        end
    end
    
elseif ~exist('week','var') && exist('day','var')
    positions = dir([filepath '/' readpath '/' day]);
    positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
    for cycle_positions = 1:numel(positions)
        exp_path = [day '/' positions(cycle_positions).name];
        stimuli = dir([filepath '/' readpath '/' exp_path]);
        stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
        for cycle_stimuli = 1:numel(stimuli)
            stimulus = stimuli(cycle_stimuli).name;
            [data_temp, exp_pars] = read_intan_2021(stimulus, filepath, readpath, exp_path, trial_params);
            
%             PSVT_test(data_temp, exp_pars)
% 
%             resume = 'Continue?  ';
%             resume = input(resume,'s');
% 
%             if resume == 'y'
%                 close all;
                rms_construct_2021(stimulus, data_temp, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)
%             else
%                 break
%             end
        end
    end
    
elseif ~exist('week','var') && ~exist('day','var')
    experiments = dir([filepath '/' readpath]);
    experiments = experiments(~ismember({experiments.name},{'.','..','.DS_Store'}));
    for cycle_dates = 1:numel(experiments)
        positions = dir([experiments(cycle_dates).folder '/' experiments(cycle_dates).name]);
        positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));

        if any(cycle_dates == 1:5)
            trial_params = 2;
        else
            trial_params = 4;
        end

        for cycle_positions = 1:numel(positions)
            exp_path = [experiments(cycle_dates).name '/' positions(cycle_positions).name];
            stimuli = dir([positions(cycle_positions).folder '/' positions(cycle_positions).name]);
            stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
            for cycle_stimuli = 1:numel(stimuli)
                stimulus = stimuli(cycle_stimuli).name;
                fprintf(['Processing ' experiments(cycle_dates).name '\n'])
                [data_temp, exp_pars] = read_intan_2021(stimulus, filepath, readpath, exp_path, trial_params);
                rms_construct_2021(stimulus, data_temp, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)  
            end
        end
    end
end
toc;