clear; close all; clc; tic;
% week = 3;
% day = '07_22_2021';

rms_window = 500;
smooth_window = 500;
bin_size = 10;  % in msecs
time = [0.5 0.75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data';
readpath = 'RHD_files';

if exist('day','var')
    positions = dir([basepath '/' readpath '/' day]);
    positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
    for cycle_positions = 1:numel(positions)
        exp_path = [day '/' positions(cycle_positions).name];
        stimuli = dir([basepath '/' readpath '/' exp_path]);
        stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
        for cycle_stimuli = 1:numel(stimuli)
            stimulus = stimuli(cycle_stimuli).name;
            if exist([basepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                load([basepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                data_filt = eval([stimulus(1:end-3) '_data_filt']);
            else
                [data_filt, stim_on, stim_off, sample_rate] = read_intan(stimulus, basepath, readpath, exp_path);
            end
            mat_to_igor(basepath, stimulus, day, cycle_positions, data_filt, sample_rate)
            rms(stimulus, data_filt, stim_on, stim_off, sample_rate, rms_window, smooth_window, bin_size, time, basepath, exp_path)
        end
    end
elseif ~exist('day','var')
    days = dir([basepath '/' readpath]);
    days = days(~ismember({days.name},{'.','..','.DS_Store'}));
    for day = 1:numel(days)
        positions = dir([days(day).folder '/' days(day).name]);
        positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
        for cycle_positions = 1:numel(positions)
            exp_path = [days(day).name '/' positions(cycle_positions).name];
            stimuli = dir([positions(cycle_positions).folder '/' positions(cycle_positions).name]);
            stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
            for cycle_stimuli = 1:numel(stimuli)
                fprintf(['Processing ' days(day).name ' ' positions(cycle_positions).name '\n'])
                stimulus = stimuli(cycle_stimuli).name;
                if exist([basepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                    load([basepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                    data_filt = eval([stimulus(1:end-3) '_data_filt']);
                    exp_pars = [stim_on, stim_off, sample_rate];
                else
                    [data_filt, stim_on, stim_off, sample_rate] = read_intan(stimulus, basepath, readpath, exp_path);
                end
                mat_to_igor(basepath, stimulus, days(day).name, cycle_positions, data_filt, sample_rate)
                rms(stimulus, data_filt, stim_on, stim_off, sample_rate, rms_window, smooth_window, bin_size, time, basepath, exp_path)
            end
        end
    end
    toc;
end
toc;
