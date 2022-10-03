clear; close all; clc; tic;
filepath = '/Users/Xander/Documents/MATLAB/Honeybee/LC';
readpath = 'RHD_files';
% day = '05_24_2022';

rms_window = 500;
smooth_window = 500;
bin_size = 10;  % in msecs
time = [0 4];

if exist('day','var')
    bees = dir([filepath '/' readpath '/' day]);
    bees = bees(~ismember({bees.name},{'.','..','.DS_Store'}));
    for cycle_bees = 1:numel(bees)
        positions = dir([bees(cycle_bees).folder '/' bees(cycle_bees).name]);
        positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
        for cycle_positions = 1:numel(positions)
            exp_path = [day '/' bees(cycle_bees).name '/' positions(cycle_positions).name];
            stimuli = dir([filepath '/' readpath '/' exp_path]);
            stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
            for cycle_stimuli = 1:numel(stimuli)
                stimulus = stimuli(cycle_stimuli).name;

                if exist([filepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                    load([filepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                    data_filt = eval([stimulus(1:end-3) '_data_filt']);
                    data_processed = eval([stimulus(1:end-3) '_data_processed']);
                    exp_pars = [stim_on, stim_off, sample_rate];
                else
                    [data_filt, data_processed, exp_pars] = read_intan_tetrodes(stimulus, filepath, readpath, exp_path);
                end
                
                % Comment in to plot test PSTH
%                 PSTH_fig = test_PSTH(data, exp_pars, 3); pause(10); close(PSTH_fig);

                rms_construct_tetrodes(stimulus, data_filt, data_processed, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)
            end
        end
    end
    
elseif ~exist('day','var')
    experiments = dir([filepath '/' readpath]);
    experiments = experiments(~ismember({experiments.name},{'.','..','.DS_Store'}));
    for cycle_dates = 2:5 %1:numel(experiments)
        bees = dir([experiments(cycle_dates).folder '/' experiments(cycle_dates).name]);
        bees = bees(~ismember({bees.name},{'.','..','.DS_Store'}));
        for cycle_bees = 1:numel(bees)
            positions = dir([bees(cycle_bees).folder '/' bees(cycle_bees).name]);
            positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
            for cycle_positions = 1:numel(positions)
                exp_path = [experiments(cycle_dates).name '/' bees(cycle_bees).name '/' positions(cycle_positions).name];
                stimuli = dir([positions(cycle_positions).folder '/' positions(cycle_positions).name]);
                stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
                for cycle_stimuli = 1:numel(stimuli)
                    fprintf(['Processing ' experiments(cycle_dates).name ' ' bees(cycle_bees).name ' ' positions(cycle_positions).name '\n'])
                    stimulus = stimuli(cycle_stimuli).name;

                    if exist([filepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                        load([filepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                        data_filt = eval([stimulus(1:end-3) '_data_filt']);
                        data_processed = eval([stimulus(1:end-3) '_data_processed']);
                        exp_pars = [stim_on, stim_off, sample_rate];
                    else
                        [data_filt, data_processed, exp_pars] = read_intan_tetrodes(stimulus, filepath, readpath, exp_path);
                    end
                    
                    rms_construct_tetrodes(stimulus, data_filt, data_processed, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)
                end
            end
        end
    end
end
toc;
