clear; close all; clc; tic;
exp_type = 'hb';     % l1: Locust 1% biomarkers  
                     % lc: Locust cell culture
                     % h1: Honeybee 1% biomarkers  
                     % hb: Honeybee breath mixture

% week = 3;         % Locust Cell Culture ONLY
day = '09_29_2022';

rms_window = 500;
smooth_window = 500;
bin_size = 10;  % in msecs
time = [0 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = '/Users/alexanderfarnum/Documents/MATLAB';
readpath = 'RHD_files';

if strcmp(exp_type,'l1')
    matpath = 'Locust/1%';
elseif strcmp(exp_type,'lc')
    matpath = 'Locust/CellCulture'; 
elseif strcmp(exp_type,'h1')
    matpath = 'Honeybee/1%'; 
elseif strcmp(exp_type,'hb')
    matpath = 'Honeybee/Breath';
    % only include specific stimuli: idx = [1, 2, 4, 5, 9, 10]
end

if exist('day','var')
    if strcmp(exp_type,'l1') || strcmp(exp_type,'lc')
        positions = dir([basepath '/' matpath '/' readpath '/' day]);
        positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
        for cycle_positions = 1:numel(positions)
            exp_path = [day '/' positions(cycle_positions).name];
            stimuli = dir([basepath '/' matpath '/' readpath '/' exp_path]);
            stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
            for cycle_stimuli = 1:numel(stimuli)
                stimulus = stimuli(cycle_stimuli).name;
                if exist([basepath '/' matpath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                    load([basepath '/' matpath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                    data_filt = eval([stimulus(1:end-3) '_data_filt']);
                    exp_pars = [stim_on, stim_off, sample_rate];
                else
                    [data_filt, exp_pars] = read_intan_2022(stimulus, [basepath '/' matpath], readpath, exp_path);
                end
                rms_construct_2022(stimulus, data_filt, exp_pars, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
            end
        end
    elseif strcmp(exp_type,'h1') || strcmp(exp_type,'hb')
        bees = dir([basepath '/' matpath '/' readpath '/' day]);
        bees = bees(~ismember({bees.name},{'.','..','.DS_Store'}));
        for cycle_bees = 1:numel(bees)
            positions = dir([bees(cycle_bees).folder '/' bees(cycle_bees).name]);
            positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
            for cycle_positions = 1:numel(positions)
                exp_path = [day '/' bees(cycle_bees).name '/' positions(cycle_positions).name];
                stimuli = dir([basepath '/' matpath '/' readpath '/' exp_path]);
                stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
                if strcmp(exp_type,'hb')
                    order = [1, 2, 4, 5, 9, 10];
                elseif strcmp(exp_type,'h1')
                    order = 1:numel(stimuli);
                end
                for cycle_stimuli = order
                    stimulus = stimuli(cycle_stimuli).name;
                    if exist([basepath '/' matpath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                        load([basepath '/' matpath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                        data_filt = eval([stimulus(1:end-3) '_data_filt']);
                        exp_pars = [stim_on, stim_off, sample_rate];
                    else
                        [data, data_filt, exp_pars] = read_intan_artifact_remove(stimulus, [basepath '/' matpath], readpath, exp_path);
                    end
                    rms_construct_2022(stimulus, data_filt, exp_pars, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
                end
            end
        end
    end

elseif ~exist('day','var')
    experiments = dir([basepath '/' readpath]);
    experiments = experiments(~ismember({experiments.name},{'.','..','.DS_Store'}));
    for cycle_dates = 1:numel(experiments)
        positions = dir([experiments(cycle_dates).folder '/' experiments(cycle_dates).name]);
        positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
        for cycle_positions = 1:numel(positions)
            exp_path = [experiments(cycle_dates).name '/' positions(cycle_positions).name];
            stimuli = dir([positions(cycle_positions).folder '/' positions(cycle_positions).name]);
            stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
            for cycle_stimuli = 1:numel(stimuli)
                fprintf(['Processing ' experiments(cycle_dates).name ' ' positions(cycle_positions).name '\n'])
                stimulus = stimuli(cycle_stimuli).name;
                if exist([basepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                    load([basepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                    data_filt = eval([stimulus(1:end-3) '_data_filt']);
                    exp_pars = [stim_on, stim_off, sample_rate];
                else
                    if strcmp(exp_type,'l1') || strcmp(exp_type,'lc')
                        [data_filt, exp_pars] = read_intan_2022(stimulus, [basepath '/' matpath], readpath, exp_path);
                    elseif strcmp(exp_type,'h1') || strcmp(exp_type,'hb')   
                        [data, data_filt, exp_pars] = read_intan_artifact_remove(stimulus, [basepath '/' matpath], readpath, exp_path);
                    end
                end
                rms_construct_2022(stimulus, data_filt, exp_pars, rms_window, smooth_window, bin_size, time, basepath, exp_path)
            end
        end
    end
    toc;
end
toc;
