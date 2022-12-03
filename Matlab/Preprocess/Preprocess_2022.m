clear; close all; clc; tic;
exp_type = 'le';     % l1: Locust 1% biomarkers  
                     % lc: Locust cell culture
                     % le: Locust endometriosis
                     % h1: Honeybee 1% biomarkers  
                     % hb: Honeybee breath mixture

% week = 3;         % Locust Cell Culture ONLY
% day = '08_13_2021';

rms_window = 500;
smooth_window = 500;
bin_size = 10;  % in msecs
time = [0.25 2.25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = '/Users/alexanderfarnum/Documents/MATLAB';
readpath = 'RHD_files';

if strcmp(exp_type,'l1')
    writepath = [basepath '/Locust/1%'];
elseif strcmp(exp_type,'lc')
    writepath = [basepath '/Locust/CellCulture'];
elseif strcmp(exp_type,'le')
    writepath = [basepath '/Locust/Endometriosis'];
elseif strcmp(exp_type,'h1')
    writepath = [basepath '/Honeybee/1%']; 
elseif strcmp(exp_type,'hb')
    writepath = [basepath '/Honeybee/Breath'];
    % only include specific stimuli: idx = [1, 2, 4, 5, 9, 10]
end

if exist('day','var')
    if strcmp(exp_type,'l1') || strcmp(exp_type,'lc') || strcmp(exp_type,'le')
        positions = dir([writepath '/' readpath '/' day]);
        positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
        for cycle_positions = 1:numel(positions)
            exp_path = [day '/' positions(cycle_positions).name];
            stimuli = dir([writepath '/' readpath '/' exp_path]);
            stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
            for cycle_stimuli = 1:numel(stimuli)
                stimulus = stimuli(cycle_stimuli).name;
                if exist([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                    load([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                    data_filt = eval([stimulus(1:end-3) '_data_filt']);
                else
                    [data_filt, stim_on, stim_off, sample_rate] = read_intan_2022(stimulus, writepath, readpath, exp_path);
                end
                mat_to_igor(writepath, stimulus, day, cycle_positions, data_filt, sample_rate)
                rms_2022(stimulus, data_filt, stim_on, stim_off, sample_rate, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
            end
        end
    elseif strcmp(exp_type,'h1') || strcmp(exp_type,'hb')
        bees = dir([writepath '/' readpath '/' day]);
        bees = bees(~ismember({bees.name},{'.','..','.DS_Store'}));
        for cycle_bees = 1:numel(bees)
            positions = dir([bees(cycle_bees).folder '/' bees(cycle_bees).name]);
            positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
            for cycle_positions = 1:numel(positions)
                exp_path = [day '/' bees(cycle_bees).name '/' positions(cycle_positions).name];
                stimuli = dir([writepath '/' readpath '/' exp_path]);
                stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
                if strcmp(exp_type,'hb')
                    order = [1, 2, 4, 5, 9, 10];
                elseif strcmp(exp_type,'h1')
                    order = 1:numel(stimuli);
                end
                for cycle_stimuli = order
                    stimulus = stimuli(cycle_stimuli).name;
                    if exist([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                        load([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                        data_filt = eval([stimulus(1:end-3) '_data_filt']);
                        exp_pars = [stim_on, stim_off, sample_rate];
                    else
                        [data, data_filt, exp_pars] = read_intan_artifact_remove(stimulus, writepath, readpath, exp_path);
                    end
                    mat_to_igor(writepath, stimulus, day, cycle_positions, data_filt, sample_rate)
                    rms_2022(stimulus, data_filt, exp_pars, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
                    rms_2022(stimulus, data_filt, stim_on, stim_off, sample_rate, rms_window, smooth_window, bin_size, time, filepath, exp_path, exp_type)
                end
            end
        end
    end

elseif ~exist('day','var')
    days = dir([writepath '/' readpath]);
    days = days(~ismember({days.name},{'.','..','.DS_Store'}));
    for day = 1:numel(days)
        if strcmp(exp_type,'l1') || strcmp(exp_type,'lc')  || strcmp(exp_type,'le')
                positions = dir([days(day).folder '/' days(day).name]);
                positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
                for cycle_positions = 1:numel(positions)
                    exp_path = [days(day).name '/' positions(cycle_positions).name];
                    stimuli = dir([positions(cycle_positions).folder '/' positions(cycle_positions).name]);
                    stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
                    for cycle_stimuli = 1:numel(stimuli)
                        fprintf(['Processing ' days(day).name ' ' positions(cycle_positions).name '\n'])
                        stimulus = stimuli(cycle_stimuli).name;
                        if exist([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                            load([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                            data_filt = eval([stimulus(1:end-3) '_data_filt']);
                            exp_pars = [stim_on, stim_off, sample_rate];
                        else
                            if strcmp(exp_type,'lc') || strcmp(exp_type,'le')
                                [data_filt, stim_on, stim_off, sample_rate] = read_intan_2022(stimulus, writepath, readpath, exp_path);
                            elseif strcmp(exp_type,'l1')
                                if strcmp(days(day).name,'04_12_2021') ||...
                                        strcmp(days(day).name,'04_15_2021') ||...
                                        strcmp(days(day).name,'04_16_2021') ||...
                                        strcmp(days(day).name,'04_26_2021')
                                    trial = 1;
                                else
                                    trial = 2;
                                end
                                [data_filt, exp_pars] = read_intan_2021v2(stimulus, writepath, readpath, exp_path, trial);
                            end
                        end
                        mat_to_igor(writepath, stimulus, day, cycle_positions, data_filt, sample_rate)
%                         rms_2022(stimulus, data_filt, exp_pars, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
                        rms_2022(stimulus, data_filt, stim_on, stim_off, sample_rate, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
                    end
                end
    
        elseif strcmp(exp_type,'h1') || strcmp(exp_type,'hb')
            bees = dir([days(day).folder '/' days(day).name]);
            bees = bees(~ismember({bees.name},{'.','..','.DS_Store'}));
            for cycle_bees = 1:numel(bees)
                positions = dir([bees(cycle_bees).folder '/' bees(cycle_bees).name]);
                positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
                for cycle_positions = 1:numel(positions)
                    exp_path = [days(day).name '/' bees(cycle_bees).name '/' positions(cycle_positions).name];
                    stimuli = dir([positions(cycle_positions).folder '/' positions(cycle_positions).name]);
                    stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
                    if strcmp(exp_type,'hb')
                        order = [1, 2, 4, 5, 9, 10];
                    elseif strcmp(exp_type,'h1')
                        order = 1:numel(stimuli);
                    end
                    for cycle_stimuli = order
                        stimulus = stimuli(cycle_stimuli).name;
                        if exist([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat'],'file')
                            load([writepath '/MAT_files/' exp_path '/' stimulus(1:end-3) '.mat']);
                            data_filt = eval([stimulus(1:end-3) '_data_filt']);
                            exp_pars = [stim_on, stim_off, sample_rate];
                        else
                            [data, data_filt, exp_pars] = read_intan_artifact_remove(stimulus, writepath, readpath, exp_path);
                        end
                        mat_to_igor(writepath, stimulus, day, cycle_positions, data_filt, sample_rate)
                        rms_2022(stimulus, data_filt, exp_pars, rms_window, smooth_window, bin_size, time, basepath, exp_path, exp_type)
                    end
                end
            end
        end
    end
    toc;
end
toc;
