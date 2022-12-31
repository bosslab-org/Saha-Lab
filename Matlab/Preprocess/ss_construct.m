clear; close all; clc; tic;
filepath = '/Users/Xander/Documents/MATLAB';
experiment = 'Honeybee_LungCancer';
readpath = 'Text_files';
writepath = 'SS_b10';

bin_size = 10;
old_stim_on = 19;
stim_on = 2;
stim_off = 6;
sample_rate = 20000;

stimuli = dir([filepath '/' experiment '/' readpath]);
stimuli = stimuli(~ismember({stimuli.name},{'.','..','.DS_Store'}));
for cycle_stimuli = 1:numel(stimuli)
    neurons = dir([filepath '/' experiment '/' readpath '/' stimuli(cycle_stimuli).name]);
    neurons = neurons(~ismember({neurons.name},{'.','..','.DS_Store'}));
    for cycle_neurons = 1:numel(neurons)
        neuron_data = load([filepath '/' experiment '/' readpath '/' stimuli(cycle_stimuli).name '/' neurons(cycle_neurons).name]);
        if nonzeros(neuron_data(end,:)) < 40
            time = 40;
        else
            time = 45;
        end
        for cycle_trials = 1:size(neuron_data,2)
            % Bin number of spikes according to specified bin_size
            spikes = histcounts(nonzeros(neuron_data(:,cycle_trials)),linspace(0,time,time*1000/bin_size+1));
            
            % Trim to 2 seconds pre-stim_on to 8 seconds post stim_on
            data(cycle_neurons,cycle_trials,:,cycle_stimuli) = spikes((old_stim_on-2)*(1000/bin_size)+1:(old_stim_on+8)*(1000/bin_size));
        end
    end
end

% Baseline Subtraction
% data = neurons, trials, samples, stimuli
data = data(:,:,(stim_on-2)*1000/bin_size+1:end,:) - ...
    mean(data(:,:,(stim_on-2)*1000/bin_size+1:stim_on*1000/bin_size,:),[2,3]);

for cycle_stimuli = 1:numel(stimuli)
    eval([stimuli(cycle_stimuli).name '_spikes = data(:,:,:,cycle_stimuli);'])
    if ~exist([filepath '/' experiment '/Master_files'],'dir')
        mkdir([filepath '/' experiment '/Master_files'])
    end
    if cycle_stimuli ~= 1
        save([filepath '/' experiment '/Master_files/' writepath], [stimuli(cycle_stimuli).name '_spikes'],'-append');
    elseif ~exist([filepath '/' experiment '/Master_files/' writepath],'file')
        fprintf('Creating new master file.\n')
        save([filepath '/' experiment '/Master_files/' writepath], [stimuli(cycle_stimuli).name '_spikes'],...
            'stim_on','stim_off','bin_size');
    end
end
