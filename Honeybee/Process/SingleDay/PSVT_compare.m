clear; close all; clc; tic;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Honeybee/LC/MAT_files';
date = '05_24_2022';
bee = 1;
position = 1;
tetrode = 1;
channel = 1;

rms_window = 500;
smooth_window = 500;
bin_size = 10; %bin size (msecs)
time = [0 4];

files = dir([readpath '/' date '/Bee_' num2str(bee) '/Position_' num2str(position)]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));
for cycle_stimuli = 1 %:numel(files)
    load([files(cycle_stimuli).folder '/' files(cycle_stimuli).name]);

    data_filt = eval([[files(cycle_stimuli).name(1:end-4)] '_data_filt']);
    data_processed = eval([[files(cycle_stimuli).name(1:end-4)] '_data_processed']);


    xtick_res = 0.5;

    subplot(1,2,1);
    for cycle_trials = 1:size(data_filt,3)
        plot(1:size(data_filt,1),squeeze(data_filt(1,channel,cycle_trials,:)),'Color',[0 0 0 0.5]);
    end
    title(stimuli{cycle_stimuli})
    xticks(0:1000/new_bin_size*xtick_res:size(data_filt,1))




    data_filt = permute(mean(reshape(data_filt,size(data_filt,1),size(data_filt,2),bin_mult,size(data_filt,3)/bin_mult),[1 3]),[4 2 1 3]);
    
    
    data_name = [files(cycle_stimuli).name(1:end-4) '_data_filt'];
    data = eval(data_name);
    
    % size(data) = (tetrodes,channels,trials,samples)
    data = preprocess_data(data(:,:,:,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)) -...
        mean(preprocess_data(data(:,:,:,(stim_on-2)*sample_rate+1:stim_on*sample_rate)),[3 4]);                
    data = permute(mean(reshape(permute(data,[4 3 1 2]),...
        (sample_rate*bin_size)/1000,size(data,4)/((sample_rate*bin_size)/1000),size(data,3),size(data,1),size(data,2)),[1 5]),[4 3 2 1]);  % takes mean of values within each bin and across tetrode channels
    % size(data) = (tetrodes,trials,binned_samples)
    
    data = permute(data(tetrode,:,:),[3 2 1]);
    
    ylim_min = min(data,[],'all');
    ylim_max = max(data,[],'all');
    
    figure('Units','normalized','Position',[0 0 0.5 1]);  % generate single figure for all Odorant PCA trajectories
    for cycle_trials = 1:size(data,2)
        subplot(size(data,2),1,cycle_trials);
        plot(1:(time(2)-time(1))*1000/bin_size,data(:,cycle_trials));
    
        if cycle_trials ~= size(data,2)
            set(gca,'XTick',[])
        else
            xticks([0:1000/bin_size:(time(2)-time(1))*1000/bin_size])
            xticklabels(time(1):time(2));
        end

        ylim([ylim_min ylim_max]);
        stim_x = [abs(time(1))*1000/bin_size abs(time(1))*1000/bin_size (stim_off-stim_on+abs(time(1)))*1000/bin_size (stim_off-stim_on+abs(time(1)))*1000/bin_size];
        stim_y = [ylim_min ylim_max ylim_max ylim_min];
        patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);    
    end 
    sgtitle(files(cycle_stimuli).name(1:end-4))
end