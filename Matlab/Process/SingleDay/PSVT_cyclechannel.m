% Creates separate figure for each stimulus and separate subplot for each channel

clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Locust/1%/MAT_files/';
date = '07_09_2021';
position = 2;
tetrode = 2;
trial = 4;

time = [-2 10];

filepath = ['Position_' num2str(position)];
files = dir([readpath '/' date '/' filepath]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));

Colors = [242,147,147 ; 189,67,67 ; 237,70,47 ; 165,245,167 ; 57,123,39 ;...
    69,148,39 ; 79,201,251 ; 49,74,251 ; 147,147,147 ; 38,38,38]./255;
    
for cycle_stimuli = 10 %1:numel(files)
    load([readpath '/' date '/' filepath '/' files(cycle_stimuli).name]);
    data = eval([files(cycle_stimuli).name(1:end-4) '_data_filt']);
    data = squeeze(data(tetrode,:,trial,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate));
    
%     data = squeeze(data(tetrode,:,trial,:));

    PSVT_fig = figure('Units','normalized','Position',[0 0 0.75 1]);
    for cycle_channels = 1:size(data,1)
        subplot(size(data,1),1,cycle_channels);
        plot(1:size(data,2),data(cycle_channels,:),'Color',Colors(cycle_stimuli,:));
        
        xticks(0:sample_rate:(time(2)-time(1))*sample_rate)
        xticklabels(time(1):time(2));
        xlim([0 size(data,2)])
        xline(-time(1)*sample_rate,'--')

        ylim([min(data,[],'all') max(data,[],'all')])
        title(files(cycle_stimuli).name(1:end-4))
    end
end