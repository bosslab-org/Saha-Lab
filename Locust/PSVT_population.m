clear; clc; close all;

readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files';
filename = 'R500_S500_b10_-1to4.mat';
include_odors = 'n';    % y or n


Colors = [246,146,30 ; 255,0,255 ; 0,200,0 ; 150,150,150 ; 0,0,0 ; 255,0,0 ; 117,76,36]/255;

load([readpath '/' filename]);
stimuli = who('*data_master');

if strcmp(include_odors,'n')
    order = [6 1 2 3 5];
elseif strcmp(include_odors,'y')
    order = [6 1 2 3 5 4 7];
end

stimuli = stimuli(order,:);
Colors = Colors(order,:);

PSVT_fig = figure('Units','normalized','Position',[0 0 0.75 1]); hold on;
for cycle_stimuli = 1:numel(stimuli)  
    data = eval(stimuli{cycle_stimuli});
    data = squeeze(mean(data,[1 2]));
    
    PSVT_plot = plot(1:size(data,1),data,'Color',Colors(cycle_stimuli,:));
end

title(Per-stimulus Voltage Plot

title(stimuli{cycle_stimuli})
xlim([0 (time(2)-time(1))*sample_rate])

stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
stim_y = [min(ylim_min(cycle_stimuli)) max(ylim_max(cycle_stimuli)) max(ylim_max(cycle_stimuli)) min(ylim_min(cycle_stimuli))];
stim_patch(cycle_stimuli) = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);    

