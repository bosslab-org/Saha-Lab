clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Honeybee_1%/Master_files';

RMS_window = 50;
smooth_window = 50;
new_bin_size = 50;
bin_size = 10;

position = 1;
tetrode = 1;
time = [0.25 2];

filename = ['R' num2str(RMS_window) '_S' num2str(RMS_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2)) '.mat']; 
load([readpath '/' filename]);
stimuli = who('*_RMS');
bin_mult = new_bin_size/bin_size;

Colors = [242,147,147 ; 189,67,67 ; 237,70,47 ; 165,245,167 ; 57,123,39 ;...
    69,148,39 ; 79,201,251 ; 49,74,251 ; 147,147,147 ; 38,38,38]./255;

PSVT_fig = figure('Units','normalized','Position',[0 0 0.75 1]); hold on;
for cycle_stimuli = 1:numel(stimuli)
    data = eval(stimuli{cycle_stimuli});
    data = permute(mean(reshape(data,size(data,1),size(data,2),bin_mult,size(data,3)/bin_mult),[1 3]),[4 2 1 3]);
    
    plots(cycle_stimuli) = subplot(ceil(numel(stimuli)/2),2,cycle_stimuli);
    plot(1:size(data,1),data,'Color',[Colors(cycle_stimuli,:) 0.75]);
    
    title(stimuli{cycle_stimuli})
    
    xtick_res = 0.5;
    xticks(0:1000/new_bin_size*xtick_res:size(data,1))
    
    if cycle_stimuli == numel(stimuli)-1 || cycle_stimuli == numel(stimuli)
        xticklabels(time(1):xtick_res:time(2))
    else
        xticklabels([]);
    end
        
    
%     linkaxes(plots,'y') % Uniforms y-scale between graphs
end

PSVT_axes = axes(PSVT_fig,'visible','off');
PSVT_axes.YLabel.Visible='on';
ylabel(PSVT_axes, 'Voltage (µV)', 'FontSize', 16);
PSVT_axes.XLabel.Visible='on';
xlabel(PSVT_axes, 'Time (s)', 'FontSize', 16);

% xlim([0 (time(2)-time(1))*sample_rate])
% 
% stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
% stim_y = [min(ylim_min(cycle_stimuli)) max(ylim_max(cycle_stimuli)) max(ylim_max(cycle_stimuli)) min(ylim_min(cycle_stimuli))];
% stim_patch(cycle_stimuli) = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);    