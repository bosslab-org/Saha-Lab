clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Honeybee/LC/Master_files';

RMS_window = 500;
smooth_window = 500;
new_bin_size = 50;
bin_size = 10;

position = 1;
tetrode = 1;
time = [0 4];

filename = ['R' num2str(RMS_window) '_S' num2str(RMS_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2)) '.mat']; 
load([readpath '/' filename]);
stimuli = who('*_RMS');
bin_mult = new_bin_size/bin_size;

Colors = [242,147,147 ; 189,67,67 ; 237,70,47 ; 165,245,167 ; 57,123,39 ;...
    69,148,39 ; 79,201,251 ; 49,74,251 ; 147,147,147 ; 38,38,38]./255;

for cycle_stimuli = 1:2:numel(stimuli)
    data_filt = eval(stimuli{cycle_stimuli});
    data_filt = permute(mean(reshape(data_filt,size(data_filt,1),size(data_filt,2),bin_mult,size(data_filt,3)/bin_mult),[1 3]),[4 2 1 3]);
    
    data_process = eval(stimuli{cycle_stimuli+1});
    data_process = permute(mean(reshape(data_process,size(data_process,1),size(data_process,2),bin_mult,size(data_process,3)/bin_mult),[1 3]),[4 2 1 3]);
    
    xtick_res = 0.5;

    subplot(1,2,1);
    plot(1:size(data_filt,1),data_filt,'Color',[0 0 0 0.5]);
    title(stimuli{cycle_stimuli})
    xticks(0:1000/new_bin_size*xtick_res:size(data_filt,1))

    subplot(1,2,2);
    plot(1:size(data_process,1),data_process,'Color',[0 0 0 0.5]);
    title(stimuli{cycle_stimuli+1})
    xticks(0:1000/new_bin_size*xtick_res:size(data_filt,1))
end  