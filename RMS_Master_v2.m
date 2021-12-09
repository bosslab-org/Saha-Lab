%% Master File for Processing RMS Data
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables below this line
filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings/';
filepath2 = 'Master_files';
master_file = 'Odor_RMS_Data_Filt';
load([filepath '/' filepath2 '/' master_file]);

time = [stimulus_on stimulus_off];
new_bin_size = 100;
PCA_smooth = 2;
PCA_line_mult = 500;      % line from origin every x msecs, MUST be set to multiple of new bin_size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables above this line

PCA_title = ['PCA: ' num2str(time(1)-stimulus_on) ' to ' num2str(time(2)-stimulus_on) ' seconds after stimulus onset'];
PCA_line_mult = PCA_line_mult/new_bin_size;
LDA_title = ['LDA: ' num2str(time(1)-stimulus_on) ' to ' num2str(time(2)-stimulus_on) ' seconds after stimulus onset'];

files = who('-file', [filepath '/' filepath2 '/' master_file]);
odors_data = contains(files, '_data_filt'); 
files = files(odors_data);
bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/(bin_size*bin_mult);

if numel(files) == 10 
    Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
end

master_ch_data = [];
master_data = [];
for cycle_odors = 1:numel(files)
    Odorants{cycle_odors} = files{cycle_odors}(1:end-10);
    odor_data = eval(files{cycle_odors});
    odor_data = odor_data(:,:,:,(time(1)*bins_per_sec)+1:time(2)*bins_per_sec);
    
    mean_ch_data = permute(mean(odor_data,2),[1 3 4 2]);
    mean_data = permute(mean(mean(odor_data,3),2),[4 1 2 3]);
    
    master_ch_data = cat(1,master_ch_data,mean_ch_data);
    master_data = vertcat(master_data,mean_data);   
end

PCA(master_data, Colors, Odorants, PCA_smooth, PCA_title, time, stimulus_on, bins_per_sec, PCA_line_mult)
LDA(master_data, Colors, Odorants, LDA_title)