% Adds new cells to master file and produces raster figure
tic; clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables below this line
Date = '08_20_2021';
Position = 1;
Tetrode = 2;
Cell = 3;

time = [-2 8];  % time relative to stimulus onset
PSTH_smooth = 3;
master_matrix = 'O_SS_Matrix.mat';      % master file in which to add/concatenate the files
%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables above this line

filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
master_path = 'Master_files';
mat_path = 'MAT_files';         % parent folder for mat_path2
mat_path2 = 'Odors';            % path extension for matlab workspace files (pre IGOR)
filepath2 = 'Spike_Sorted';     % parent folder for txt_path and fig_path
txt_path = 'Txt_files';         % path extension for text files (post IGOR)
fig_path = 'Figures';           % writes raster and PSTH figs to this path

Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
    
mat_file = [filepath '/' mat_path '/' mat_path2 '/' Date '/Position_' num2str(Position) '/Position_' num2str(Position)];
load([mat_file '.mat'], 'stim_on', 'stim_off', 'stim_on_trim', 'stim_off_trim', 'total_time', 'trim_time', 'bin_size', 'sample_rate', 'total_time');    % loads workspace file to identify sample rate, stimulus onset, stimulus offset, and total recording time in seconds
read_txt_path = dir([filepath '/' filepath2 '/' txt_path '/' Date '/Position_' num2str(Position) '/Tetrode_' num2str(Tetrode) '/*C' num2str(Cell) '.txt']);   % returns directory of all txt files pertaining to a specific cell

date = datestr(datetime(Date, 'InputFormat','MM_dd_yyyy'));
var_name = [date(4:6) '_' date(1:2) '_' date(10:11) '_Position_' num2str(Position)];
bins_per_sec = 1000/bin_size;

if strlength(num2str(Cell)) == 1                 % ensures all file names are same length
    newCell = ['0' num2str(Cell)];
elseif strlength(num2str(Cell)) == 2
    newCell = num2str(Cell);
end

experiment = [var_name '_C' newCell];
exp_path = [var_name '_C' newCell];

raster_fig = figure('Position', [40 0 600 1000]);
PSTH_fig = figure('Position', [40 0 600 1000]);
% imagesc_fig = figure('Position', [40 0 600 1000]);
for cycle_odors = 1:length(read_txt_path)
    odor = [read_txt_path(cycle_odors).name(1:end-7)];    % returns cell array of odorants
    spike_times = readmatrix([read_txt_path(cycle_odors).folder '/' read_txt_path(cycle_odors).name]);  % reads in raw spike times from txt files
    spikes = [];
    for cycle_trials = 1:size(spike_times,2)        % bins spike times and counts number of spikes per bin
        spikes_temp = nonzeros(spike_times(:,cycle_trials));
        spikes_temp = spikes_temp(spikes_temp >= stim_on+time(1) & spikes_temp <= stim_on+time(2))-(stim_on+time(1));
        
        spikes_bin(1,cycle_trials,:) = histcounts(spikes_temp,'BinLimits',[0 abs(time(1))+abs(time(2))],'BinWidth',bin_size/1000);

        if cycle_trials ~= 1 && size(spikes,1) < size(spikes_temp,1)
            spikes = padarray(spikes,size(spikes_temp,1)-size(spikes,1),'post');
        elseif cycle_trials ~= 1 && size(spikes,1) > size(spikes_temp,1)
            spikes_temp = padarray(spikes_temp,size(spikes,1)-size(spikes_temp,1),'post');
        end
        spikes = [spikes spikes_temp];
    end

    figure(raster_fig)
    subplot(length(read_txt_path),1,cycle_odors); hold on;   
    plot_Raster(spikes,time,stim_on,stim_off,odor,Cell,cycle_odors,Colors,date); hold off;
    figure(PSTH_fig)
    subplot(length(read_txt_path),1,cycle_odors);   
    plot_PSTH(spikes_bin,time,stim_on,stim_off,odor,Cell,cycle_odors,Colors,date,PSTH_smooth); hold off;   
    
    var_name = [odor '_spike_times'];
    eval([var_name '= spikes;']);

    if ~(exist([filepath '/' filepath2 '/' mat_path2 '/Trim_' txt_path '/' exp_path])==7)
        mkdir([filepath '/' filepath2 '/' mat_path2 '/Trim_' txt_path '/' exp_path]);
        save([filepath '/' filepath2 '/' mat_path2 '/Trim_' txt_path '/' exp_path '/' odor '.txt'], [odor '_spike_times']);        % save data to new master file
        fprintf([odor ' cell ' num2str(Cell) ' added to text file'])       
    elseif ~(exist([filepath '/' filepath2 '/' mat_path2 '/Trim_' txt_path '/' exp_path '/' odor '.txt'])==2)
        save([filepath '/' filepath2 '/' mat_path2 '/Trim_' txt_path '/' exp_path '/' odor '.txt'], [odor '_spike_times']);        % save data to new master file
        fprintf([odor ' cell ' num2str(Cell) ' added to text file'])       
    else
        fprintf([odor ' cell ' num2str(Cell) 'text file previously saved.\n'])
    end
    
    if ~(exist([filepath '/' master_path '/' master_matrix])==2)
        var_name = [odor '_spikes_bin'];
        eval([var_name '= spikes_bin;']); 
    
        save([filepath '/' master_path '/' master_matrix], [odor '_spikes_bin'], 'experiment');        % save data to new master file
        fprintf([' and master matrix!\n'])
    elseif ~ismember([odor '_spikes_bin'],who('-file',[filepath '/' master_path '/' master_matrix]))
        var_name = [odor '_spikes_bin'];
        eval([var_name '= spikes_bin;']); 
    
        save([filepath '/' master_path '/' master_matrix], [odor '_spikes_bin'], 'experiment', '-append');        % save data to new master file 
        fprintf([' and master matrix!\n'])
    else
        load([filepath '/' master_path '/' master_matrix],[odor '_spikes_bin'],'experiment');
        spikes = vertcat(eval([odor '_spikes_bin']),spikes_bin);
        if cycle_odors == 1
            experiment = vertcat(experiment,{exp_path});
        end
        var_name = [odor '_spikes_bin'];
        eval([var_name '= spikes;']);
        save([filepath '/' master_path '/' master_matrix], [odor '_spikes_bin'], 'experiment', '-append');        % save data to new master file 
        fprintf([' and master matrix!\n'])
    end
end

if ~(exist([filepath '/' filepath2 '/' mat_path2 '/' fig_path])==7)
    mkdir([filepath '/' filepath2 '/' mat_path2 '/' fig_path])
end
if ~(exist([filepath '/' filepath2 '/' mat_path2 '/' fig_path '/' exp_path '_Raster.fig'])==2) && ~(exist([filepath '/' filepath2 '/' mat_path2 '/' fig_path '/' exp_path '_PSTH.fig'])==2)
    fprintf('Saving Raster and PSTH figures.\n')
    saveas(raster_fig,[filepath '/' filepath2 '/' mat_path2 '/' fig_path '/' exp_path '_Raster.fig']);
    saveas(PSTH_fig,[filepath '/' filepath2 '/' mat_path2 '/' fig_path '/' exp_path '_PSTH.fig']);
else
    fprintf('Raster and PSTH figures previously saved.\n')
end
toc;