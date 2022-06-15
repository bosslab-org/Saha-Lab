clear; clc; close all; tic;

file_type = 'r';       % ss (Spike Sorted) or r (RMS)
input_files = 'm';      % m (Multiple) or s (Single) file(s)
include_odors = 'y';    % y or n
daywise_comp = 'n'; % y or n

filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
bin_size = 10;
time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;          % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
PSTH_smooth = 1;            % PSTH smoothing factor
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 50;        % line from origin every x msecs, MUST be set to multiple of new bin_size

if strcmp(input_files,'m')
    if strcmp(include_odors,'n')
        order = [6 1 2 3 5];
    elseif strcmp(include_odors,'y')
        order = [6 1 2 3 5 4 7];
    end        
    Colors = [246,145,30 ; 255,0,255 ; 0,255,0 ; 150,150,150 ; 0,0,0 ; 255,0,0 ; 117,76,36]./255; 
    Colors = Colors(order,:);
    if strcmp(file_type,'r')
        readfile = '/RMS_files/R500_S500_b10';
%         readfile = [readfile '_' num2str(time(1)) 'to' num2str(time(2))];
    elseif strcmp(file_type,'ss') 
        readfile = 'SS_b10';
    end
    if strcmp(daywise_comp,'y')
        select_day = ['Day' num2str(input('Select day (1, 2, 3, 4): ','s'))];
        readfile = [readfile '_' select_day];
    end
    var_ext = '_RMS';
    load([filepath '/' readfile '.mat']);
    files = who([filepath '/' readfile '.mat'], ['*' var_ext]);
    files = files(order);

    bin_mult = new_bin_size/bin_size;
    bins_per_sec = 1000/bin_size;
    for cycle_stimuli = 1:numel(files)
        stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
        data_temp = eval(files{cycle_stimuli});       
        data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
    end    
end         
    
set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
title_time = [stim_on+time(1) stim_on+time(2)];
PCA_line_mult = PCA_line_mult/new_bin_size;
titles = [num2str(time(1)) ' to ' num2str(time(2)) ' seconds after stimulus onset'];
    

PCA_master_fig = PCA(data, Colors, stimuli, PCA_smooth, ['PCA (' file_type '): ' titles], time, bins_per_sec, floor(PCA_line_mult));        
LDA_master_fig = LDA(data, Colors, stimuli, ['LDA (' file_type '): ' titles]);

[bin_preds, trial_preds] = TestData(data);

Conf_fig_bin = Confusion(bin_preds, stimuli, ['Confusion (' file_type '): ' titles]);
Conf_fig_trial = Confusion(trial_preds, stimuli, ['Confusion (' file_type '): ' titles]);



%     PSTH_fig = PSTH(data,Colors,stimuli,time,bins_per_sec,['PSTH (' file_type '): ' titles]);
%     Del_PSTH_fig = Delta_PSTH(data,Colors,stimuli,time,bins_per_sec,['PSTH (' file_type '): ' titles]);


    
%%%%%%%%%%%%%%%%%%% Single day analysis %%%%%%%%%%%%%%%%%%%%%c
% elseif input_files == 's'
%     Date = 'Specify date (MM_DD_YYYY): ';
%     Date = input(Date,'s');
%     Position = 'Specify position: ';
%     Position = input(Position,'s');
%     Tetrode = 'Specify tetrode: ';
%     Tetrode = input(Tetrode,'s');
%         
%     if file_type == 'r'
%         if stimulus == 'o'
%             readpath2 = '/MAT_files/Odors';    
%             Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
%         elseif stimulus == 'c'
%             readpath2 = 'MAT_files/Cell_Culture';
%         end
%         Odorant = 'Specify odorant name or all: ';
%         Odorant = input(Odorant,'s');
%         files = dir([readpath '/' readpath2 '/' Date '/Position_' Position]);
%         if strcmp(Odorant, 'all')
%             files = files(~ismember({files.name},{'.','..',['Position_' Position '.mat']}));
%         else
%             files = files(ismember({files.name},{[Odorant '.mat']})); 
%         end
%         for cycle_odors = 1:numel(files)
%             Odorant = files(cycle_odors).name(1:end-4);
%             if cycle_odors == 1
%                 load([readpath '/' readpath2 '/' Date '/Position_' Position '/' files(cycle_odors).name],[Odorant '_data_filt'],'stim_on','stim_off','sample_rate','total_time');
%                 load([readpath '/' readpath2 '/' Date '/Position_' Position '/Position_' num2str(Position) '.mat'],[Odorant '_RMS_data_filt'],'stim_on_trim','stim_off_trim','bin_size','trim_time');
% 
%                 trim_time = [stim_on_trim+time(1) stim_on_trim+time(2)];
%                 time = [stim_on+time(1) stim_on+time(2)];
%                 
%                 bin_mult = new_bin_size/bin_size;
%                 bins_per_sec = 1000/(bin_size*bin_mult);
%             else
%                 load([readpath '/' readpath2 '/' Date '/Position_' Position '/' files(cycle_odors).name],[Odorant '_data_filt']);
%                 load([readpath '/' readpath2 '/' Date '/Position_' Position '/Position_' num2str(Position) '.mat'],[Odorant '_RMS_data_filt']);
%             end
%             odor_data = eval([Odorant '_data_filt']);
%             odor_RMS_data = eval([files(cycle_odors).name(1:end-4) '_RMS_data_filt']);
%             
%             odor_data = odor_data(str2double(Tetrode),:,:,time(1)*sample_rate+1:(time(2)*sample_rate));            
%             
%             odor_RMS_data_base = mean(odor_RMS_data(str2double(Tetrode),:,:,(stim_on_trim-2)*bins_per_sec+1:stim_on_trim*bins_per_sec),4);
%             odor_RMS_data_temp = odor_RMS_data(str2double(Tetrode),:,:,trim_time(1)*bins_per_sec+1:trim_time(2)*bins_per_sec);
%             odor_RMS_data_temp = odor_RMS_data_temp-odor_RMS_data_base;
% 
%             exp_pars = struct('Date',Date,'Odorant',Odorant,'Position',Position,'stim_on',stim_on,'stim_on_trim',stim_on_trim,'stim_off',stim_off,'stim_off_trim',stim_off_trim,'sample_rate',sample_rate,'bins_per_sec',bins_per_sec,'total_time',total_time,'trim_time',trim_time);
%             for cycle_chans = 1:size(odor_RMS_data_temp,2)            
%                 PSVT(exp_pars,odor_data(:,cycle_chans,:,:),str2double(Tetrode),cycle_chans,time)
%                 RMS_PSVT(exp_pars,odor_RMS_data_temp(:,cycle_chans,:,:),str2double(Tetrode),cycle_chans,trim_time)
%     %             RMS_PSTH(exp_pars,file,cycle_tetrodes,cycle_channels,time)
%     %             RMS_Raster(exp_pars,file,cycle_tetrodes,cycle_channels,time)
%             end
%         end
% 
%     elseif file_type == 'ss'
%         Cell = 'Specify cell: ';
%         Cell = input(Cell,'s');
%         readpath2 = 'Spike_Sorted';
%         if stimulus == 'o'
%             readpath2 = [readpath2 '/Odors/Trim_Txt_files'];
%             Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ; 255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
%         elseif stimulus == 'c'
%             readpath2 = [readpath2 '/Cell_Culture/Trim_Txt_files'];
%                 Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192]./255;           
%         end
%         
%         files = dir([readpath '/' readpath2 '/' Date '/Position_' Position '/Tetrode_' num2str(Tetrode) '/Cell' num2str(Cell) '/*.txt']);        
%         load([readpath '/MAT_files/Odors/' Date '/Position_' Position '/Position_' num2str(Position) '.mat'],'stim_on','stim_on_trim','stim_off','stim_off_trim','trim_time','bin_size');
%         set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
% 
%         trim_time = [stim_on_trim+time(1) stim_on_trim+time(2)];
%               
%         Raster_title = ['Raster: ' num2str(trim_time(1)-stim_on_trim) ' to ' num2str(trim_time(2)-stim_on_trim) ' seconds after stimulus onset'];
%         PSTH_title = ['PSTH: ' num2str(trim_time(1)-stim_on_trim) ' to ' num2str(trim_time(2)-stim_on_trim) ' seconds after stimulus onset'];
%   
%         raster_fig = figure('Position', [40 0 600 1000]);
%         PSTH_fig = figure('Position', [40 0 600 1000]);
%         for cycle_odors = 1:numel(files)
%             odor = files(cycle_odors).name(1:end-4);
%             spike_times = readmatrix([files(cycle_odors).folder '/' files(cycle_odors).name]);  % reads in raw spike times from txt files
%             spikes = [];
%             for cycle_trials = 1:size(spike_times,2)        % bins spike times and counts number of spikes per bin
%                 spikes_temp = nonzeros(spike_times(:,cycle_trials));
%                 spikes_temp = spikes_temp(spikes_temp >= trim_time(1) & spikes_temp <= trim_time(2));
%                 spikes_bin(cycle_trials,:) = histcounts(spikes_temp,'BinLimits',[trim_time(1) trim_time(2)],'BinWidth',bin_size/1000);
% 
%                 if cycle_trials ~= 1 && size(spikes,1) < size(spikes_temp,1)
%                     spikes = padarray(spikes,size(spikes_temp,1)-size(spikes,1),'post');
%                 elseif cycle_trials ~= 1 && size(spikes,1) > size(spikes_temp,1)
%                     spikes_temp = padarray(spikes_temp,size(spikes,1)-size(spikes_temp,1),'post');
%                 end
%                 spikes = [spikes spikes_temp];
%             end
%             figure(raster_fig)
%             Raster(cycle_odors) = subplot(numel(files),1,cycle_odors);   
%             plot_Raster_Master(spikes,time,stim_on_trim,stim_off_trim,odor,Cell,cycle_odors,numel(files),Colors,date,Raster_title);
% 
%             figure(PSTH_fig)
%             PSTH(cycle_odors) = subplot(numel(files),1,cycle_odors); 
%             [PSTH_max(cycle_odors), PSTH_patch(cycle_odors)] = plot_PSTH_Master(spikes_bin,time,stim_on_trim,stim_off_trim,odor,Cell,cycle_odors,numel(files),Colors,date,PSTH_smooth,1000/new_bin_size,PSTH_title);  
% 
%             var_name = [odor '_spike_times'];
%             eval([var_name '= spikes;']);
%         end
%         linkaxes(PSTH,'y');
%         PSTH(1,1).YLim = [0 max(PSTH_max)];
%         linkprop(PSTH_patch, 'YData');
%         PSTH_patch(1,1).YData = [0 max(PSTH_max) max(PSTH_max) 0];
%     end
% end
% toc;
% 
% save_figs = 'Save all figures (y / n)? ';
% save_figs = input(save_figs,'s');
% 
% if save_figs == 'y'   
%     writepath = [writepath '_' writefigs '/' num2str(time(1)) 'to' num2str(time(2)) '/Bin_size_' num2str(new_bin_size)];
%     if ~exist(writepath)
%         mkdir(writepath)
%     end
%     for cycle_figs = 1:numel(figs)
%         figure(figs{cycle_figs});
%         if ~exist([writepath '/' figs_titles{cycle_figs}])
%             saveas(figs{cycle_figs},[writepath '/' figs_titles{cycle_figs}])
%         end
%     end
%         
%     LDA_writepath = [writepath '/Norm' num2str(LDA_norm) '_' num2str(lambda) 'lambda'];
%     if ~exist(LDA_writepath)
%         mkdir(LDA_writepath)
%     end
%     for cycle_norm_figs = 1:numel(LDA_figs)
%         figure(LDA_figs{cycle_norm_figs});
%         if ~exist([LDA_writepath '/' LDA_figs_titles{cycle_norm_figs}])
%             saveas(LDA_figs{cycle_norm_figs},[LDA_writepath '/' LDA_figs_titles{cycle_norm_figs}])
%         end
%     end
%     
%     all_writepath = [writepath '/Norm' num2str(HD_norm)];
%     if ~exist(all_writepath)
%         mkdir(all_writepath)
%     end
%     for cycle_norm_figs = 1:numel(all_figs)
%         figure(all_figs{cycle_norm_figs});
%         if ~exist([all_writepath '/' all_figs_titles{cycle_norm_figs}])
%             saveas(all_figs{cycle_norm_figs},[all_writepath '/' all_figs_titles{cycle_norm_figs}])
%         end
%     end
%     close all;
% end
% toc;
%     
