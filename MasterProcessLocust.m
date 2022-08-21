clear; clc; close all; tic;

file_type = 'r';        % s (Spike Sorted) or r (RMS)
include_odors = 'y';    % y or n
daywise_comp = 'n';     % y or n
save_figs = 'n';        % y or n

time = [0 4];               % time relative to stimulus onset
new_bin_size = 100;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 1;               % p-vector norm for class prediction
train_trials = [1 2];

bin_size = 10;
PSTH_smooth = 1;            % PSTH smoothing factor
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix

syspath = '/Users/Xander/Documents/MATLAB/Paper_final';
codepath = '/Users/Xander/Documents/MATLAB/Paper_final/Code_final';
var_ext = '_data_master';
filepath = 'Data/Master_files';

% Process only cell culture or cell culture and odors
if strcmp(include_odors,'n')
    order = [6 1 2 3 5];
    writepath = 'Figures/NoOdors';
elseif strcmp(include_odors,'y')
    order = [6 1 2 3 5 4 7];
    writepath = 'Figures/Odors';
end

% Specify color scheme
Colors = [246,145,30 ; 255,0,255 ; 0,255,0 ; 150,150,150 ; 0,0,0 ; 255,0,0 ; 117,76,36]./255; 

% Adjust color order in accordance with stimulus order 
Colors = Colors(order,:);

% RMS-transformed or spike-sorted data?
if strcmp(file_type,'r')
    filename = 'R500_S500';
    read_filename = [filename '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    write_filename = [filename '_b' num2str(new_bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
elseif strcmp(file_type,'s') 
    filename = 'SS';
    read_filename = [filename '_b' num2str(bin_size)];
    write_filename = [filename '_b' num2str(new_bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
end

% Perform daywise comparisons
if strcmp(daywise_comp,'y')
    select_day = ['Day' num2str(input('Select day (1, 2, 3, 4): ','s'))];
    read_filename = [read_filename '_' select_day];
    write_filename = [write_filename '_' select_day];
%     readfile = [read_filename '_' select_day];
end

load([syspath '/' filepath '/' read_filename '.mat']);
files = who([syspath '/' filepath '/' read_filename '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});

    % Compile all stimuli data into master tensor
    if strcmp(file_type,'s')
        data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
            - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
        data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
    elseif strcmp(file_type,'r')
        data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
    end
end
% data (RMS) = positions, binned samples, trials, stimuli
% data (SS) = neurons, binned samples, trials, stimuli 


% Generate figures
addpath(codepath)
set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
titles = [num2str(time(1)) ' to ' num2str(time(2)) ' seconds after stimulus onset'];

PCA_master_fig = PCA(data, Colors, stimuli, PCA_smooth, ['PCA (' file_type '): ' titles], time, new_bin_size, bins_per_sec, PCA_line_mult);        
LDA_master_fig = LDA(data, Colors, stimuli, ['LDA (' file_type '): ' titles]);
figs = {PCA_master_fig,LDA_master_fig};
figs_titles = {'PCA','LDA'};

[bin_preds, trial_preds] = TrainTestLocust(data,normtype,train_trials);
Conf_fig_bin = Confusion(bin_preds, stimuli, ['Confusion (' file_type '): ' titles],conf_res,1);
Conf_fig_trial = Confusion(trial_preds, stimuli, ['Confusion (' file_type '): ' titles],conf_res,1);
norm_figs = {Conf_fig_bin,Conf_fig_trial};
norm_figs_titles = {'Bin_Conf','Trial_Conf'};

[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, ['CV Confusion (' file_type '): ' titles],conf_res,1);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, ['CV Confusion (' file_type '): ' titles],conf_res,1);
norm_figs = cat(2,norm_figs,{Conf_fig_cv_bin,Conf_fig_cv_trial});
norm_figs_titles = cat(2,norm_figs_titles, {'CV_Bin_Conf','CV_Trial_Conf'});

% Save figures if specified
if strcmp(save_figs,'y')
%     write_filename = [filename '_b' num2str(new_bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    fig_format = {'fig', 'png'};
     
    for cycle_fig_formats = 1:numel(fig_format)
        for cycle_figs = 1:numel(figs)
            figure(figs{cycle_figs});
            
            % Make directory and save norm independent figures
            if ~exist([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename],'dir')
                mkdir([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename]);
                saveas(figs{cycle_figs},[syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}])
            elseif exist([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename],'dir')... 
                    && ~exist([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}],'file')
                saveas(figs{cycle_figs},[syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}])
            else
                fprintf([figs_titles{cycle_figs} '_' num2str(time(1)) 'to' num2str(time(2)) ' not saved. File with same name already exists\n'])
            end
        end
            
        for cycle_norm_figs = 1:numel(norm_figs)
            figure(norm_figs{cycle_norm_figs});
            
            % Make directory and save norm dependent figures
            if ~exist([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype)],'dir')
                mkdir([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename  '/Norm' num2str(normtype) ]);
                saveas(norm_figs{cycle_norm_figs},[syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename  '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}])
            elseif exist([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype) ],'dir')... 
                    && ~exist([syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}],'file')
                saveas(norm_figs{cycle_norm_figs},[syspath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}])
            else
                fprintf([norm_figs_titles{cycle_norm_figs} '_' num2str(time(1)) 'to' num2str(time(2)) ' not saved. File with same name already exists\n'])
            end 
        end
    end
    close all;
    fprintf(['Completed ' write_filename ' Norm ' num2str(normtype) '\n'])
else
    close_figs = 'Close all figures (y / n)? ';
    close_figs = input(close_figs,'s');
    if close_figs == 'y'
        close all;
    end
end
toc;


% PSTH comparison for each neuron
% angular distance inclusion, see Deb's code
