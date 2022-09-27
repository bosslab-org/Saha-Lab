clear; clc; close all; tic;

file_type = 'r';        % s (Spike Sorted) or r (RMS)
save_figs = 'n';        % y or n
exp_type = '1%';        % 1% or LC (Lung Cancer)

time = [0 4];               % time relative to stimulus onset
new_bin_size = 200;          % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction

bin_size = 10;
PSTH_smooth = 1;            % PSTH smoothing factor
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;        % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix

if strcmp(exp_type,'LC') 
    readpath = '/Users/Xander/Documents/MATLAB/Honeybee/LC';
    
    % Specify color scheme
    Colors = [242,147,147 ; 189,67,67 ; 237,70,47 ; 165,245,167 ; 57,123,39 ;...
        69,148,39 ; 79,201,251 ; 49,74,251 ; 147,147,147 ; 38,38,38]./255;
elseif strcmp(exp_type,'1%')
    readpath = '/Users/Xander/Documents/MATLAB/Honeybee/1%';
    
    % Specify color scheme
    Colors = [89,190,73 ; 131,195,232 ; 255,60,45 ; 44,46,53 ; 200,206,205 ;...
        255,210,214 ; 95,209,213 ; 213,226,127 ; 150,58,158 ; 209,44,28]./255;
end
filepath = 'Master_files';
writepath = 'Figures';
    
% % Adjust color order in accordance with stimulus order 
% Colors = Colors(order,:);

% RMS-transformed or spike-sorted data?
if strcmp(file_type,'r')
    filename = 'R500_S500';
    read_filename = [filename '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    write_filename = [filename '_b' num2str(new_bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    var_ext = '_RMS';
    train_ext = '_1_RMS';
elseif strcmp(file_type,'s') 
    filename = 'SS';
    read_filename = [filename '_b' num2str(bin_size)];
    write_filename = [filename '_b' num2str(new_bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    var_ext = '_spikes';
    train_ext = '_1_spikes';
end

load([readpath '/' filepath '/' read_filename '.mat']);
files = who([readpath '/' filepath '/' read_filename '.mat'], ['*' var_ext]);
% files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});

    % Compile all stimuli data into master tensor
    if strcmp(file_type,'s')
        data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
            - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
        data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
    elseif strcmp(file_type,'r')
        data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
    end
end
% data.shape (RMS) = positions, binned samples, trials, stimuli
% data.shape (SS) = neurons, binned samples, trials, stimuli 

% Generate figures
set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
titles = [num2str(time(1)) ' to ' num2str(time(2)) ' seconds after stimulus onset'];

PCA_master_fig = PCA_Honeybee(data, Colors, stimuli, PCA_smooth, ['PCA (' file_type '): ' titles], time, new_bin_size, bins_per_sec, PCA_line_mult);        
LDA_master_fig = LDA_Honeybee(data, Colors, stimuli, ['LDA (' file_type '): ' titles]);
figs = {PCA_master_fig,LDA_master_fig};
figs_titles = {'PCA','LDA'};

[loto_bin_preds, loto_trial_preds] = LotoBee(data,normtype);
Conf_loto_bin = Confusion_Honeybee(loto_bin_preds, stimuli, stimuli, ['LOTO Confusion (' file_type '): ' titles],conf_res,1);
Conf_loto_trial = Confusion_Honeybee(loto_trial_preds, stimuli, stimuli, ['LOTO Confusion (' file_type '): ' titles],conf_res,1);
norm_figs = {Conf_loto_bin,Conf_loto_trial};
norm_figs_titles = {'LOTO_Bin_Conf','LOTO_Trial_Conf'};

if strcmp(exp_type,'LC') 
    [tt_bin_preds, tt_trial_preds, train_stimuli, test_stimuli] = TtBee(data,normtype,files,train_ext);
    Conf_tt_bin = Confusion_Honeybee(tt_bin_preds, train_stimuli, test_stimuli, ['TT Confusion (' file_type '): ' titles],conf_res,1);
    Conf_tt_trial = Confusion_Honeybee(tt_trial_preds, train_stimuli, test_stimuli, ['TT Confusion (' file_type '): ' titles],conf_res,1);
    norm_figs = cat(2,norm_figs,{Conf_tt_bin,Conf_tt_trial});
    norm_figs_titles = cat(2,norm_figs_titles, {'TT_Bin_Conf','TT_Trial_Conf'});
end

% Save figures if specified
if strcmp(save_figs,'y')
    fig_format = {'fig', 'png'};
    for cycle_fig_formats = 1:numel(fig_format)
        for cycle_figs = 1:numel(figs)
            figure(figs{cycle_figs});
            
            % Make directory and save norm independent figures
            if ~exist([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename],'dir')
                mkdir([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename]);
                saveas(figs{cycle_figs},[readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}])
            elseif exist([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename],'dir')... 
                    && ~exist([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}],'file')
                saveas(figs{cycle_figs},[readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}])
            else
                fprintf([figs_titles{cycle_figs} '_' num2str(time(1)) 'to' num2str(time(2)) ' not saved. File with same name already exists\n'])
            end
        end
            
        for cycle_norm_figs = 1:numel(norm_figs)
            figure(norm_figs{cycle_norm_figs});
            
            % Make directory and save norm dependent figures
            if ~exist([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype)],'dir')
                mkdir([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename  '/Norm' num2str(normtype) ]);
                saveas(norm_figs{cycle_norm_figs},[readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename  '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}])
            elseif exist([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype) ],'dir')... 
                    && ~exist([readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}],'file')
                saveas(norm_figs{cycle_norm_figs},[readpath '/' writepath '/' upper(fig_format{cycle_fig_formats}) '/' write_filename '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}])
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
