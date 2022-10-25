clear; clc; close all; tic;

exp_type = 'l1';            % l1: Locust 1% biomarkers  
                            % lc: Locust cell culture
                            % h1: Honeybee 1% biomarkers  
                            % hb: Honeybee breath mixture

file_type = 'r';            % s (Spike Sorted) or r (RMS)
save_figs = 'n';            % y or n

time = [0 4];               % time relative to stimulus onset
bin_size = 10;              % bin size of saved data (msecs)
new_bin_size = 50;          % new bin size (msecs) -- must be multiple of bin_size
normtype = 2;               % Lp-vector norm for class prediction
test_lambda = 0;           % cycle lambda values and generate LDA projection plots

smooth_window = 500;        % RMS ONLY: smooth window of saved data (number of samples)
rms_window = 500;           % RMS ONLY: rms window of saved data (number of samples)

inc_odors = 'y';            % Locust cell culture ONLY: include 1% odorants?
select_day = 'n';             % Locust cell culture ONLY: process individual day? 1, 2, 3, 4, n


% Plotting parameters
PSTH_smooth = 1;            % PSTH smoothing factor
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;        % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix


basepath = '/Users/alexanderfarnum/Documents/MATLAB';
filepath = 'Master_files';

if strcmp(file_type,'r')
    filename = ['R' num2str(rms_window) '_S' num2str(rms_window) '_' exp_type '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    write_filename = ['R' num2str(rms_window) '_S' num2str(rms_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    var_ext = '_RMS';
%     var_ext = '_data_master';
elseif strcmp(file_type,'s') 
    filename = ['SS_' exp_type '_b' num2str(bin_size)];
    write_filename = ['SS_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
    var_ext = '_SS';
end

if strcmp(exp_type,'l1') || strcmp(exp_type,'h1')
    Colors = [89,190,73 ; 131,195,232 ; 255,60,45 ; 44,46,53 ; 200,206,205 ;...
        255,210,214 ; 95,209,213 ; 213,226,127 ; 150,58,158 ; 209,44,28]./255;
    order = 1:10;
elseif strcmp(exp_type,'lc')
    Colors = [255,0,0 ; 0,255,0 ; 150,150,150 ; 255,255,0 ; 0,0,0 ; 0,0,255 ; 0,255,255]./255; 
    if strcmp(inc_odors,'y')
        order = [6 1 2 3 5 4 7];
    elseif strcmp(inc_odors,'n')
        order = [6 1 2 3 5];
    end

    switch select_day
        case '1', filename = cat(2, filename, '_Day1');
        case '2', filename = cat(2, filename, '_Day2');
        case '3', filename = cat(2, filename, '_Day3');
        case '4', filename = cat(2, filename, '_Day4');
    end
elseif strcmp(exp_type,'hb')
    Colors = [180,0,0 ; 0,220,0 ; 102,178,255 ; 255,10,10 ; 192,192,192 ;...
        255,200,200 ; 0,225,225 ; 204,255,103 ; 153,0,153 ; 50,50,50]./255;
    order = 1:10;
end

load([basepath '/' filepath '/' filename '.mat']);
files = who([basepath '/' filepath '/' filename '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    if strcmp(file_type,'s')
        % Spike Sorted: subtract 2 second baseline and compile master tensor
        data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
            - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
        data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
    elseif strcmp(file_type,'r')
        % RMS: compile master tensor (baseline already subtracted in preprocessing)
        data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
    end
end
% size(data) = neurons, binned samples, trials, stimuli 


% Generate figures
set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
titles = [num2str(time(1)) ' to ' num2str(time(2)) ' seconds after stimulus onset'];

PCA_master_fig = Pca(data, Colors, stimuli, PCA_smooth, ['PCA (' file_type '): ' titles], time, new_bin_size, bins_per_sec, PCA_line_mult);        

if any(test_lambda)
    LDA_lambdas_fig = Lda_lambdas(data, Colors, stimuli, ['LDA (' file_type '): ' titles], test_lambda);
    lambda = 'Select lambda value: ';
    lambda = input(lambda,'s');
    close(LDA_lambdas_fig)
else
    lambda = 0;
end

LDA_master_fig = Lda(data, Colors, stimuli, ['LDA (' file_type '): ' titles], lambda);

figs = {PCA_master_fig,LDA_master_fig};
figs_titles = {'PCA','LDA'};



[loto_bin_preds, loto_trial_preds] = Loto(data,normtype);
Conf_loto_bin = Confusion(loto_bin_preds, stimuli, stimuli, ['LOTO Confusion (' file_type '): ' titles],conf_res,1);
Conf_loto_trial = Confusion(loto_trial_preds, stimuli, stimuli, ['LOTO Confusion (' file_type '): ' titles],conf_res,1);
norm_figs = {Conf_loto_bin,Conf_loto_trial};
norm_figs_titles = {'LOTO_Bin_Conf','LOTO_Trial_Conf'};


if strcmp(exp_type,'hb')
    [tt_bin_preds, tt_trial_preds] = TrainTest(data,normtype);
    Conf_tt_bin = Confusion(tt_bin_preds, train_stimuli, test_stimuli, ['TrainTest Confusion (' file_type '): ' titles],conf_res,1);
    Conf_tt_trial = Confusion(tt_trial_preds, train_stimuli, test_stimuli, ['TrainTest Confusion (' file_type '): ' titles],conf_res,1);
    norm_figs = cat(2, norm_figs, {Conf_tt_bin,Conf_tt_trial});
    norm_figs_titles = cat(2, norm_figs_titles, {'TrainTest_Bin_Conf','TrainTest_Trial_Conf'});
end


% Save figures
if strcmp(save_figs,'y')
    if strcmp(exp_type,'l1')
        writepath = 'Master_figs/Locust/1%';
    elseif strcmp(exp_type,'lc')
        writepath = 'Master_figs/Locust/CellCulture'; 
    elseif strcmp(exp_type,'h1')
        writepath = 'Master_figs/Honeybee/1%'; 
    elseif strcmp(exp_type,'hb')
        writepath = 'Master_figs/Honeybee/Breath';
    end
    fig_format = {'fig', 'png'};
    for cycle_fig_formats = 1:numel(fig_format)
        for cycle_figs = 1:numel(figs)
            figure(figs{cycle_figs});
            % Make directory and save norm independent figures
            if ~exist([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats})],'dir')
                mkdir([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats})]);
                saveas(figs{cycle_figs},[basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats}) '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}])
            elseif exist([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats})],'dir')... 
                    && ~exist([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats}) '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}],'file')
                saveas(figs{cycle_figs},[basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats}) '/' figs_titles{cycle_figs} '.' fig_format{cycle_fig_formats}])
            else
                fprintf([figs_titles{cycle_figs} '_' num2str(time(1)) 'to' num2str(time(2)) ' not saved. File with same name already exists\n'])
            end
        end
            
        for cycle_norm_figs = 1:numel(norm_figs)
            figure(norm_figs{cycle_norm_figs});
            % Make directory and save norm dependent figures
            if ~exist([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats}) '/Norm' num2str(normtype)],'dir')
                mkdir([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats})  '/Norm' num2str(normtype) ]);
                saveas(norm_figs{cycle_norm_figs},[basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats})  '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}])
            elseif ~exist([basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats}) '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}],'file')
                saveas(norm_figs{cycle_norm_figs},[basepath '/' writepath '/' write_filename '/' upper(fig_format{cycle_fig_formats}) '/Norm' num2str(normtype) '/' norm_figs_titles{cycle_norm_figs} '.' fig_format{cycle_fig_formats}])
            else
                fprintf([norm_figs_titles{cycle_norm_figs} '_' num2str(time(1)) 'to' num2str(time(2)) ' not saved. File with same name already exists\n'])
            end 
        end
    end
    close all;
    fprintf(['Finished ' write_filename ' Norm ' num2str(normtype) '\n'])
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
