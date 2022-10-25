% Harnessing...

%% Figure 1 - Peri-Stimulus Voltage Traces (PSVTs)
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data';
rmspath = 'RMS/Day_files';

% Example position 1
date1 = '08_05_2021';
position1 = 1;
tetrode1 = 1;
channel1 = 3;
trial1 = [1 1 1 1 1 1 1];

% Example position 2
date2 = '07_22_2021';
position2 = 3;
tetrode2 = 1;
channel2 = 4;
trial2 = [1 1 1 1 1 1 1];

PSVT([readpath '/' rmspath],date1,position1,tetrode1,channel1,trial1)
PSVT([readpath '/' rmspath],date2,position2,tetrode2,channel2,trial2)

masterpath = 'Master_files/Cell_Culture_SS_Master';
var_ext = '_spikes_bin';

pvalue = 0.05;
base_stim_idx = 1;
test_stim_idx = [2:7];
bin_size = 50;

Neuron_ANOVA([readpath '/' masterpath],pvalue,base_stim_idx,test_stim_idx,bin_size)

%% Figure 2 - Spike Sorted without Odors
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files/SS_b10';

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 

order = [6 1 2 3 5];
var_ext = '_data_master';

load([readpath '.mat']);
files = who([readpath '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
        - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
    data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
end

PCA_master_fig = PCA(data, colors, stimuli, PCA_smooth, 'PCA', time, new_bin_size, bins_per_sec, PCA_line_mult);        
LDA_master_fig = LDA(data, colors, stimuli, 'LDA');
[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

%% Figure 3 - Spike Sorted without Odors: Daywise Comparison
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files/SS_b10';

% Specify Day 1, 2, 3 or 4
Day = 1;

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 

order = [6 1 2 3 5];
var_ext = '_data_master';

load([readpath '_Day' num2str(Day)]);
files = who(readpath, ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
        - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
    data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
end

PCA_master_fig = PCA(data, colors, stimuli, PCA_smooth, 'PCA', time, new_bin_size, bins_per_sec, PCA_line_mult);        
LDA_master_fig = LDA(data, colors, stimuli, 'LDA');
[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

%% Figure 4 - RMS time windows without Odors
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files';

% Specify time window of interest realtive to stimulus onset
time = [0 4];         % [0.5 0.75] , [0.75 1] , [2 2.25] , [2.25 2.5]

new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 
order = [6 1 2 3 5];

filename = ['R500_S500_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
var_ext = '_data_master';

load([readpath '/' filename '.mat']);
files = who([readpath '/' filename '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
end

PCA_master_fig = PCA(data, colors, stimuli, PCA_smooth, 'PCA', time, new_bin_size, bins_per_sec, PCA_line_mult);        
[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

% TO ADD: DELTA RMS PSVT

%% Figure S1 - Spike Sorted with Odors
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files/SS_b10';

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5 4 7];

var_ext = '_data_master';

load([readpath '.mat']);
files = who([readpath '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
        - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
    data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
end

PCA_master_fig = PCA(data, colors, stimuli, PCA_smooth, 'PCA', time, new_bin_size, bins_per_sec, PCA_line_mult);        
LDA_master_fig = LDA(data, colors, stimuli, 'LDA');
[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

%% Figure S5 - PSVT and RMS transformed data
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/RMS/Day_files';

% Example position 1
date = '08_05_2021';
position = 1;
tetrode = 1;
channel = 3;
trial = [1 1 1 1 1 1 1];

PSTH_smooth = 1;            % PSTH smoothing factor
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
bin_size = 10;
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5 4 7];

PSVT(readpath,date,position,tetrode,channel,trial)
rms_PSVT(readpath,date,position,tetrode,channel,trial)


%% Figure S6 - RMS with Odors
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files';

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5 4 7];

filename = ['R500_S500_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
var_ext = '_data_master';

load([readpath '/' filename '.mat']);
files = who([readpath '/' filename '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
end

PCA_master_fig = PCA(data, colors, stimuli, PCA_smooth, 'PCA', time, new_bin_size, bins_per_sec, PCA_line_mult);        
LDA_master_fig = LDA(data, colors, stimuli, 'LDA');
[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

%% Figure S7 - RMS baseline, transient and steady state
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files';

% Specify time window of interest realtive to stimulus onset
time = [0 4];         % [-1 0.5] , [0.5 2] , [2 2.25] , [2 3.5]

new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 
order = [6 1 2 3 5];

filename = ['R500_S500_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
var_ext = '_data_master';

load([readpath '/' filename '.mat']);
files = who([readpath '/' filename '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
end

PCA_master_fig = PCA(data, colors, stimuli, PCA_smooth, 'PCA', time, new_bin_size, bins_per_sec, PCA_line_mult);        
[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

% TO ADD: DELTA RMS PSVT

%% Figure S8 - Spike Sorted Leave-trial-out vs. Train-Test set confusions
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files/SS_b10';

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 
train_trials = [1 2];
order = [6 1 2 3 5];
var_ext = '_data_master';

load([readpath '.mat']);
files = who([readpath '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
        - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
    data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
end

[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Leave-Trial-Out Confusion Bin' ,conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Leave-Trial-Out Confusion Trial' ,conf_res,0);

[bin_preds, trial_preds] = TrainTestLocust(data,normtype,train_trials);
Conf_fig_bin = Confusion(bin_preds, stimuli, 'Train-Test Confusion Bin' ,conf_res,0);
Conf_fig_trial = Confusion(trial_preds, stimuli, 'Train-Test Confusion Trial' ,conf_res,0);

%% Figure S9 - Spike Sorted Varying bin sizes
clear; clc; close all;
readpath = '/Users/Xander/Documents/MATLAB/Paper_final/Data/Master_files/SS_b10';

time = [0 4];               % time relative to stimulus onset
new_bin_size = [50 100 200];         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 
train_trials = [1 2];
order = [6 1 2 3 5];
var_ext = '_data_master';

load([readpath '.mat']);
files = who([readpath '.mat'], ['*' var_ext]);
files = files(order);

for cycle_bin_size = 1:numel(new_bin_size)
    clear data
    bin_mult = new_bin_size(cycle_bin_size)/bin_size;
    bins_per_sec = 1000/new_bin_size(cycle_bin_size);
    for cycle_stimuli = 1:numel(files)
        stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
        data_temp = eval(files{cycle_stimuli});
        data_norm = data_temp(:,:,(time(1)+stim_on)*bins_per_sec*bin_mult+1:(time(2)+stim_on)*bins_per_sec*bin_mult)...
            - mean(data_temp(:,:,(stim_on-2)*bins_per_sec*bin_mult+1:stim_on*bins_per_sec*bin_mult),[2,3]);
        data(:,:,:,cycle_stimuli) = permute(sum(reshape(data_norm,size(data_temp,1),size(data_temp,2),bin_mult,numel((time(1)+stim_on)*bins_per_sec+1:(time(2)+stim_on)*bins_per_sec)),3),[1 4 2 3]);
    end

    [cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
    Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, ['Confusion Bin: ' num2str(new_bin_size(cycle_bin_size)) ' msec binsize'] ,conf_res,0);
    Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, ['Confusion Trial: ' num2str(new_bin_size(cycle_bin_size)) ' msec binsize'] ,conf_res,0);
end