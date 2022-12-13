% Harnessing Insect Olfactory Circuits for Detecting and Discriminating Human Cancers

%% Figure 1 - Individual projection neurons respond differentially to the oral cancer vs. control VOCs
clear; clc; close all;
basepath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data';
rmspath = 'MAT_files';

% Example position 1
date1 = '08_05_2021';
position1 = 1;
tetrode1 = 1;
channel1 = 3;
trial1 = [1 3 1 2 3 4 3];

% Example position 2
date2 = '07_22_2021';
position2 = 3;
tetrode2 = 1;
channel2 = 4;
trial2 = [5 5 5 5 5 5 5];

colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5 4 7];

PSVT([basepath '/' rmspath],date1,position1,tetrode1,channel1,trial1,colors,order)
% PSVT([basepath '/' rmspath],date2,position2,tetrode2,channel2,trial2,colors,order)


masterpath = 'Master_files/Cell_Culture_SS_Master';

pvalue = 0.05;
base_stim_idx = 1;
test_stim_idx = 2:7;
bin_size = 50;

Neuron_ANOVA([basepath '/' masterpath],pvalue,base_stim_idx,test_stim_idx,bin_size)

%% Figure 2 - Cancer vs. non-cancer VOCs are distinguished by spatiotemporal PN responses
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files/SS_b10';

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

%% Figure 3 - Neural response-based classification of oral cancer is robust at different time-matched points of growth of cancer vs. non-cancer cells
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files/SS_b10';

% Specify Day 1, 2, 3 or 4
Day = 2;

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

%% Figure 4 - Rapid classification of oral cancer VOC profiles using neural voltage responses
clear; clc; close all;
basepath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files';

% Specify time window of interest relative to stimulus onset
time = [0 4];  %[0.75 1] , [2 2.25] , [2.25 2.5]

new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 
order = [6 1 2 3 5];

filename = ['R500_S500_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
% var_ext = '_RMS';
var_ext = '_data_master';


load([basepath '/' filename '.mat']);
files = who([basepath '/' filename '.mat'], ['*' var_ext]);
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


load([basepath '/R500_S500_b' num2str(bin_size) '_0to4.mat']);
files = who([basepath '/R500_S500_b' num2str(bin_size) '_0to4.mat'], ['*' var_ext]);
files = files(order);

for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
end

pairwise_distance_master_fig = pairwise_distance(data);


%% Figure 5- Principles of neuronal response-based noninvasive cancer detection
% Made in Adobe Illustrator

%% Figure S1 - Electrophysiological recording setup
% Made in Adobe Illustrator

%% Figure S2 - PN population response-based classification of the entire stimulus panel
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files/SS_b10';

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

%% Figure S3 - Cancer vs. non-cancer classification using non-overlapping train-test PN response datasets
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files/SS_b10';

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
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

train_set = 1:2;
[cv_bin_preds, cv_trial_preds] = Train_test(data,normtype,train_set);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

%% Figure S4 - Spatiotemporal PN response-based VOC mixture classification is robust
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files/SS_b10';

time = [0 4];               % time relative to stimulus onset
new_bin_size = [50 100 200];         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5];

var_ext = '_data_master';

load([readpath '.mat']);
files = who([readpath '.mat'], ['*' var_ext]);
files = files(order);

for cycle_bin_size = 1:numel(new_bin_size)
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
    Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
    Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);
    clear data
end
%% Figure S5 - Representative images of cell cultures over days
% Made in Adobe Illustrator

%% Figure S6 - Depiction of the cell counting procedure
% Made in Adobe Illustrator

%% Figure S7 - Cell growth curves of all four cell lines
% Made in Adobe Illustrator

%% Figure S8 - Neural response-based classification of oral cancer remains robust with non-overlapping training-testing datasets
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files/SS_b10';

% Specify Day 1, 2, 3 or 4
Day = 2;

time = [0 4];               % time relative to stimulus onset
new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5];

var_ext = '_data_master';

load([readpath '_Day' num2str(Day)]);
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

train_set = 1:2;
[cv_bin_preds, cv_trial_preds] = Train_test(data,normtype,train_set);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

%% Figure S9 - Root-mean squared (R.M.S) transformation largely preserved stimulus-specific spiking dynamics.
clear; clc; close all;
readpath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/MAT_files';

% Example position 1
date = '08_05_2021';
position = 1;
tetrode = 1;
channel = 3;
trial = [1 3 1 2 3 4 3];

PSTH_smooth = 1;            % PSTH smoothing factor
new_bin_size = 10000;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
bin_size = 10;
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0 ; 150,150,150 ; 117,76,36]./255; 
order = [6 1 2 3 5 4 7];

bin_mult = new_bin_size/bin_size;

PSVT(readpath,date,position,tetrode,channel,trial,colors,order)
rms_PSVT(readpath,date,position,tetrode,channel,trial,bin_mult,colors,order)

%% Figure S10 - Using R.M.S transformed population PN voltage responses to classify the stimulus panel
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

%% Figure S11 - Cancer VOC classification during transient vs. steady state response periods of PN response
clear; clc; close all;
basepath = '/Users/alexanderfarnum/Documents/MATLAB/Biosensors_code/Data/Master_files';

% Specify time window of interest relative to stimulus onset
time = [0 4];  %[-1 0.5] , [0.5 2] , [2 3.5]

new_bin_size = 50;         % new bin size in msecs: mod(size(data_temp,3),new_bin_size) must equal 0
normtype = 2;               % p-vector norm for class prediction
bin_size = 10;
PCA_smooth = 3;             % PCA smoothing factor
PCA_line_mult = 200;         % line from origin every x msecs, MUST be set to multiple of new bin_size
conf_res = 20;              % resolution (step size) of confusion matrix
colors = [255,0,0 ; 246,145,30 ; 255,0,255 ; 0,255,0 ; 0,0,0]./255; 
order = [6 1 2 3 5];

filename = ['R500_S500_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];
var_ext = '_RMS';

load([basepath '/' filename '.mat']);
files = who([basepath '/' filename '.mat'], ['*' var_ext]);
files = files(order);

bin_mult = new_bin_size/bin_size;
bins_per_sec = 1000/new_bin_size;
for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
end

[cv_bin_preds, cv_trial_preds] = LeaveTrialOutLocust(data,normtype);
Conf_fig_cv_bin = Confusion(cv_bin_preds, stimuli, 'Confusion',conf_res,0);
Conf_fig_cv_trial = Confusion(cv_trial_preds, stimuli, 'Confusion',conf_res,0);

clear data

load([basepath '/R500_S500_b' num2str(bin_size) '_0to4.mat']);
var_ext = '_data_master';
files = who([basepath '/R500_S500_b' num2str(bin_size) '_0to4.mat'], ['*' var_ext]);
files = files(order);

for cycle_stimuli = 1:numel(files)
    stimuli{cycle_stimuli} = files{cycle_stimuli}(1:end-numel(var_ext));
    data_temp = eval(files{cycle_stimuli});
    data(:,:,:,cycle_stimuli) = permute(mean(reshape(data_temp,size(data_temp,1),size(data_temp,2),bin_mult,size(data_temp,3)/bin_mult),3),[1 4 2 3]);
end

delta_rms_PSVT(data,colors,stimuli)
