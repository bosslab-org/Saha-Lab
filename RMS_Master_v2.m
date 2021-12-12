%% Master File for Processing RMS Data
clear; clc; close all;
stimulus = 'Is/Are the file(s) of odorants ("o"), cell culture ("c") or honeybee ("h")?  ';
stimulus = input(stimulus,'s');
input_files = 'Process multiple ("m") or individual ("s") file(s)? ';
input_files = input(input_files,'s');

filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables below this line
time = [-2 4];          % time relative to stimulus onset
new_bin_size = 100;

if input_files == 'm'
    filepath2 = 'Master_files';
    master_file = 'Odor_RMS_Data_Filt';
    load([filepath '/' filepath2 '/' master_file]);
    PCA_smooth = 2;
    PCA_line_mult = 500;      % line from origin every x msecs, MUST be set to multiple of new bin_size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edit variables above this line
    time = [stimulus_on+time(1) stimulus_on+time(2)];
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

elseif input_files == 's'
    Date = 'Specify date to be processed (MM_DD_YYYY): ';
    Date = input(Date,'s');
    Odorant = 'Specify odorant to be processed: ';
    Odorant = input(Odorant,'s');
    Position = 'Specify position to be processed: ';
    Position = input(Position,'s');
    
    if stimulus == 'o'
        filepath2 = 'MAT_files/Odors';
    elseif stimulus == 'c'
        filepath2 = 'Cell_Culture';
    elseif stimulus == 'h'
        filepath2 = 'Honeybee';
    end

    odor_file = dir([filepath '/' filepath2 '/' Date '/Position_' Position '/' Odorant '.mat']);
    set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
           
    load([odor_file.folder '/' odor_file.name]);
    
    time = [stimulus_onset+time(1) stimulus_onset+time(2)];
    exp_pars = struct('Date',Date,'Odorant',Odorant,'Position',Position, 'Time',time,'stim_on',stimulus_onset,'stim_off',stimulus_offset,'sample_rate',sample_rate,'new_bin_size',new_bin_size);
    odor_data = eval([Odorant '_data_filt']);
    for cycle_tets = 1:size(odor_data,1)
        tet_data = odor_data(cycle_tets,:,:,time(1)*sample_rate+1:time(2)*sample_rate);
        for cycle_chans = 1:size(tet_data,2)            
%             RMS_PSTH_file(exp_pars,file,cycle_tetrodes,cycle_channels,time)s
%             RMS_Raster_file(exp_pars,file,cycle_tetrodes,cycle_channels,time)
            RMS_PSVT_file(exp_pars,tet_data(:,cycle_chans,:,:),cycle_tets,cycle_chans,time)
        end
    end
else
    fprintf(1,'Please enter a valid option and run again \n');
end