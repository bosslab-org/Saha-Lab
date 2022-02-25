clear; clc; close all; tic;

rms_window = 2000;  %bin size in number of samples (msecs = (sample_rate/(rms_window*1000))
smooth_window = 2000;
bin_size = 10; %bin size in msecs
filepath_temp = 'Cell_Culture'; % Cell_Culture Odors Honeybee

filepath = ['/Users/Xander/Documents/MATLAB/Neural_Recordings/RMS/' filepath_temp];
readpath = 'Day_files';
writepath = ['R' num2str(rms_window) '_S' num2str(smooth_window) '_b' num2str(bin_size)];

dates = cat(1,dir([filepath '/' readpath '/*2021']),dir([filepath '/' readpath '/*2022']));
for cycle_dates = 1:numel(dates)
    date = dates(cycle_dates).name;
    if ~exist([filepath '/' writepath '/' date],'dir')
        mkdir([filepath '/' writepath '/' date]);
    end
    writepath_temp = [filepath '/' writepath '/' date];
    positions = dir([filepath '/' readpath '/' date]);
    positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
    for cycle_positions = 1:numel(positions)
        experiment = strings;
        data = [];
        position = positions(cycle_positions).name;
        fprintf(['Processing ' date ' ' position '...']);
        odors = dir([filepath '/' readpath '/' date '/' position]);
        odors = odors(~ismember({odors.name},{'.','..','.DS_Store',[position '.mat']}));
        for cycle_odors = 1:numel(odors)
            odor{cycle_odors} = odors(cycle_odors).name(1:end-4);
            load([filepath '/' readpath '/' date '/' position '/' odors(cycle_odors).name]);
            data_temp = eval([odor{cycle_odors} '_data_filt']);
            for cycle_tetrodes = 1:size(data_temp,1)
                if cycle_odors == 1
                    experiment(cycle_tetrodes,:) = [date '_' position '_Tetrode_' num2str(cycle_tetrodes)];
                end
                for cycle_trials = 1:size(data_temp,3)
                    data(cycle_tetrodes,cycle_odors,cycle_trials,:) = mean(smoothdata(sqrt(movmean(data_temp(cycle_tetrodes,:,cycle_trials,(stim_on-2)*sample_rate+1:(stim_off+4)*sample_rate).^ 2, rms_window, 4)),4,'movmean',smooth_window),2);             
                end
            end
        end        
        rms_data = permute(rms(reshape(data,size(data,1),size(data,2),size(data,3),(sample_rate*bin_size)/1000,[]),4),[1 2 3 5 4]);
        mean_data = permute(mean(reshape(data,size(data,1),size(data,2),size(data,3),(sample_rate*bin_size)/1000,[]),4),[1 2 3 5 4]);
        
        fprintf(['Saving ' date ' ' position ' data... \n']);
        for cycle_odors = 1:numel(odors)
            mean_var_name = [odor{cycle_odors} '_mean_data'];            
            eval([mean_var_name '= permute(mean_data(:,cycle_odors,:,:),[1 3 4 2]);'])
            if cycle_odors == 1
                stim_on = 2; stim_off = 6; total_time = 10; sample_rate = 20000;
                save([writepath_temp '/' position], [odor{cycle_odors} '_mean_data'],'experiment','stim_on','stim_off','total_time','sample_rate','bin_size');
            else
                save([writepath_temp '/' position], [odor{cycle_odors} '_mean_data'],'-append');
            end
        end
    end
    fprintf([num2str(cycle_dates/numel(dates)*100) '%% complete \n\n']);
end
fprintf(['Finished ' writepath '\n']);
toc;