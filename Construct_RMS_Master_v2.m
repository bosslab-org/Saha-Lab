clear; clc; close all;

% select_day = 'Day1'; % comment out for aggregate master, comment in for daywise master
filename = 'R2000_S2000_b10';
exp_type = 'Cell_Culture'; % Cell_Culture Odors Honeybee

Day1 = {'05_25_2021';'06_24_2021';'06_29_2021';'07_20_2021';'08_03_2021';'08_10_2021'};
Day2 = {'05_26_2021';'06_30_2021';'07_14_2021';'08_04_2021';'08_11_2021'};
Day3 = {'05_27_2021';'06_26_2021';'07_01_2021';'07_15_2021';'07_22_2021';'08_05_2021';'08_12_2021'};
Day4 = {'05_28_2021';'07_02_2021';'07_16_2021';'07_23_2021';'08_06_2021'};

filepath = '/Users/Xander/Documents/MATLAB/Neural_Recordings';
odorpath = [filepath '/RMS/' exp_type '/Day_files'];
readpath = [filepath '/RMS/' exp_type '/' filename];
 
dates = cat(1,dir([readpath '/*2021']),dir([readpath '/*2022'])); %pull all files with specified year suffix

if exist('select_day','var')
    writepath = [filepath '/Master_files/' filename '_' exp_type '_' select_day '.mat'];
    dates = dates(ismember({dates.name},eval(select_day)));
else
    writepath = [filepath '/Master_files/' filename '_' exp_type '.mat'];
end
for cycle_dates = 1:numel(dates)
    date = dates(cycle_dates).name;
    fprintf(1, 'Processing %s...', date);
    positions = dir([readpath '/' date]);
    positions = positions(~ismember({positions.name},{'.','..','.DS_Store'}));
    for cycle_positions = 1:numel(positions)
        position = positions(cycle_positions).name(1:end-4);
        if exist('bin_size','var')
            bin_size_prev = bin_size;
        end
        load([readpath '/' date '/' position '.mat']);
        if exist('bin_size_prev','var') && bin_size_prev ~= bin_size
            fprintf('Ensure bin sizes are equivalent!')
            return
        end
        odors = dir([odorpath '/' date '/' position]);
        odors = odors(~ismember({odors.name},{'.','..','.DS_Store'}));
        for cycle_odors = 1:numel(odors)
            odor = odors(cycle_odors).name(1:end-4);
            mean_data = eval([odor '_mean_data']);            
            mean_var_name = [odor '_mean_data_master'];
            if exist(mean_var_name,'var') && exist(rms_var_name,'var')
                mean_data = cat(1,eval([odor '_mean_data_master']),eval([odor '_mean_data']));
                eval([mean_var_name '= mean_data;'])
                if cycle_odors ~= 1
                    save(writepath, [odor '_mean_data_master'],'-append');
                else
                    experiment_master = cat(1,experiment_master,experiment);
                    save(writepath, [odor '_mean_data_master'],'experiment_master','-append');
                end
            elseif ~exist(mean_var_name,'var')
                eval([mean_var_name '= mean_data;'])                
                if cycle_odors ~= 1
                    save(writepath, [odor '_mean_data_master'],'-append');
                else 
                    experiment_master = experiment;
                    save(writepath, [odor '_mean_data_master'],'experiment_master','stim_on','stim_off','total_time','sample_rate','bin_size');
                end
            end
                        
%             switch Daywise
%                 case ismember(date,Day1)
%                     if ~exist([filepath '/Master_files/' filetype '_CC_Day1.mat'])
%                         experiment_master = experiment;
%                         save(writepath, [odor '_mean_data_master'],[odor '_rms_data_master'],'experiment_master','stim_on','stim_off','total_time','sample_rate','bin_size');
%                     elseif exist([filepath '/Master_files/' filetype '_CC_' select_day '.mat']) && ~ismember(date, who('-file', [filepath '/Master_files/' filetype '_CC_' select_day '.mat']))
%                         save(writepath, [odor '_mean_data_master'],[odor '_rms_data_master'],'-append');
%                     elseif exist([filepath '/Master_files/' filetype '_CC_' select_day '.mat']) && cycle_odors == 1
%                        
%                         mean_data = cat(1,eval([odor '_mean_data_master']),eval([odor '_mean_data']));
%                         eval([mean_var_name '= mean_data;'])
% 
%                         rms_data = cat(1,eval([odor '_rms_data_master']),eval([odor '_rms_data']));
%                         eval([rms_var_name '= rms_data;'])      
%                     end
%             end
        
        end
    end
    fprintf([num2str(cycle_dates/numel(dates)*100) '%% complete \n']);
end
fprintf(['Finished ' filename '_' exp_type '\n']);