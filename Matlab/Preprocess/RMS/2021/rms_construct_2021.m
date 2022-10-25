function rms_construct_2021(stimulus, data, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)

writepath = [filepath '/Master_files/']; % date '/Position_' position];
filename = ['R' num2str(rms_window) '_S' num2str(smooth_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];

stim_on = exp_pars(1);
stim_off = exp_pars(2);
sample_rate = exp_pars(5);

preprocess_data = @(x) movmean(sqrt(movmean(x.^ 2, rms_window, 4)),smooth_window,4);

% size(data) = (tetrodes,channels,trials,samples)
data = preprocess_data(data(:,:,:,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)) -...
    mean(preprocess_data(data(:,:,:,(stim_on-2)*sample_rate+1:stim_on*sample_rate)),[3 4]);                
data = permute(mean(reshape(permute(data,[4 3 1 2]),...
    (sample_rate*bin_size)/1000,size(data,4)/((sample_rate*bin_size)/1000),size(data,3),size(data,1),size(data,2)),[1 5]),[4 3 2 1]);  % takes mean of values within each bin and across tetrode channels
% size(data) = (tetrodes,trials,binned_samples)

if exist([writepath '/' filename '.mat'],'file') %check for file existence
    load([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'], 'experiment')
    exp = ismember([exp_path ' Tetrode_' num2str(size(data,1))], experiment, 'rows');

    if exp==1 && ~exist([stimulus(1:end-3) '_RMS'],'var')
        eval([stimulus(1:end-3) '_RMS = data;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'],'-append');
    elseif exp~=1 && exist([stimulus(1:end-3) '_RMS'],'var')
        for cycle_tetrodes = 1:size(data,1)
            exp_new(cycle_tetrodes,:) = [exp_path ' Tetrode_' num2str(cycle_tetrodes)];
        end
        experiment = cat(1,experiment,exp_new);
        data = cat(1,data,eval([stimulus(1:end-3) '_RMS']));
        eval([stimulus(1:end-3) '_RMS = data;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'], 'experiment','-append');
    elseif exp==1 && exist([stimulus(1:end-3) '_RMS'],'var')
        data = cat(1,eval([stimulus(1:end-3) '_RMS']),data);
        eval([stimulus(1:end-3) '_RMS = data;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'], 'experiment','-append');
    end
else
    if ~exist(writepath,'dir')
        mkdir(writepath);
    end
    eval([stimulus(1:end-3) '_RMS = data;'])
    for cycle_tetrodes = 1:size(data,1)
        experiment(cycle_tetrodes,:) = [exp_path ' Tetrode_' num2str(cycle_tetrodes)];
    end
    save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'],'experiment','stim_on','stim_off','sample_rate','bin_size');
end

fprintf([exp_path ' ' stimulus ' RMS-transformed. \n\n'])

% Day1 = {'05_25_2021';'06_24_2021';'06_29_2021';'07_20_2021';'08_03_2021';'08_10_2021'};
% Day2 = {'05_26_2021';'06_30_2021';'07_14_2021';'08_04_2021';'08_11_2021'};
% Day3 = {'05_27_2021';'06_26_2021';'07_01_2021';'07_15_2021';'07_22_2021';'08_05_2021';'08_12_2021'};
% Day4 = {'05_28_2021';'07_02_2021';'07_16_2021';'07_23_2021';'08_06_2021'};
% 
% for ii = 1:4
%     if ismember(date,strcat('Day',num2str(eval(num2str(2)))),'rows')
%         day = strcat('Day',num2str(eval(num2str(ii))))
%         break
%     end
% end
% 
% daypath = [writepath '/Day_files' filename Day(ii)];

return