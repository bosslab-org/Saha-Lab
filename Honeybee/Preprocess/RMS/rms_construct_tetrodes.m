function rms_construct_tetrodes(stimulus, data_filt, data_processed, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)

writepath = [filepath '/Master_files/'];
filename = ['R' num2str(rms_window) '_S' num2str(smooth_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];

stim_on = exp_pars(1);
stim_off = exp_pars(2);
sample_rate = exp_pars(3);

preprocess_data = @(x) movmean(sqrt(movmean(x.^ 2, rms_window, 4)),smooth_window,4);
% input size = (tetrodes,channels,trials,samples)
% output size = (tetrodes,trials,binned_samples)

data_filt = preprocess_data(data_filt(:,:,:,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)) -...
    mean(preprocess_data(data_filt(:,:,:,(stim_on-2)*sample_rate+1:stim_on*sample_rate)),[3 4]);                
data_filt = permute(mean(reshape(permute(data_filt,[4 3 1 2]),...
    (sample_rate*bin_size)/1000,size(data_filt,4)/((sample_rate*bin_size)/1000),size(data_filt,3),size(data_filt,1),size(data_filt,2)),[1 5]),[4 3 2 1]);  % takes mean of values within each bin and across tetrode channels

data_processed = preprocess_data(data_processed(:,:,:,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)) -...
    mean(preprocess_data(data_processed(:,:,:,(stim_on-2)*sample_rate+1:stim_on*sample_rate)),[3 4]);                
data_processed = permute(mean(reshape(permute(data_processed,[4 3 1 2]),...
    (sample_rate*bin_size)/1000,size(data_processed,4)/((sample_rate*bin_size)/1000),size(data_processed,3),size(data_processed,1),size(data_processed,2)),[1 5]),[4 3 2 1]);  % takes mean of values within each bin and across tetrode channels

if exist([writepath '/' filename '.mat'],'file') %check for file existence
    load([writepath '/' filename '.mat'], [stimulus(1:end-3) '_filt_RMS'], [stimulus(1:end-3) '_processed_RMS'], 'experiment')
    exp = ismember([exp_path ' Tetrode_' num2str(size(data_filt,1))], experiment, 'rows');

    if exp==1 && ~exist([stimulus(1:end-3) '_filt_RMS'],'var')
        eval([stimulus(1:end-3) '_filt_RMS = data_filt;'])
        eval([stimulus(1:end-3) '_processed_RMS = data_processed;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_filt_RMS'], [stimulus(1:end-3) '_processed_RMS'],'-append');
    elseif exp~=1 && exist([stimulus(1:end-3) '_filt_RMS'],'var')
        for cycle_tetrodes = 1:size(data_filt,1)
            exp_new(cycle_tetrodes,:) = [exp_path ' Tetrode_' num2str(cycle_tetrodes)];
        end
        experiment = cat(1,experiment,exp_new);
        data_filt = cat(1,eval([stimulus(1:end-3) '_filt_RMS']),data_filt);
        eval([stimulus(1:end-3) '_filt_RMS = data_filt;'])
        data_processed = cat(1,eval([stimulus(1:end-3) '_processed_RMS']),data_processed);
        eval([stimulus(1:end-3) '_processed_RMS = data_processed;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_filt_RMS'], [stimulus(1:end-3) '_processed_RMS'], 'experiment','-append');
    elseif exp==1 && exist([stimulus(1:end-3) '_filt_RMS'],'var')       
        data_filt = cat(1,eval([stimulus(1:end-3) '_filt_RMS']),data_filt);
        eval([stimulus(1:end-3) '_filt_RMS = data_filt;'])
        data_processed = cat(1,eval([stimulus(1:end-3) '_processed_RMS']),data_processed);
        eval([stimulus(1:end-3) '_processed_RMS = data_processed;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_filt_RMS'], [stimulus(1:end-3) '_processed_RMS'], 'experiment','-append');
    end
else
    if ~exist(writepath,'dir')
        mkdir(writepath);
    end
    eval([stimulus(1:end-3) '_filt_RMS = data_filt;'])
    eval([stimulus(1:end-3) '_processed_RMS = data_processed;'])
    for cycle_tetrodes = 1:size(data_filt,1)
        experiment(cycle_tetrodes,:) = [exp_path ' Tetrode_' num2str(cycle_tetrodes)];
    end
    save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_filt_RMS'], [stimulus(1:end-3) '_processed_RMS'],'experiment','stim_on','stim_off','sample_rate','bin_size');
end
fprintf([exp_path ' ' stimulus ' RMS-transformed. \n\n'])
return