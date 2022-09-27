function rms_construct_tetrodes(stimulus, data, exp_pars, rms_window, smooth_window, bin_size, time, filepath, exp_path)

writepath = [filepath '/Master_files/']; % date '/Position_' position];
filename = ['R' num2str(rms_window) '_S' num2str(smooth_window) '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2))];

stim_on = exp_pars(1);
stim_off = exp_pars(2);
sample_rate = exp_pars(3);

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

    if exp==1 && ~exist([stimulus(1:end-3) '_RMS'])
        eval([stimulus(1:end-3) '_RMS = data;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'],'-append');
    elseif exp~=1 && exist([stimulus(1:end-3) '_RMS'])
        for cycle_tetrodes = 1:size(data,1)
            exp_new(cycle_tetrodes,:) = [exp_path ' Tetrode_' num2str(cycle_tetrodes)];
        end
        experiment = cat(1,experiment,exp_new);
        data = cat(1,eval([stimulus(1:end-3) '_RMS']),data);
        eval([stimulus(1:end-3) '_RMS = data;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'], 'experiment','-append');
    elseif exp==1 && exist([stimulus(1:end-3) '_RMS'])
        data = cat(1,eval([stimulus(1:end-3) '_RMS']),data);
        eval([stimulus(1:end-3) '_RMS = data;'])
        save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'], 'experiment','-append');
    end
else
    if ~exist(writepath)
        mkdir(writepath);
    end
    eval([stimulus(1:end-3) '_RMS = data;'])
    for cycle_tetrodes = 1:size(data,1)
        experiment(cycle_tetrodes,:) = [exp_path ' Tetrode_' num2str(cycle_tetrodes)];
    end
    save([writepath '/' filename '.mat'], [stimulus(1:end-3) '_RMS'],'experiment','stim_on','stim_off','sample_rate','bin_size');
end

fprintf([exp_path ' ' stimulus ' RMS-transformed. \n\n'])

return