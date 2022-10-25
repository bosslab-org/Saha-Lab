function rms_PSVT(readpath,date,position,tetrode,channel,trial)

time = [-1 6];
rms_window = 500;
smooth_window = 500;
bin_size = 50;

filepath = ['Position_' num2str(position)];
files = dir([readpath '/' date '/' filepath]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));
sample_rate = 20000;

colors = [246,146,30 ; 255,0,255 ; 0,200,0 ; 150,150,150 ; 0,0,0 ; 255,0,0 ; 117,76,36]/255;

preprocess_data = @(x) movmean(sqrt(movmean(x.^ 2, rms_window, 4)),smooth_window,4);

rms_PSVT_fig = figure('Units','normalized','Position',[0 0 0.75 1]);
set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
sgtitle(date);
for cycle_stimuli = 1:numel(files)
    load([readpath '/' date '/' filepath '/' files(cycle_stimuli).name]);
    
    data = eval([files(cycle_stimuli).name(1:end-4) '_data_filt']);
    
    % size(data) = (tetrodes,channels,trials,samples)
    data = preprocess_data(data(:,:,:,(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate)) -...
        mean(preprocess_data(data(:,:,:,(stim_on-2)*sample_rate+1:stim_on*sample_rate)),[3 4]);
    
    data = permute(mean(reshape(permute(data,[4 3 1 2]),...
        (sample_rate*bin_size)/1000,size(data,4)/((sample_rate*bin_size)/1000),size(data,3),size(data,1),size(data,2)),1),[4 5 3 2 1]); % takes mean of values within each bin
    % size(data) = (tetrodes,channels,trials,binned_samples)
    
    plots(cycle_stimuli) = subplot(numel(files),1,cycle_stimuli);
    rms_PSVT_subplots(cycle_stimuli) = plot(1:size(data,4),squeeze(data(tetrode,channel,trial(cycle_stimuli),:)),'Color',colors(cycle_stimuli,:),'LineWidth',2);
     
%     xlim([0 (time(2)-time(1))*sample_rate])
    
    if cycle_stimuli ~= numel(files)
        set(gca,'XTick',[])
    else
        xticks([0:size(data,4):size(data,4)/(abs(time(1))+abs(time(2)))]);
        xticklabels(time(1):time(2));
    end
%     
%     ylim_min(cycle_stimuli) = min(data);
%     ylim_max(cycle_stimuli) = max(data);
%     
%     stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
%     stim_y = [min(ylim_min(cycle_stimuli)) max(ylim_max(cycle_stimuli)) max(ylim_max(cycle_stimuli)) min(ylim_min(cycle_stimuli))];
%     stim_patch(cycle_stimuli) = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);    
% 
%     
%     title(files(cycle_stimuli).name(1:end-4))
%     linkaxes(plots,'y')
    %link = linkprop(stim_patch);
end

% rms_PSVT_axes = axes(rms_PSVT_fig,'visible','off');
% rms_PSVT_axes.YLabel.Visible='on';
% ylabel(rms_PSVT_axes, 'Voltage (µV)', 'FontSize', 16);
% rms_PSVT_axes.XLabel.Visible='on';
% xlabel(rms_PSVT_axes, 'Time (s)', 'FontSize', 16);