function rms_PSVT(readpath,date,position,tetrode,channel,trial,bin_mult,colors,order)

time = [-1 6];
rms_window = 500;
smooth_window = 500;

filepath = ['Position_' num2str(position)];
files = dir([readpath '/' date '/' filepath]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));
files = files(order);
colors = colors(order,:);

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
        bin_mult,size(data,4)/bin_mult,size(data,3),size(data,1),size(data,2)),1),[4 5 3 2 1]); % takes mean of values within each bin
    % size(data) = (tetrodes,channels,trials,binned_samples)

    plots(cycle_stimuli) = subplot(numel(files),1,cycle_stimuli);
    plot(1:size(data,4),squeeze(data(tetrode,channel,trial(cycle_stimuli),:)),'Color',colors(cycle_stimuli,:),'LineWidth',2);
      
    xlim([0 size(data,4)])
    xticks(0:size(data,4)/(time(2)-time(1)):size(data,4));
    if cycle_stimuli ~= numel(files)
        set(gca,'XTick',[])
    else
        xticklabels(time(1):time(2));
    end

    title(files(cycle_stimuli).name(1:end-4))
    linkaxes(plots,'y')
end

rms_PSVT_axes = axes(rms_PSVT_fig,'visible','off');
rms_PSVT_axes.YLabel.Visible='on';
ylabel(rms_PSVT_axes, 'Voltage (µV)', 'FontSize', 16);
rms_PSVT_axes.XLabel.Visible='on';
xlabel(rms_PSVT_axes, 'Time (s)', 'FontSize', 16);