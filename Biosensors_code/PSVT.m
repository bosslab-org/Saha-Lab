function PSVT(readpath,date,position,tetrode,channel,trial,colors,order)

time = [-1 6];

filepath = ['Position_' num2str(position)];
files = dir([readpath '/' date '/' filepath]);
files = files(~ismember({files.name},{'.','..','.DS_Store'}));
files = files(order);

PSVT_fig = figure('Units','normalized','Position',[0 0 0.75 1]);
set(0, 'DefaultTextInterpreter', 'none');    % format underscores in title to avoid subscripting
sgtitle(date);
for cycle_stimuli = 1:numel(files)
    load([readpath '/' date '/' filepath '/' files(cycle_stimuli).name]);
    
    data = eval([files(cycle_stimuli).name(1:end-4) '_data_filt']);
    data = squeeze(data(tetrode,channel,trial(cycle_stimuli),(stim_on+time(1))*sample_rate+1:(stim_on+time(2))*sample_rate));
    
    plots(cycle_stimuli) = subplot(numel(files),1,cycle_stimuli);
    PSVT_subplots(cycle_stimuli) = plot(1:(time(2)-time(1))*sample_rate,data,'Color',colors(cycle_stimuli,:));
     
    xlim([0 (time(2)-time(1))*sample_rate])
    
    if cycle_stimuli ~= numel(files)
        set(gca,'XTick',[])
    else
        xticks(0:sample_rate:(time(2)-time(1))*sample_rate)
        xticklabels(time(1):time(2));
    end
    
    ylim_min(cycle_stimuli) = min(data);
    ylim_max(cycle_stimuli) = max(data);

    ylim([min(ylim_min) max(ylim_max)])
    
    stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
    stim_y = [min(ylim_min(cycle_stimuli)) max(ylim_max(cycle_stimuli)) max(ylim_max(cycle_stimuli)) min(ylim_min(cycle_stimuli))];
    stim_patch(cycle_stimuli) = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);    

    
    title(files(cycle_stimuli).name(1:end-4))
    linkaxes(plots,'y')
    linkprop(stim_patch,{'XData','YData'});
end

PSVT_axes = axes(PSVT_fig,'visible','off');
PSVT_axes.YLabel.Visible='on';
ylabel(PSVT_axes, 'Voltage (µV)', 'FontSize', 16);
PSVT_axes.XLabel.Visible='on';
xlabel(PSVT_axes, 'Time (s)', 'FontSize', 16);
