function RMS_PSVT(exp_pars,data,cycle_tets,cycle_chans,time)

RMS_PSVT_fig = figure('Units', 'normalized', 'Position', [.5 0 .5 1]);
sgtitle(RMS_PSVT_fig,[exp_pars.Date ': ' exp_pars.Odorant '_Position_' exp_pars.Position ' Tetrode ' num2str(cycle_tets) ', Channel ', num2str(cycle_chans)], 'FontSize', 18, 'FontWeight', 'bold');

for cycle_tris = 1:size(data,3)
    subplot(size(data,3),1,cycle_tris); 
    hold on;
    data_temp = permute(data(:,:,cycle_tris,:),[3 4 1 2]);
    plot(1:(time(2)-time(1))*exp_pars.bins_per_sec, data_temp)
    
    stim_x = [(exp_pars.stim_on_trim-time(1))*exp_pars.bins_per_sec (exp_pars.stim_on_trim-time(1))*exp_pars.bins_per_sec (exp_pars.stim_off_trim-time(1))*exp_pars.bins_per_sec (exp_pars.stim_off_trim-time(1))*exp_pars.bins_per_sec];
    stim_y = [min(data_temp) max(data_temp) max(data_temp) min(data_temp)]; 
    stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);
    
    filtered_trace = gca;
    filtered_trace.FontSize = 12;
    title(['Trial ' num2str(cycle_tris)], 'FontSize', 16);
    xlim([0 (time(2)-time(1))*exp_pars.bins_per_sec]);
    xticks([0:exp_pars.bins_per_sec:length(data_temp)])
    xticklabels(0:time(end)-time(1));
    ylim([min(data_temp) max(data_temp)]);
    hold off;
end
RMS_PSVT_han=axes(RMS_PSVT_fig,'visible','off'); 
RMS_PSVT_han.XLabel.Visible='on';
RMS_PSVT_han.YLabel.Visible='on';
ylabel(RMS_PSVT_han,'Voltage (µV)', 'FontSize', 16);
xlabel(RMS_PSVT_han,'Time (s)', 'FontSize', 16);
end