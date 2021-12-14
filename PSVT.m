function PSVT(exp_pars,data,cycle_tets,cycle_chans,time)

PSVT_fig = figure('Units', 'normalized', 'Position', [0 0 .5 1]);
sgtitle(PSVT_fig,[exp_pars.Date ': ' exp_pars.Odorant '_Position_' exp_pars.Position ' Tetrode ' num2str(cycle_tets) ', Channel ', num2str(cycle_chans)], 'FontSize', 18, 'FontWeight', 'bold');

for cycle_tris = 1:size(data,3)
    subplot(size(data,3),1,cycle_tris); 
    hold on;
    data_temp = permute(data(:,:,cycle_tris,:),[3 4 1 2]);
    plot(1:(time(2)-time(1))*exp_pars.sample_rate, data_temp)

    stim_x = [(exp_pars.stim_on_trim-time(1))*exp_pars.sample_rate (exp_pars.stim_on_trim-time(1))*exp_pars.sample_rate (exp_pars.stim_off_trim-time(1))*exp_pars.sample_rate (exp_pars.stim_off_trim-time(1))*exp_pars.sample_rate];
    stim_y = [min(data_temp) max(data_temp) max(data_temp) min(data_temp)]; 
    stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0);

    filtered_trace = gca;
    filtered_trace.FontSize = 12;
    title(['Trial ' num2str(cycle_tris)], 'FontSize', 16);
    xlim([0 time(2)*exp_pars.sample_rate-time(1)*exp_pars.sample_rate]);
    xticks([0:exp_pars.sample_rate:length(data_temp)])
    xticklabels(0:time(end)-time(1));
    ylim([min(data_temp) max(data_temp)]);
    hold off;
end
PSVT_handle=axes(PSVT_fig,'visible','off'); 
PSVT_handle.XLabel.Visible='on';
PSVT_handle.YLabel.Visible='on';
ylabel(PSVT_handle,'Voltage (µV)', 'FontSize', 16);
xlabel(PSVT_handle,'Time (s)', 'FontSize', 16);
end