function static_fig = PsthSem(data, Colors, stimuli, smooth, PCA_title, time, bin_size)

sem = @(x) std(x,0,2)/size(x,2);
sample_rate = 20000;
static_fig = figure; hold on;
for cycle_stimuli = 1:numel(stimuli)
    data_temp = permute(mean(data(:,:,:,cycle_stimuli)),[2 3 1]);

    spline_interpolates = linspace(0,size(data_temp,1),size(data_temp,1)*smooth);    % calculate number of spline interpolates
    spline_mean = makima(1:size(mean(data_temp,2),1),mean(data_temp,2),spline_interpolates);
    spline_sem = makima(1:size(data_temp,1),sem(data_temp),spline_interpolates);
    psvt_plot(cycle_stimuli) = plot(1:size(data_temp,1)*smooth,spline_mean,'Color',Colors(cycle_stimuli,:));
    patch([1:size(data_temp,1)*smooth, fliplr(1:size(data_temp,1)*smooth)], [(spline_mean - spline_sem), fliplr(spline_mean + spline_sem)], Colors(cycle_stimuli,:),'FaceAlpha',0.25,'LineStyle','none');   

    max_psvt(cycle_stimuli) = max(spline_mean + spline_sem,[],'all');
end
xtick_res = 0.5;
xticks(0:1000/bin_size*xtick_res*smooth:size(data,2)*smooth)
xticklabels(time(1):xtick_res:time(2))
% legend(psvt_plot,stimuli);
ylim([0 max(max_psvt)]);

bps = 1000/bin_size;

xlim([0 size(data,2)*smooth])
% xlim([1.5*bps*smooth+1 3.5*bps*smooth])
% xlim([5.5*bps*smooth+1 7.5*bps*smooth])


% stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (4+abs(time(1)))*sample_rate (4+abs(time(1)))*sample_rate];
% stim_y = [min(ylim) max(ylim) max(ylim) min(ylim)];
% stim_patch = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0); 
return