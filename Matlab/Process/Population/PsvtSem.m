clear; clc; close all;
exp_type = 'l1';            % l1: Locust 1% biomarkers  
                            % lc: Locust cell culture
                            % h1: Honeybee 1% biomarkers  
                            % hb: Honeybee breath mixture

RMS_window = 500;
smooth_window = 500;
new_bin_size = 50;
bin_size = 10;
time = [0 4];
smooth = 3;

basepath = '/Users/alexanderfarnum/Documents/MATLAB';
filename = ['R' num2str(RMS_window) '_S' num2str(RMS_window) '_' exp_type '_b' num2str(bin_size) '_' num2str(time(1)) 'to' num2str(time(2)) '.mat']; 
load([basepath '/Master_files/' filename]);
stimuli = who('*_RMS');
bin_mult = new_bin_size/bin_size;

Colors = [242,147,147 ; 189,67,67 ; 237,70,47 ; 165,245,167 ; 57,123,39 ;...
    69,148,39 ; 79,201,251 ; 49,74,251 ; 147,147,147 ; 38,38,38]./255;

sem = @(x) std(x,0,2)/size(x,2);

figure; hold on;
for cycle_stimuli = 1:numel(stimuli)
    data_filt = eval(stimuli{cycle_stimuli});
    data_filt = permute(mean(reshape(data_filt,size(data_filt,1),size(data_filt,2),bin_mult,size(data_filt,3)/bin_mult),[1 3]),[4 2 1 3]);
    
%     plot(1:size(data_filt,1),mean(data_filt,2),'Color',Colors(cycle_stimuli,:));
%     patch([1:size(data_filt,1), fliplr(1:size(data_filt,1))], [(mean(data_filt,2) - sem(data_filt)); flipud(mean(data_filt,2) + sem(data_filt))], Colors(cycle_stimuli,:),'FaceAlpha',0.25,'LineStyle','none');   

    spline_interpolates = linspace(0,size(data_filt,1),size(data_filt,1)*smooth);    % calculate number of spline interpolates
    spline_mean = makima(1:size(mean(data_filt,2),1),mean(data_filt,2),spline_interpolates);
    spline_sem = makima(1:size(data_filt,1),sem(data_filt),spline_interpolates);
    plot(1:size(data_filt,1)*smooth,spline_mean,'Color',Colors(cycle_stimuli,:));
    patch([1:size(data_filt,1)*smooth, fliplr(1:size(data_filt,1)*smooth)], [(spline_mean - spline_sem), fliplr(spline_mean + spline_sem)], Colors(cycle_stimuli,:),'FaceAlpha',0.25,'LineStyle','none');   
end
xtick_res = 0.5;
xticks(0:1000/new_bin_size*xtick_res*smooth:size(data_filt,1)*smooth)
xticklabels(0:xtick_res:4)
% legend(stimuli)
%     stim_x = [abs(time(1))*sample_rate abs(time(1))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate (stim_off-stim_on+abs(time(1)))*sample_rate];
%     stim_y = [min(ylim_min(cycle_stimuli)) max(ylim_max(cycle_stimuli)) max(ylim_max(cycle_stimuli)) min(ylim_min(cycle_stimuli))];
%     stim_patch(cycle_stimuli) = patch(stim_x,stim_y,'k', 'FaceAlpha', 0.05, 'EdgeAlpha', 0.0); 
