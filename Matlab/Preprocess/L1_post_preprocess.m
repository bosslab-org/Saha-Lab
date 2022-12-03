clear; close all;

load('/Users/alexanderfarnum/Documents/MATLAB/Master_files/R500_S500_l1_b10_-2to6.mat')
save('/Users/alexanderfarnum/Documents/MATLAB/Master_files/R500_S500_l1_b10_-2to6original.mat')

idx = [2,6,9,10,13,16:37,40,41];
% idx = [16:37,40,41];

experiment = experiment(idx,:);
Control_filt_RMS = Control_filt_RMS(idx,:,:);
Decane_filt_RMS = Decane_filt_RMS(idx,:,:);
Hexanal_filt_RMS = Hexanal_filt_RMS(idx,:,:);
Methylheptane_filt_RMS = Methylheptane_filt_RMS(idx,:,:);
Nonanal_filt_RMS = Nonanal_filt_RMS(idx,:,:);
Pentamethylheptane_filt_RMS = Pentamethylheptane_filt_RMS(idx,:,:);
Pentanal_filt_RMS = Pentanal_filt_RMS(idx,:,:);
Propylbenzene_filt_RMS = Propylbenzene_filt_RMS(idx,:,:);
Trichloroethylene_filt_RMS = Trichloroethylene_filt_RMS(idx,:,:);
Undecane_filt_RMS = Undecane_filt_RMS(idx,:,:);
clear idx ans
save('/Users/alexanderfarnum/Documents/MATLAB/Master_files/R500_S500_l1_b10_-2to6.mat')


% %% MAT
% clear; close all;
% load('/Users/alexanderfarnum/Documents/MATLAB/Locust/Odorants/MAT_files/08_20_2021/Position_2/Hexanal.mat')
% 
% for cycle_tetrodes = 1:size(Hexanal_data_filt,1)
%     figure;
%     for cycle_trials = 1:size(Hexanal_data_filt,3)
%         subplot(size(Hexanal_data_filt,3),1,cycle_trials);
%         plot(1:size(Hexanal_data_filt,4),squeeze(Hexanal_data_filt(cycle_tetrodes,1,cycle_trials,:)));
%         xticks(0:sample_rate:total_time*sample_rate)
%         xticklabels(0:total_time)
%     end
%     close;
% end



% %% RMS
% load('/Users/alexanderfarnum/Documents/MATLAB/Locust/Odorants/Master_files/R500_S500_b10_0to4.mat')
% 
% for cycle_tetrodes = 1:41
%     figure;
%     for cycle_trials = 1:size(Hexanal_RMS,2)
%         subplot(size(Hexanal_RMS,2),1,cycle_trials);
%         plot(1:size(Hexanal_RMS,3),squeeze(Hexanal_RMS(cycle_tetrodes,cycle_trials,:)));
%     end
%     close;
% end

