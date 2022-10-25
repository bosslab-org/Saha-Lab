clear; close all;

load('/Users/alexanderfarnum/Documents/MATLAB/Locust/Odorants/Master_files/R500_S500_b10_0to2.mat')
save('/Users/alexanderfarnum/Documents/MATLAB/Locust/Odorants/Master_files/R500_S500_b10_0to2original.mat')

% idx = [2,6,9,10,13,16:37,40,41];
idx = [6:29];

experiment = experiment(idx,:);
Control_RMS = Control_RMS(idx,:,:);
Decane_RMS = Decane_RMS(idx,:,:);
Hexanal_RMS = Hexanal_RMS(idx,:,:);
Methylheptane_RMS = Methylheptane_RMS(idx,:,:);
Nonanal_RMS = Nonanal_RMS(idx,:,:);
Pentamethylheptane_RMS = Pentamethylheptane_RMS(idx,:,:);
Pentanal_RMS = Pentanal_RMS(idx,:,:);
Propylbenzene_RMS = Propylbenzene_RMS(idx,:,:);
Trichloroethylene_RMS = Trichloroethylene_RMS(idx,:,:);
Undecane_RMS = Undecane_RMS(idx,:,:);
clear idx ans
save('/Users/alexanderfarnum/Documents/MATLAB/Master_files/R500_S500_l1_b10_0to4.mat')


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

