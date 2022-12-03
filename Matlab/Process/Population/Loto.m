function [bin_preds, trial_preds] = Loto(data,normtype)
bin_preds = zeros(size(data,4),size(data,4),size(data,3)); % zero padded array
trial_preds = zeros(size(data,4),size(data,4)); % zero padded array

for cycle_trials = 1:size(data,3)
    
%     Time bin average was used previously
%     train_data_mean = permute(mean(data(:,:,1:size(data,3) ~= cycle_trials,:),[2,3]),[1 2 4 3]);

%     No time bin average seems to give far better results
    train_data_mean = permute(mean(data(:,:,1:size(data,3) ~= cycle_trials,:),3),[1 2 4 3]);
    
    test_data = permute(data(:,:,cycle_trials,:),[1 2 4 3]);        % testing data
    for cycle_stimuli= 1:size(data,4)
        bin_norm = squeeze(vecnorm(test_data(:,:,cycle_stimuli)-train_data_mean,normtype,1)); % Specified norm along dimension 1
        [~,bin_pred_class] = min(bin_norm,[],2);
        bin_preds(:,cycle_stimuli,cycle_trials) = cat(1,accumarray(bin_pred_class,1),zeros(size(data,4)-max(bin_pred_class),1));
        trial_preds(mode(bin_pred_class),cycle_stimuli) = trial_preds(mode(bin_pred_class),cycle_stimuli) + 1;
    end
end

% total binwise predictions = number of time bins * number of trials
bin_preds = sum(bin_preds,3)/(size(data,2)*size(data,3))*100;

% total trialwise predictions = number of trials
trial_preds = trial_preds/size(data,3)*100;