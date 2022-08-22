function [bin_preds, trial_preds] = TrainTestLocust(data,normtype,train_trials)

train_data = permute(mean(data(:,:,train_trials,:),3),[1 2 4 3]);
test_data = data(:,:,setxor(1:size(data,3),train_trials),:);

% preallocate arrays with zeros for padding
bin_preds = zeros(size(test_data,4),size(test_data,4),size(test_data,3));
trial_preds = zeros(size(test_data,4),size(test_data,4));

for cycle_trials = 1:size(test_data,3)    
    test_data_temp = permute(test_data(:,:,cycle_trials,:),[1 2 4 3]);
    for cycle_stimuli= 1:size(test_data,4)
        bin_norm = squeeze(vecnorm(test_data_temp(:,:,cycle_stimuli)-train_data,normtype,1)); % Specified norm along dimension 1
        [~,bin_pred_class] = min(bin_norm,[],2);
        bin_preds(:,cycle_stimuli,cycle_trials) = cat(1,accumarray(bin_pred_class,1),zeros(size(test_data,4)-max(bin_pred_class),1));
        trial_preds(mode(bin_pred_class),cycle_stimuli) = trial_preds(mode(bin_pred_class),cycle_stimuli) + 1;
    end
end

% total binwise predictions = number of time bins * number of trials
bin_preds = sum(bin_preds,3)/(size(test_data,2)*size(test_data,3))*100;

% total trialwise predictions = number of trials
trial_preds = trial_preds/size(test_data,3)*100;
