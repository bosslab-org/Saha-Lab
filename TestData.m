function [bin_preds, trial_preds] = TestData(data)

bin_preds = zeros(size(data,4),size(data,4),size(data,3));    % preallocate array with zeros for padding
trial_preds = zeros(size(data,4),size(data,4));                     % preallocate array with zeros for padding 

for cycle_trials = 1:size(data,3)
    train_data_mean = permute(mean(data(:,:,1:size(data,3) ~= cycle_trials,:),3),[1 4 2 3]);
    test_data = permute(data(:,:,cycle_trials,:),[1 4 2 3]);        % testing data        
    for cycle_stimuli=1:size(data,4)
        bin_norm = squeeze(vecnorm(test_data(:,cycle_stimuli,:)-train_data_mean,2,1)); % 2-norm along dimension 1
        [~,bin_pred_class] = min(bin_norm,[],1);
        bin_preds(:,cycle_stimuli,cycle_trials) = accumarray(bin_pred_class',1);        
        
        trial_predst = mode(bin_pred_class);
        
        trial_preds(mode(bin_pred_class),cycle_stimuli) = trial_preds(mode(bin_pred_class),cycle_stimuli) + 1;
    end
end


bin_preds = sum(bin_preds,3)/(size(data,2)*size(data,3))*100;
trial_preds = trial_preds/size(data,3)*100;

% bin_class_acc = diag(sum(bin_pred_counts,3))/(size(data,2)*size(data,3)*size(data,4))*100;
% bin_total_acc = mean(bin_class_acc);

% trial_class_acc = (sum(trial_preds,2)/size(data,3))*100;
% trial_total_acc = mean(trial_class_acc);