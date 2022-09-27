function [bin_preds, trial_preds, train_stimuli, test_stimuli] = TtBee(data,normtype,stimuli,train_ext)

% train_stimuli = stimuli(contains(stimuli,train_ext));
% test_stimuli = stimuli(~contains(stimuli,train_ext));

% % With mineral oil
% train_idx = [1 4 7 9];
% test_idx = [2 5 8 10];

% Without mineral oil
train_idx = [1 4 7];
test_idx = [2 5 8];

train_stimuli = stimuli(train_idx,:);
test_stimuli = stimuli(test_idx,:);

%     Time bin average was used previously
% train_data = permute(mean(data(:,:,:,train_idx),[2 3]),[1 2 4 3]);

%     No time bin average seems to give far better results
train_data = permute(mean(data(:,:,:,train_idx),3),[1 2 4 3]);


test_data = permute(data(:,:,:,test_idx),[1 2 4 3]);

bin_preds = zeros(size(test_data,3),size(test_data,3),size(test_data,4));    % preallocate array with zeros for padding
trial_preds = zeros(size(test_data,3),size(test_data,3));                    % preallocate array with zeros for padding 

for cycle_trials = 1:size(test_data,4)    
    test_data_temp = test_data(:,:,:,cycle_trials);        % testing data
    for cycle_stimuli= 1:size(test_data,3)
        bin_norm = squeeze(vecnorm(test_data_temp(:,:,cycle_stimuli)-train_data,normtype,1)); % Specified norm along dimension 1
        [~,bin_pred_class] = min(bin_norm,[],2);
        bin_preds(:,cycle_stimuli,cycle_trials) = cat(1,accumarray(bin_pred_class,1),zeros(size(test_data,3)-max(bin_pred_class),1));
        trial_preds(mode(bin_pred_class),cycle_stimuli) = trial_preds(mode(bin_pred_class),cycle_stimuli) + 1;
    end
end

% total binwise predictions = number of time bins * number of trials
bin_preds = sum(bin_preds,3)/(size(data,2)*size(data,3))*100;

% total trialwise predictions = number of trials
trial_preds = trial_preds/size(data,3)*100;