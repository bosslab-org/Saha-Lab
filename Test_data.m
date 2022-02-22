function [trial_pred_counts, bin_pred_counts] = Test_data(test_data,train_data_mean,pred_counts_zeros,norm_type)

if norm_type == 1
    pred_data = squeeze(vecnorm(test_data-permute(train_data_mean,[3 1 2]),1,2));    % compute 1-norm for all time bins
elseif norm_type == 2
    pred_data = squeeze(vecnorm(test_data-permute(train_data_mean,[3 1 2]),2,2));    % compute 2-norm for all time bins
end
[~,bin_pred_class] = min(pred_data,[],2);

bin_pred_counts_temp = accumarray(bin_pred_class,1)';
bin_pred_counts = cat(2,bin_pred_counts_temp,zeros(1,size(pred_counts_zeros,2)-size(bin_pred_counts_temp,2)));

trial_pred_counts = mode(bin_pred_class);

return