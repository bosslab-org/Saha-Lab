function [HD_trial_pred_counts_temp, HD_bin_pred_counts] = Test_data(test_data,train_data_mean, files, pred_counts_zeros, cycle_classes)

pred_data = squeeze(vecnorm(test_data(cycle_classes,:,:)-permute(train_data_mean,[1 3 2]),2,3))';    % compute L2 norm for all time bins

[~,HD_bin_pred_class] = min(pred_data,[],2);

HD_bin_pred_counts_temp = accumarray(HD_bin_pred_class,1)';
if numel(HD_bin_pred_counts_temp) ~= numel(files)
    HD_bin_pred_counts_temp(numel(pred_counts_zeros)) = 0;     % pad vector with 0 events if necessary
end
HD_bin_pred_counts = HD_bin_pred_counts_temp;

HD_trial_pred_counts_temp = mode(HD_bin_pred_class);