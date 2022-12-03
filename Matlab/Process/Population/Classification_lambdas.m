function lambda_fig = Classification_lambdas(data, lambda, lambda_res)

mu = mean(data,[2 3 4]);
class_mu = squeeze(mean(data,[2 3]));

odor_data_norm = squeeze(mean(data,3));
all_data = reshape(permute(odor_data_norm,[2 3 1]),size(data,2)*size(data,4),size(data,1));

lambda_fig = figure('Units','normalized','Position',[0 0 0.75 0.75]); LDA_axes = axes(lambda_fig);
for cycle_lambda = 1:lambda*lambda_res
    Sw = zeros(size(data,1));
    Sb = zeros(size(data,1));
    for cycle_classes = 1:size(data,4)
        Sw = Sw + (odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))*(odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))';
        Sb = Sb + (class_mu(:,cycle_classes)-mu)*(class_mu(:,cycle_classes)-mu)';
    end
    Sw = Sw/size(data,4);

    mat_inv = ((1-cycle_lambda/lambda)*Sw+cycle_lambda/lambda*eye(size(Sw)))\Sb;
    [eigvec, eigval] = eig(mat_inv,'vector');
    [eigval, idx] = sort(eigval,'descend');
    eigvec = real(eigvec(:,idx));
    var_exp = round(10000*eigval/sum(eigval))/100;
    LDs = eigvec(:,1:3);
    LDA_proj = all_data*LDs;
    LDA_projs = permute(reshape(LDA_proj,size(data,2),size(data,4),3),[1 3 2]);

%     Generate plots
%     for cycle_classes = 1:numel(Odorants)
%         LDA_handle(cycle_classes) = scatter3(LDA_projs(:,1,cycle_classes),LDA_projs(:,2,cycle_classes),LDA_projs(:,3,cycle_classes),60,Colors(cycle_classes,:),'filled','MarkerEdgeColor','k');
%     end
    
bin_preds = zeros(size(data,4),size(data,4),size(data,3)); % zero padded array
trial_preds = zeros(size(data,4),size(data,4)); % zero padded array

for cycle_trials = 1:size(data,3)
    
%     Time bin average was used previously
    train_data_mean = permute(mean(data(:,:,1:size(data,3) ~= cycle_trials,:),[2,3]),[1 2 4 3]);

    mean(LDA_projs,1)

%     No time bin average seems to give far better results
%     train_data_mean = permute(mean(data(:,:,1:size(data,3) ~= cycle_trials,:),3),[1 2 4 3]);
    
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

end

return