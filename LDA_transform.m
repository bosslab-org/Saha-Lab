function [LDA_projs,LDA_means] = LDA_transform(data, files, lambda, LDA_dims)
mu = mean(data,[3 4]);
data_norm = data - mu;
odor_data_norm = squeeze(mean(data_norm,2));
mu = mean(odor_data_norm,[2 3]);
class_mu = squeeze(mean(odor_data_norm,2)); 

all_data = reshape(permute(odor_data_norm,[2 3 1]),size(odor_data_norm,2)*size(odor_data_norm,3),size(odor_data_norm,1));
Sw = zeros(size(data,1));
Sb = zeros(size(data,1));
for cycle_classes = 1:length(files)
    Sw = Sw + (odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))*(odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))';
    Sb = Sb + (class_mu(:,cycle_classes)-mu)*(class_mu(:,cycle_classes)-mu)';
end
Sw = Sw/length(files);  % normalize Sw based on the number of classes

mat_inv = ((1-lambda)*Sw+lambda*eye(size(Sw)))\Sb;
[eigvec, eigval] = eig(mat_inv,'vector');
[eigval, idx] = sort(eigval,'descend');
eigvec = real(eigvec(:,idx));
var_exp = round(10000*eigval/sum(eigval))/100;
LDs = eigvec(:,1:LDA_dims);
LDA_proj = all_data*LDs;
LDA_projs = permute(reshape(LDA_proj,size(odor_data_norm,2),length(files),size(LDs,2)),[2 1 3]);
LDA_means = mean(LDA_projs,2);







% 
%timespan = size(odor_data_norm,1);
% num_channels = size(odor_data_norm,2);
% num_trials = size(odor_data_norm,3);
% num_odors = size(odor_data_norm,4);
% odor_idx = 1:num_odors;
% 
% % odor_data_train = cell(1,num_odors);
% % for cycle_odors = 1:num_odors % condition data for LDA_train use
% %     odor_data_train{1,cycle_odors}=data(:,:,:,cycle_odors);
% % end
% 
% Xtrain = zeros(num_trials*num_odors*timespan,num_channels);
% Y = zeros(num_trials*num_odors,1);
% for cycle_odors = 1:num_odors
%     odor_data = squeeze(odor_data_norm(:,:,:,cycle_odors));
% %     tmp = odor_data_train{cycle_odors};
%     tmp2 = permute(odor_data,[2,1,3]);
%     tmp3 = reshape(tmp2,num_channels,[]);
%     tmp4 = tmp3';
%     idx_start = timespan*num_trials*(cycle_odors-1)+1;
%     idx = idx_start:idx_start+timespan*num_trials-1;
%     Xtrain(idx,:) = tmp4;
%     Y(idx,1) = odor_idx(cycle_odors);
% end
% mu = mean(Xtrain);
% X_norm = Xtrain-mu;

end