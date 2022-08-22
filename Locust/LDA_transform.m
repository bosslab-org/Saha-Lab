function [LDA_projs,LDA_means,LDs] = LDA_transform(data)

mu = mean(data,[2 3 4]);
class_mu = squeeze(mean(data,[2 3]));

odor_data_norm = squeeze(mean(data,2));
all_data = reshape(permute(odor_data_norm,[2 3 1]),size(data,3)*size(data,4),size(data,1));

Sw = zeros(size(data,1));
Sb = zeros(size(data,1));
for cycle_classes = 1:size(data,4)
    Sw = Sw + (odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))*(odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))';
    Sb = Sb + (class_mu(:,cycle_classes)-mu)*(class_mu(:,cycle_classes)-mu)';
end
Sw = Sw/size(data,4);

% mat_inv = ((1-lambda)*Sw+lambda*eye(size(Sw)))\Sb;
% mat_inv = pinv(Sw)*Sb;
[eigvec, eigval] = eig(pinv(Sw)*Sb,'vector');
[eigval, idx] = sort(eigval,'descend');
eigvec = real(eigvec(:,idx));
var_exp = round(10000*eigval/sum(eigval))/100;
LDs = eigvec(:,1:3);
LDA_proj = all_data*LDs;

LDA_projs = permute(reshape(LDA_proj,size(data,3),size(data,4),3),[3 1 2]);
LDA_means = squeeze(mean(LDA_projs,2));
end
