function LDA(data, Colors, Odorants, LDA_title)

data = permute(reshape(data,size(data,1)/numel(Odorants),numel(Odorants),size(data,2)),[1 3 2]);

mu = mean(data,[1 3]);
data_norm = data-mu;
all_data_norm = reshape(permute(data_norm,[1 3 2]),size(data_norm,1)*size(data_norm,3),size(data_norm,2));
mu_norm = mean(data_norm,[1 3]);
odor_mu = mean(data_norm,1);

Sw = zeros(size(data,2));
Sb = zeros(size(data,2));
for cycle_odors = 1:numel(Odorants)    
    Sw = Sw + (data_norm(:,:,cycle_odors)-ones(size(data_norm,1),1)*odor_mu(:,:,cycle_odors))'*(data_norm(:,:,cycle_odors)-ones(size(data_norm,1),1)*odor_mu(:,:,cycle_odors));
    Sb = Sb+(odor_mu(:,:,cycle_odors)-mu_norm)'*(odor_mu(:,:,cycle_odors)-mu_norm);
end

[eigvec, eigval] = eig(pinv(Sw)*Sb,'vector');
[eigval, t] = sort(eigval,'descend');
% eigvec = real(eigvec(:,t));
var_exp = round(10000*eigval/sum(eigval))/100;
num_dims = 3;
LDA_vecs = eigvec(:,1:num_dims);
LDA_proj = all_data_norm*LDA_vecs;
LDA_projs = permute(reshape(LDA_proj,size(data_norm,1),size(data_norm,3),num_dims),[1 3 2]);

static_fig = figure('Units','normalized','Position',[0 0 0.5 0.39]);
static_axes = axes(static_fig);
hold on; grid on; box off; view(3);

for cycle_odors = 1:numel(Odorants)
    LDA_handle(cycle_odors) = scatter3(LDA_projs(:,1,cycle_odors),LDA_projs(:,2,cycle_odors),LDA_projs(:,3,cycle_odors),60,Colors(cycle_odors,:),'filled','MarkerEdgeColor','k');
end
title(static_axes,LDA_title,'FontSize',20,'FontWeight','Bold')

xlim(static_axes,[min(LDA_proj(:,1)) max(LDA_proj(:,1))])
ylim(static_axes,[min(LDA_proj(:,2)) max(LDA_proj(:,2))])
zlim(static_axes,[min(LDA_proj(:,3)) max(LDA_proj(:,3))])
xlabel(static_axes,['LDA1 (' num2str(var_exp(1)) '%)'],'FontSize',18,'FontWeight','Bold')
ylabel(static_axes,['LDA2 (' num2str(var_exp(2)) '%)'],'FontSize',18,'FontWeight','Bold')
zlabel(static_axes,['LDA3 (' num2str(var_exp(3)) '%)'],'FontSize',18,'FontWeight','Bold')
  
num_files = annotation(static_fig,'textbox', [0.2, 0.9, 0, 0], 'String', ['n = ' num2str(size(data,2))],'Units','normalized','Color','k','FontSize',14,'FontWeight','Bold','FitBoxToText','on','HorizontalAlignment','center');

[~, marker] = legend(LDA_handle,Odorants,'location','eastoutside','FontSize',16,'FontWeight','Bold');    
set(findobj(marker,'-property','MarkerSize'),'MarkerSize',16) 
return