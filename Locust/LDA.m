function LDA_fig = LDA(data, Colors, Odorants, LDA_title)

mu = mean(data,[2 3 4]);
class_mu = squeeze(mean(data,[2 3]));

odor_data_norm = squeeze(mean(data,3));
all_data = reshape(permute(odor_data_norm,[2 3 1]),size(data,2)*size(data,4),size(data,1));

Sw = zeros(size(data,1));
Sb = zeros(size(data,1));
for cycle_classes = 1:size(data,4)
    Sw = Sw + (odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))*(odor_data_norm(:,:,cycle_classes)-class_mu(:,cycle_classes))';
    Sb = Sb + (class_mu(:,cycle_classes)-mu)*(class_mu(:,cycle_classes)-mu)';
end
Sw = Sw/size(data,4);

[eigvec, eigval] = eig(pinv(Sw)*Sb,'vector');
[eigval, idx] = sort(eigval,'descend');
eigvec = real(eigvec(:,idx));
var_exp = round(10000*eigval/sum(eigval))/100;
LDs = eigvec(:,1:3);
LDA_proj = all_data*LDs;
LDA_projs = permute(reshape(LDA_proj,size(data,2),size(data,4),3),[1 3 2]);

LDA_fig = figure('Units','normalized','Position',[0 0 0.5 0.39]); LDA_axes = axes(LDA_fig);

hold on; grid on; box off; view(3);
for cycle_classes = 1:numel(Odorants)
    LDA_handle(cycle_classes) = scatter3(LDA_projs(:,1,cycle_classes),LDA_projs(:,2,cycle_classes),LDA_projs(:,3,cycle_classes),60,Colors(cycle_classes,:),'filled','MarkerEdgeColor','k');
end
xlim(gca,[min(LDA_proj(:,1)) max(LDA_proj(:,1))])
ylim(gca,[min(LDA_proj(:,2)) max(LDA_proj(:,2))])
zlim(gca,[min(LDA_proj(:,3)) max(LDA_proj(:,3))])
xlabel(LDA_axes,['LDA1 (' num2str(var_exp(1)) '%)'],'FontSize',18,'FontWeight','Bold')
ylabel(LDA_axes,['LDA2 (' num2str(var_exp(2)) '%)'],'FontSize',18,'FontWeight','Bold')
zlabel(LDA_axes,['LDA3 (' num2str(var_exp(3)) '%)'],'FontSize',18,'FontWeight','Bold')
title(LDA_axes,LDA_title,'FontSize',20,'FontWeight','Bold')
annotation(LDA_fig,'textbox', [0.2, 0.9, 0, 0], 'String', ['n = ' num2str(size(data,1))],'Units','normalized','Color','k','FontSize',14,'FontWeight','Bold','FitBoxToText','on','HorizontalAlignment','center');

[~, marker] = legend(LDA_handle,Odorants,'location','eastoutside','FontSize',16,'FontWeight','Bold');
set(findobj(marker,'-property','MarkerSize'),'MarkerSize',16) 
end
