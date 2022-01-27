function static_fig = LDA_v2(data, Colors, Odorants, files, lambda, LDA_title)
mu = mean(data,[2 3]);
class_mu = squeeze(mean(data,2)); 

all_data = reshape(permute(data,[2 3 1]),size(data,2)*size(data,3),size(data,1));
Sw = zeros(size(data,1));
Sb = zeros(size(data,1));
for cycle_classes = 1:length(files)
    Sw = Sw + (data(:,:,cycle_classes)-class_mu(:,cycle_classes))*(data(:,:,cycle_classes)-class_mu(:,cycle_classes))';
    Sb = Sb + (class_mu(:,cycle_classes)-mu)*(class_mu(:,cycle_classes)-mu)';
end
Sw = Sw/length(files);  % normalize Sw based on the number of classes


%%%%%LDA-Nearest Centroid Interpolation: Ridge Regularization (for overfitting/high variance)%%%%%
figure;
num_graphs = 6; % must be even number
for cycle_lambda = 1:num_graphs
    subplot(2,num_graphs/2,cycle_lambda)
    test_lambda = (cycle_lambda-1)/(num_graphs-1);     % Vanilla LDA when lambda is zero; spherical LDA when lambda is one
    mat_inv = ((1-test_lambda)*Sw+test_lambda*eye(size(Sw)))\Sb;
    [eigvec, eigval] = eig(mat_inv,'vector');
    [eigval, idx] = sort(eigval,'descend');
    eigvec = real(eigvec(:,idx));
    LDs = eigvec(:,1:3);
    LDA_proj = all_data*LDs;
    LDA_projs = permute(reshape(LDA_proj,size(data,2),length(files),size(LDs,2)),[1 3 2]);

    hold on; grid on; box off; view(3);
    for cycle_classes = 1:length(files)
        LDA_handle(cycle_classes) = scatter3(LDA_projs(:,1,cycle_classes),LDA_projs(:,2,cycle_classes),LDA_projs(:,3,cycle_classes),60,Colors(cycle_classes,:),'filled'); %,'MarkerEdgeColor','k');
    end
    xlim(gca,[min(LDA_proj(:,1)) max(LDA_proj(:,1))])
    ylim(gca,[min(LDA_proj(:,2)) max(LDA_proj(:,2))])
    zlim(gca,[min(LDA_proj(:,3)) max(LDA_proj(:,3))])
    title(['Lambda = ' num2str(test_lambda)]);
    hold off;
end
sgtitle('0: Vanilla LDA, 1: Spherical LDA','FontSize',20,'FontWeight','Bold');


% Vanilla LDA when lambda is zero; spherical LDA when lambda is one
mat_inv = ((1-lambda)*Sw+lambda*eye(size(Sw)))\Sb;
[eigvec, eigval] = eig(mat_inv,'vector');
[eigval, idx] = sort(eigval,'descend');
eigvec = real(eigvec(:,idx));
var_exp = round(10000*eigval/sum(eigval))/100;
LDs = eigvec(:,1:3);
LDA_proj = all_data*LDs;
LDA_projs = permute(reshape(LDA_proj,size(data,2),length(files),size(LDs,2)),[1 3 2]);

static_fig = figure('Units','normalized','Position',[0 0 0.5 0.39]);
static_axes = axes(static_fig);
hold on; grid on; box off; view(3);
for cycle_classes = 1:length(files)
    LDA_handle(cycle_classes) = scatter3(LDA_projs(:,1,cycle_classes),LDA_projs(:,2,cycle_classes),LDA_projs(:,3,cycle_classes),60,Colors(cycle_classes,:),'filled'); %,'MarkerEdgeColor','k');
end
xlim(static_axes,[min(LDA_proj(:,1)) max(LDA_proj(:,1))])
ylim(static_axes,[min(LDA_proj(:,2)) max(LDA_proj(:,2))])
zlim(static_axes,[min(LDA_proj(:,3)) max(LDA_proj(:,3))])
xlabel(static_axes,['LDA1 (' num2str(var_exp(1)) '%)'],'FontSize',18,'FontWeight','Bold')
ylabel(static_axes,['LDA2 (' num2str(var_exp(2)) '%)'],'FontSize',18,'FontWeight','Bold')
zlabel(static_axes,['LDA3 (' num2str(var_exp(3)) '%)'],'FontSize',18,'FontWeight','Bold')
title(static_axes,LDA_title,'FontSize',20,'FontWeight','Bold')
num_files = annotation(static_fig,'textbox', [0.2, 0.9, 0, 0], 'String', ['n = ' num2str(size(data,1))],'Units','normalized','Color','k','FontSize',14,'FontWeight','Bold','FitBoxToText','on','HorizontalAlignment','center');

% [~, marker] = legend(LDA_handle,Odorants,'location','eastoutside','FontSize',16,'FontWeight','Bold');
% set(findobj(marker,'-property','MarkerSize'),'MarkerSize',16) 

end