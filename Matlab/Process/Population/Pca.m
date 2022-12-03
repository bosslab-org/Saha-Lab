function static_fig = Pca(data, Colors, Odorants, smooth, PCA_title, time, new_bin_size, bps, line_mult)

data_norm = reshape(permute(mean(data,3),[2,4,1,3]),size(data,2)*size(data,4),size(data,1));
[~,scores,~,~,var_exp] = pca(data_norm);

my_filter = designfilt('lowpassiir','FilterOrder',smooth,'HalfPowerFrequency',0.15,'DesignMethod','butter'); % design custom filter
scores = filtfilt(my_filter,scores);
scores = permute(reshape(scores,size(scores,1)/numel(Odorants),numel(Odorants),size(scores,2)),[3 1 2]);
scores_norm = scores-scores(:,1,:);

static_fig = figure('Units','normalized','Position',[0 0.5 0.5 0.39]);  % generate single figure for all Odorant PCA trajectories
static_axes = axes(static_fig);
hold on; grid on; box off; view(42,5);

plot3(0,0,0,'.','Color',[.5,.5,.5,.5],'MarkerSize',40);
for cycle_classes = 1:numel(Odorants)
    for cycle_points = 0:line_mult/new_bin_size:size(scores_norm,2)
        if cycle_points ~= 0
            plot3([0,scores_norm(1,cycle_points,cycle_classes)], [0,scores_norm(2,cycle_points,cycle_classes)],[0,scores_norm(3,cycle_points,cycle_classes)],...
                'Color',[Colors(cycle_classes,:),0.2],'Linewidth', .5);
            
            switch cycle_points % works for 50 msec time bins, need to update for compatability with other bin sizes
                case bps*.25, text(scores_norm(1,cycle_points,cycle_classes), scores_norm(2,cycle_points,cycle_classes), scores_norm(3,cycle_points,cycle_classes),...
                        num2str(time(1) + .25),...
                        'FontSize',14,'Color',Colors(cycle_classes,:)) %,'BackgroundColor',[Colors(cycle_classes,:),0.25])                
                case bps*.5, text(scores_norm(1,cycle_points,cycle_classes), scores_norm(2,cycle_points,cycle_classes), scores_norm(3,cycle_points,cycle_classes),...
                        num2str(time(1) + .5),...
                        'FontSize',14,'Color',Colors(cycle_classes,:)) %,'BackgroundColor',[Colors(cycle_classes,:),0.5])   
                case bps*1, text(scores_norm(1,cycle_points,cycle_classes), scores_norm(2,cycle_points,cycle_classes), scores_norm(3,cycle_points,cycle_classes),...
                        num2str(time(1) + 1),...
                        'FontSize',14,'Color',Colors(cycle_classes,:)) %,'BackgroundColor',[Colors(cycle_classes,:),0.25])                
                case bps*1.5, text(scores_norm(1,cycle_points,cycle_classes), scores_norm(2,cycle_points,cycle_classes), scores_norm(3,cycle_points,cycle_classes),...
                        num2str(time(1) + 1.5),...
                        'FontSize',14,'Color',Colors(cycle_classes,:)) %,'BackgroundColor',[Colors(cycle_classes,:),0.5])   
                case bps*2, text(scores_norm(1,cycle_points,cycle_classes), scores_norm(2,cycle_points,cycle_classes), scores_norm(3,cycle_points,cycle_classes),...
                        num2str(time(1) + 2),...
                        'FontSize',14,'Color',Colors(cycle_classes,:)) %,'BackgroundColor',[Colors(cycle_classes,:),0.25])                
                case bps*4, text(scores_norm(1,cycle_points,cycle_classes), scores_norm(2,cycle_points,cycle_classes), scores_norm(3,cycle_points,cycle_classes),...
                        num2str(time(1) + 4),...
                        'FontSize',14,'Color',Colors(cycle_classes,:)) %,'BackgroundColor',[Colors(cycle_classes,:),0.5])   
            end
        end
    end
    PCA_handle(cycle_classes) = plot3(scores_norm(1,1:size(scores_norm,2),cycle_classes),scores_norm(2,1:size(scores_norm,2),cycle_classes),scores_norm(3,1:size(scores_norm,2),cycle_classes),'Color',[Colors(cycle_classes,:),0.85],'Linewidth', 2.5);
end
title(static_axes,PCA_title,'FontSize',20,'FontWeight','Bold');
lims = [min(scores_norm(1,:,:),[],'all') max(scores_norm(1,:,:),[],'all'); min(scores_norm(2,:,:),[],'all') max(scores_norm(2,:,:),[],'all'); min(scores_norm(3,:,:),[],'all') max(scores_norm(3,:,:),[],'all')]; 

xlim(static_axes,lims(1,:));
ylim(static_axes,lims(2,:));
zlim(static_axes,lims(3,:));
xlabel(static_axes,['PC1 (' num2str(round(var_exp(1)*100)/100) '%)'],'FontName','Arial','FontSize',18,'FontWeight','Bold');
ylabel(static_axes,['PC2 (' num2str(round(var_exp(2)*100)/100) '%)'],'FontName','Arial','FontSize',18,'FontWeight','Bold');
zlabel(static_axes,['PC3 (' num2str(round(var_exp(3)*100)/100) '%)'],'FontName','Arial','FontSize',18,'FontWeight','Bold');
annotation(static_fig,'textbox', [0.2, 0.9, 0, 0], 'String', ['n = ' num2str(size(data,1))],'Units','normalized','Color','k','FontSize',14,'FontWeight','Bold','FitBoxToText','on','HorizontalAlignment','center');
% legend(PCA_handle,Odorants,'location','eastoutside','FontSize',16,'FontWeight','Bold');
return