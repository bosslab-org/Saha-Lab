function static_fig = PSTH(data,Colors,Odorants,time,bins_per_sec,PSTH_title)
files = size(data,1);
data = squeeze(mean(data,[1,2]));

static_fig = figure('Units','normalized','Position',[0 0.5 0.5 0.39]);  % generate single figure for all Odorant PCA trajectories
static_axes = axes(static_fig);
hold on; grid on; box off; %view(-66,57);

for cycle_classes = 1:size(data,2)
%     plot3(repmat(cycle_classes,1,size(data,1)), 1:size(data,1), flip(data(:,cycle_classes)),'-','Color',[Colors(cycle_classes,:),0.8],'LineWidth',3);
%     subplot(size(data,2),1,cycle_classes)
    plot(1:size(data,1),data(:,cycle_classes),'-','Color',[Colors(cycle_classes,:),0.3],'LineWidth',3)
%     ylim([-.5 .5]);
%     xlim([0 40])
end


% x_patch = [0 0 numel(Odorants) numel(Odorants)];
% y_patch = [-((time(1)*bins_per_sec)+1) -((time(1)-4)*bins_per_sec) -((time(1)-4)*bins_per_sec) -((time(1)*bins_per_sec)+1)];
% z_patch = [min(data,[],'all') min(data,[],'all') min(data,[],'all') min(data,[],'all')];
% stim_patch = patch(x_patch,size(data,1)-y_patch,z_patch,'k', 'FaceAlpha', 0.05, 'EdgeColor','k', 'EdgeAlpha', 0.05);
% 
% 
% title(static_axes,PSTH_title,'FontSize',20,'FontWeight','Bold');
% xticklabels(static_axes,[0 Odorants]);
% yticks(static_axes,0:size(data,1)/(time(2)-time(1)):size(data,1))
% yticklabels(static_axes,fliplr(time(1):time(2)));
% 
% xlim(static_axes,[1 numel(Odorants)]);
% ylim(static_axes,[0 size(data,1)]);
% zlim(static_axes,[min(data,[],'all') max(data,[],'all')]);
% 
% xlabel(static_axes,['Stimulus'],'FontName','Arial','FontSize',18,'FontWeight','Bold');
% ylabel(static_axes,['Time'],'FontName','Arial','FontSize',18,'FontWeight','Bold');
% zlabel(static_axes,['Voltage (uV)'],'FontName','Arial','FontSize',18,'FontWeight','Bold');
% annotation(static_fig,'textbox', [0.2, 0.9, 0, 0], 'String', ['n = ' num2str(files)],'Units','normalized','Color','k','FontSize',14,'FontWeight','Bold','FitBoxToText','on','HorizontalAlignment','center');
return