function Conf_fig = Confusion(data,class,Confusion_title)


Conf_fig = figure('Units','normalized','Position', [0.5 0.5 0.5 0.5]);

imagesc(data);

xticks(1:numel(class));
yticks(1:numel(class));
xticklabels(class);
yticklabels(class);
xtickangle(45);
ytickangle(45);
ax=gca;
ax.FontSize = 16;
xlabel('Target Odor', 'FontSize', 18)
ylabel('Predicted Odor', 'FontSize', 18)
set(text(1.1, .42, 'Accuracy', 'Units', 'normalized', 'FontSize', 16), 'Rotation', 90)
caxis([0,100]);
colormap([linspace(1,0.5,20); linspace(1,0,20); linspace(1,1,20)]')
% text(.83,.95,['{\it n} = ' num2str(size(CombinedTrainData,1))],'Units','normalized', 'Color', 'black', 'FontSize', 18)

total_acc = mean(diag(data));
title([Confusion_title '- Accuracy = ' num2str(total_acc) '%'])
colorbar('Ticks', [0, 100], 'TickLabels', {'0%', '100%'}, 'FontSize', 16);
return