function Confusion(data,Class,Confusion_title)

% figure('Units','normalized','Position',[0.5 0.5 0.5 0.39]);
% conf_chart = confusionchart(y,yhat);

figure('Position', [50 50 800 700]);
imagesc(data');

xticks([1:numel(Class)]);
yticks([1:numel(Class)]);
xticklabels(Class);
yticklabels(Class);
xtickangle(45);
ytickangle(45);
ax=gca;
ax.FontSize = 16;
xlabel('Target Odor', 'FontSize', 18)
ylabel('Predicted Odor', 'FontSize', 18)
set(text(1.1, .42, 'Accuracy', 'Units', 'normalized', 'FontSize', 16), 'Rotation', 90)
caxis([0,1]);
colormap([linspace(1,0.5,20); linspace(1,0,20); linspace(1,1,20)]')
% text(.83,.95,['{\it n} = ' num2str(size(CombinedTrainData,1))],'Units','normalized', 'Color', 'black', 'FontSize', 18)
meanAccuracy = mean(diag(data))*100;
title([Confusion_title '- Accuracy = ' num2str(meanAccuracy) '%'])
colorbar('Ticks', [0, 1], 'TickLabels', {'0%', '100%'}, 'FontSize', 16);