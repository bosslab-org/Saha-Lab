function Conf_fig = Confusion(data,class,Confusion_title,res,acc)

Conf_fig = figure('Units','normalized','Position', [0.5 0.5 0.5 0.5]);
Conf_axes = axes(Conf_fig);
imagesc(data);

if acc == 1
    for cycle_rows = 1:size(data,1)
        for cycle_cols = 1:size(data,2)
            if data(cycle_rows,cycle_cols) >= 100/res
                text(cycle_cols,cycle_rows,[num2str(round(data(cycle_rows,cycle_cols),1)) '%'],...
                    'HorizontalAlignment', 'center',...
                    'Color','k',...
                    'FontSize',12,...
                    'FontWeight','bold');
            end
        end
    end
end

xticks(1:numel(class));
yticks(1:numel(class));
xticklabels(class);
yticklabels(class);
xtickangle(45);
ytickangle(45);
Conf_axes.FontSize = 16;
xlabel('Target Odor', 'FontSize', 18)
ylabel('Predicted Odor', 'FontSize', 18)
set(text(1.1, .42, 'Accuracy', 'Units', 'normalized', 'FontSize', 16), 'Rotation', 90)
caxis([0,100]);
colormap([linspace(1,0.5,res); linspace(1,0,res); linspace(1,1,res)]')

total_acc = mean(diag(data));
title([Confusion_title '- Accuracy = ' num2str(round(total_acc,2)) '%'])
colorbar('Ticks', [0, 100], 'TickLabels', {'0%', '100%'}, 'FontSize', 16);
return
