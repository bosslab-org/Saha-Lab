function delta_rms_PSVT(data, Colors, Odorants)

figure('Units','normalized','Position',[0 0 0.75 1]); hold on;
for cycle_stimuli = 1:size(data,4)
    plot(1:size(data,2),mean(data(:,:,:,cycle_stimuli),[1 3]),'Color',Colors(cycle_stimuli,:))
end
legend(Odorants,'location','eastoutside','FontSize',16,'FontWeight','Bold');

return