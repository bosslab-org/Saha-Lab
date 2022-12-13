function pairwise_fig = pairwise_distance(data)

pairwise_fig = figure('Units','normalized','Position',[0 0 0.5 0.39]); %LDA_axes = axes(LDA_fig);
hold on; grid on;
data_average = [];
for cycle_stimuli = 1:size(data,4)-1
    data_base = mean(data(:,:,:,cycle_stimuli),[1 3]);
    data_test = mean(data(:,:,:,cycle_stimuli+1:size(data,4)),[1 3]);
    data_temp = abs(data_base-data_test);
    for cycle_plot_stim = 1:size(data_temp,4)
        plot(1:size(data,2),squeeze(data_temp(:,:,:,cycle_plot_stim)),'Color',[0,0,0,0.25])
    end
    data_average = cat(4,data_average,data_temp);
end

plot(1:size(data,2),mean(data_average,4),'Color',[0,0,0,0.75])

xticks(0:size(data,2)/8:size(data,2))
xticklabels(0:0.5:4)

return