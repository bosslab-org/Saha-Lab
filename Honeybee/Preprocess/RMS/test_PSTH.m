function PSTH_fig = test_PSTH(data, exp_pars, channel)

PSTH_fig = figure('Units','normalized','Position',[0 0 0.75 0.75]);
for cycle_trials = 1:size(data,3)
    subplot(size(data,3),1,cycle_trials)
    plot(1:size(data,4),squeeze(data(:,channel,cycle_trials,:)));
    xticks([0:exp_pars(3)*5:size(data,4)])
    xticklabels([0:5:size(data,4)/exp_pars(3)])
end