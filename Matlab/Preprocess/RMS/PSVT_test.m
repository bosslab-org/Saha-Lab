function PSVT_test(data_filt,exp_pars)

for cycle_tetrodes = 1:size(data_filt,1)
    figure;
    for cycle_trials = 1:size(data_filt,3)
        subplot(size(data_filt,3),1,cycle_trials);
        plot(1:size(data_filt,4),squeeze(data_filt(cycle_tetrodes,1,cycle_trials,:)));
        xticks(0:exp_pars(5):exp_pars(3)*exp_pars(5))
        xticklabels(0:exp_pars(3))
    end
end
