function Neuron_ANOVA(readpath,pvalue,base_stim_idx,test_stim_idx,bin_size)

test_duration = 4;

stim_on = 2;
bins_per_sec = 1000/bin_size;
var_ext = '_spikes_bin';

load(readpath);
stimuli = who(readpath, ['*' var_ext]);


for cycle_odors = 1:length(stimuli)
    eval(sprintf('temp = %s;',stimuli{cycle_odors}));
    AllPN_test(:,:,cycle_odors) = sum(temp(:,:,stim_on*bins_per_sec+1:(stim_on+test_duration)*bins_per_sec),3);
    
    % Average over trials and pre-stim samples
    AllPN_base(:,:,cycle_odors) = mean(temp(:,:,(stim_on-2)*bins_per_sec+1:stim_on*bins_per_sec),[2 3]);   
end

% PN_Diff_test = AllPN_test;
PN_Diff_test = AllPN_test-AllPN_base;

for cycle_figs = 1:size(PN_Diff_test,3)-1
    figure; hold on; grid on;
    for cycle_neurons = 1:size(PN_Diff_test,1)
        if cycle_figs == 1
            data(:,:,cycle_neurons) = permute(PN_Diff_test(cycle_neurons,:,:),[2 3 1]);
            [mu,~,stats] = anova1(data(:,:,cycle_neurons),stimuli,'off');
            stats_comp(:,:,cycle_neurons) = multcompare(stats,'alpha',pvalue,'display','off','ctype','bonferroni');
            base_stim_data(cycle_neurons) = mean(data(:,base_stim_idx,cycle_neurons));
            test_stim_data(:,cycle_neurons) = mean(data(:,test_stim_idx,cycle_neurons));
            base_SEM(cycle_neurons) = std(data(:,base_stim_idx,cycle_neurons))/sqrt(5);
            test_SEM(:,cycle_neurons) = std(data(:,test_stim_idx,cycle_neurons))/sqrt(5);
        end
        h = plot(base_stim_data(cycle_neurons),test_stim_data(cycle_figs,cycle_neurons),'o','MarkerSize',3);
        l1 = line([base_stim_data(cycle_neurons),base_stim_data(cycle_neurons)],...
            [test_stim_data(cycle_figs,cycle_neurons)-test_SEM(cycle_figs,cycle_neurons),test_stim_data(cycle_figs,cycle_neurons)+test_SEM(cycle_figs,cycle_neurons)]);
        l2 = line([base_stim_data(cycle_neurons)-base_SEM(cycle_neurons),base_stim_data(cycle_neurons)+base_SEM(cycle_neurons)],...
            [test_stim_data(cycle_figs,cycle_neurons),test_stim_data(cycle_figs,cycle_neurons)]);
        if stats_comp(cycle_figs,3,cycle_neurons)>0    % spike count decreases significantly, shown in blue
            h.Color = [0 0 1];
            l1.Color = [0 0 1];
            l2.Color = [0 0 1];
        elseif  stats_comp(cycle_figs,5,cycle_neurons)<0   % spike count increases significantly, shown in red
            h.Color = [1 0 0];
            l1.Color = [1 0 0];
            l2.Color = [1 0 0];
        else    % no significant changes, shown in grey
            h.Color = [0.5 0.5 0.5];
            l1.Color = [0.6 0.6 0.6];
            l2.Color = [0.6 0.6 0.6];
        end
    end
    xlabel(stimuli{base_stim_idx});
    ylabel(stimuli{test_stim_idx(cycle_figs)});
    xlim([0 max(base_stim_data,[],'all')+max(base_SEM,[],'all')]);
    ylim([0 max(test_stim_data,[],'all')+max(test_SEM,[],'all')]);
end