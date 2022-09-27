function [data, SS_PSVT_fig] = process_artifacts(data, tref, thresh, stim_on, stim_off, sample_rate, gen_plot)
[pos_pks,pos_locs] = findpeaks(data,1:size(data,1),'MinPeakDistance',tref);
pos_artifacts = pos_locs(pos_pks > thresh*std(data));
[neg_pks,neg_locs] = findpeaks(-data,1:size(data,1),'MinPeakDistance',tref);
neg_artifacts = neg_locs(neg_pks > thresh*std(data));
artifacts = cat(2, pos_artifacts, neg_artifacts);

if gen_plot == 1
    ylim_min = min(data);
    ylim_max = max(data);

    SS_PSVT_fig = figure('Position', [0 0 1000 750]);
    subplot(2,1,1)
    plot(1:size(data,1), data,'-*','MarkerIndices', artifacts,'MarkerSize', 2, 'MarkerFaceColor', 'red','MarkerEdgeColor','red')

    yline(thresh*std(data),'--','Color', 'k','LineWidth',.5);
    yline(-thresh*std(data),'--','Color', 'k','LineWidth',.5);
    patch_x = [stim_on stim_on stim_off stim_off];
    patch_y = [min(data) max(data) max(data) min(data)]; 
    patch(patch_x,patch_y,'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.0);
    xticks(0:sample_rate:size(data,1))
    xticklabels(0:(size(data,1)/sample_rate))
    % xlim([21.45*sample_rate 21.65*sample_rate])

    xlim([0 size(data,1)])
    ylim([ylim_min ylim_max])
end

data_mean = mean(data);
for ii = 1:numel(artifacts)
    data(artifacts(ii)-19:artifacts(ii)+20) = data_mean;
end

if gen_plot == 1
    subplot(2,1,2)
    plot(1:size(data,1), data,'-*','MarkerIndices', artifacts,'MarkerSize', 2, 'MarkerFaceColor', 'red','MarkerEdgeColor','red')

    yline(thresh*std(data),'--','Color', 'k','LineWidth',.5);
    yline(-thresh*std(data),'--','Color', 'k','LineWidth',.5);
    patch_x = [stim_on stim_on stim_off stim_off];
    patch_y = [ylim_min ylim_max ylim_max ylim_min]; 
    patch(patch_x,patch_y,'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.0);
    xticks(0:sample_rate:size(data,1))
    xticklabels(0:(size(data,1)/sample_rate))
    % xlim([21.45*sample_rate 21.65*sample_rate])

    xlim([0 size(data,1)])
    ylim([ylim_min ylim_max])
end
    
return