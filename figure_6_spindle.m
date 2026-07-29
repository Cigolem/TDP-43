
%% Time Frequency Plots — PFC ECoG Light ON — WT vs HOM - Figure 6 a and b

% Wavelet_CG_relative returns relative power in linear (%) and log10 form.
% TF arrays store the linear output; log conversion (10*log10) is applied here at plotting.

clear


load('D:\TDP\submission_sleep\final_spindle_scripts\sp_swa_assoc_results_5_withTF_PFC_EcOG.mat')

clearvars -except all_results

channel = 'PFC_EcOG';
light   = 'light_on';

tf_types  = {'tf_sol_swa_norm', 'tf_sol_sp_norm', 'tf_coupled_norm'};
tf_titles = {'SWA', 'Spindle', 'SWA-Spindle'};

% Get indices
groups  = {all_results.(channel).(light).group};
WT_idx  = find(strcmp(groups, 'WT'));
HOM_idx = find(strcmp(groups, 'Hom'));

% Get axes
time_axis = all_results.(channel).(light)(1).time_axis;
freqs     = all_results.(channel).(light)(1).freqs;
n_freqs   = length(freqs);
n_time    = length(time_axis);

n_WT  = length(WT_idx);
n_HOM = length(HOM_idx);

figure('Position', [100 100 900 600]);

for tf = 1:length(tf_types)
    tf_field = tf_types{tf};

    wt_mats  = zeros(n_freqs, n_time, n_WT);
    hom_mats = zeros(n_freqs, n_time, n_HOM);

    for i = 1:n_WT
        wt_mats(:,:,i)  = all_results.(channel).(light)(WT_idx(i)).(tf_field);
    end
    for i = 1:n_HOM
        hom_mats(:,:,i) = all_results.(channel).(light)(HOM_idx(i)).(tf_field);
    end

    wt_mean  = mean(wt_mats,  3);
    hom_mean = mean(hom_mats, 3);

    if tf == 1
        clim_range = [0.8 2.5];
    elseif tf == 2
        clim_range = [0.8 6.5];
    else 
        clim_range = [0.8 5.5];
    end
    

    if tf == 1
        ylim_range = [1 10];
    else
        ylim_range = [1 25];
    end


    % WT
    subplot(3, 2, (tf-1)*2 + 1);
    contourf(time_axis, freqs, 10*log10(wt_mean), 40, 'LineColor','none');
    colormap jet; colorbar;
    hold on;
    plot([0 0], [freqs(1) freqs(end)], 'k--', 'LineWidth', 2);
    if tf == 3, xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); end
    title(sprintf('WT — %s', tf_titles{tf}), 'FontSize', 13);
    
    xlim([-2 2]); ylim(ylim_range);
    clim(clim_range);
    set(gca, 'FontSize', 18, 'TickDir', 'out', 'box', 'off');

    % HOM
    subplot(3, 2, (tf-1)*2 + 2);
    contourf(time_axis, freqs, 10*log10(hom_mean), 40, 'LineColor','none');
    colormap jet; colorbar;
    hold on;
    plot([0 0], [freqs(1) freqs(end)], 'k--', 'LineWidth', 2);
    if tf == 3, xlabel('Time (s)', 'FontSize', 18, 'FontWeight', 'bold'); end
    title(sprintf('Q331K — %s', tf_titles{tf}), 'FontSize', 20, 'FontWeight', 'bold');
    xlim([-2 2]); ylim(ylim_range);
    clim(clim_range);
    set(gca, 'FontSize', 18, 'TickDir', 'out', 'box', 'off');

    % Shared y-label per column
    han1 = axes(gcf, 'Visible', 'off'); han1.Position = [0.09 0.11 0.01 0.815];
    han1.YLabel.Visible = 'on';
    ylabel(han1, 'Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');

    han2 = axes(gcf, 'Visible', 'off'); han2.Position = [0.53 0.11 0.01 0.815];
    han2.YLabel.Visible = 'on';
    ylabel(han2, 'Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');

end

% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_6_a_b.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_6_a_b.tiff'), 'Resolution', 300);



%% Power violin plots — PFC ECoG Light ON - Figure 6 c-f

clearvars -except all_results

color_wt  = [0.25 0.65 0.55];
color_hom = [0.9 0.2 0.7];

% Extract data
WT_sol_slow   = []; HOM_sol_slow   = [];
WT_coup_full  = []; HOM_coup_full  = [];
WT_coup_slow  = []; HOM_coup_slow  = [];
WT_coup_fast  = []; HOM_coup_fast  = [];

for i = 1:length(all_results.PFC_EcOG.light_on)
    entry = all_results.PFC_EcOG.light_on(i);
    if strcmp(entry.group, 'WT')
        WT_sol_slow(end+1)  = entry.power_solitary_slow;
        WT_coup_full(end+1) = entry.power_coupled_full;
        WT_coup_slow(end+1) = entry.power_coupled_slow;
        WT_coup_fast(end+1) = entry.power_coupled_fast;
    else
        HOM_sol_slow(end+1)  = entry.power_solitary_slow;
        HOM_coup_full(end+1) = entry.power_coupled_full;
        HOM_coup_slow(end+1) = entry.power_coupled_slow;
        HOM_coup_fast(end+1) = entry.power_coupled_fast;
    end
end

% Stats
[~, p1, ~, st1] = ttest2(WT_sol_slow,  HOM_sol_slow,  'Vartype','unequal');
[~, p2, ~, st2] = ttest2(WT_coup_full, HOM_coup_full, 'Vartype','unequal');
[~, p3, ~, st3] = ttest2(WT_coup_slow, HOM_coup_slow, 'Vartype','unequal');
[~, p4, ~, st4] = ttest2(WT_coup_fast, HOM_coup_fast, 'Vartype','unequal');

fprintf('power_solitary_slow:  t(%.1f)=%.3f, p=%.4f\n', st1.df, st1.tstat, p1);
fprintf('power_coupled_full:   t(%.1f)=%.3f, p=%.4f\n', st2.df, st2.tstat, p2);
fprintf('power_coupled_slow:   t(%.1f)=%.3f, p=%.4f\n', st3.df, st3.tstat, p3);
fprintf('power_coupled_fast:   t(%.1f)=%.3f, p=%.4f\n', st4.df, st4.tstat, p4);

wt_data_all  = {WT_sol_slow,  WT_coup_full,  WT_coup_slow,  WT_coup_fast};
hom_data_all = {HOM_sol_slow, HOM_coup_full, HOM_coup_slow, HOM_coup_fast};
titles       = {'Solitary - Slow band', 'Coupled - Full band', ...
                'Coupled - Slow band',  'Coupled - Fast band'};
pvals        = [p1, p2, p3, p4];

pos_wt  = 1;
pos_hom = 2;

figure('Position', [100 100 800 700]);

for panel = 1:4
    subplot(2, 2, panel);
    hold on;

    wt_data  = wt_data_all{panel};
    hom_data = hom_data_all{panel};

    simple_violin({wt_data},  pos_wt,  color_wt,  0.35);
    simple_violin({hom_data}, pos_hom, color_hom, 0.35);

    scatter(pos_wt  * ones(1,length(wt_data)),  wt_data,  80, 'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor',color_wt, ...
        'LineWidth',1.5,'MarkerFaceAlpha',0.8);
    scatter(pos_hom * ones(1,length(hom_data)), hom_data, 80, 'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor',color_hom, ...
        'LineWidth',1.5,'MarkerFaceAlpha',0.8);

    plot([pos_wt-0.2  pos_wt+0.2],  [mean(wt_data)  mean(wt_data)],  'k-','LineWidth',2);
    plot([pos_hom-0.2 pos_hom+0.2], [mean(hom_data) mean(hom_data)], 'k-','LineWidth',2);

    p = pvals(panel);
    pstr = sprintf('p = %.3f', p);
    text(1.6, 11, pstr, 'HorizontalAlignment','center','FontSize',16);
    ylim([2 13]);
    xlim([0.4 2.6]);
    xticks([1 2]);
    xticklabels({'WT', 'Q331K'});
    if panel == 1 || panel == 2
        xticklabels({});
    end

    if panel == 1 || panel == 3
        ylabel('Relative Power (%)', 'FontSize', 20, 'FontWeight', 'bold');
    end
    title(titles{panel}, 'FontSize', 20, 'FontWeight', 'bold');
    box off;
    set(gca, 'FontSize', 20, 'LineWidth', 1.2, 'FontWeight', 'bold');
end


% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_6_c_f.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_6_c_f.tiff'), 'Resolution', 300);

%% end

