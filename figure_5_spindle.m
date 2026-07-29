


%% three panel figure - Figure 5

load('D:\TDP\submission_sleep\final_spindle_scripts\sp_swa_assoc_results_5_on_and_off.mat')
load('D:\TDP\submission_sleep\final_spindle_scripts\summary_ver5_on_off.mat')

% Panel 1: Spindle density
% Panel 2: Solitary spindle duration
% Panel 3: Coupled spindle duration

clearvars -except all_results summary_table_ver5_on_off

color_wt  = [0.25 0.65 0.55];
color_hom = [0.9 0.2 0.7];

% Extract Light ON data
WT_density   = [];
HOM_density  = [];
WT_sol_on    = [];
HOM_sol_on   = [];
WT_coup_on   = [];
HOM_coup_on  = [];

for i = 1:length(all_results.PFC_EcOG.light_on)
    entry = all_results.PFC_EcOG.light_on(i);

    if strcmp(entry.group, 'WT')
        WT_density(end+1)  = entry.spindle_density;
        WT_sol_on(end+1)   = entry.dur_solitary;
        WT_coup_on(end+1)  = entry.dur_coupled;
    else
        HOM_density(end+1) = entry.spindle_density;
        HOM_sol_on(end+1)  = entry.dur_solitary;
        HOM_coup_on(end+1) = entry.dur_coupled;
    end
end

WT_count  = length(WT_density);
HOM_count = length(HOM_density);

% Stats
[~, p_density,  ~, st_density]  = ttest2(WT_density,  HOM_density,  'Vartype', 'unequal');
[~, p_sol_on,   ~, st_sol_on]   = ttest2(WT_sol_on,   HOM_sol_on,   'Vartype', 'unequal');
[~, p_coup_on,  ~, st_coup_on]  = ttest2(WT_coup_on,  HOM_coup_on,  'Vartype', 'unequal');

fprintf('Spindle density:        t(%.1f)=%.3f, p=%.4f\n', st_density.df,  st_density.tstat,  p_density);
fprintf('Solitary duration ON:   t(%.1f)=%.3f, p=%.4f\n', st_sol_on.df,   st_sol_on.tstat,   p_sol_on);
fprintf('Coupled duration ON:    t(%.1f)=%.3f, p=%.4f\n', st_coup_on.df,  st_coup_on.tstat,  p_coup_on);

% Organize data by panel
wt_data_all  = {WT_density, WT_sol_on, WT_coup_on};
hom_data_all = {HOM_density, HOM_sol_on, HOM_coup_on};
pvals        = [p_density, p_sol_on, p_coup_on];
ylabels      = {'Spindles/min', 'Duration (s)', 'Duration (s)'};
titles       = {'Spindle Density', 'Solitary Spindle Duration', 'Coupled Spindle Duration'};

pos_wt  = 1;
pos_hom = 2;

figure('Position', [100 100 1300 500]);

for panel = 1:3
    subplot(1,3,panel);
    hold on;

    wt_data  = wt_data_all{panel};
    hom_data = hom_data_all{panel};

    simple_violin({wt_data},  pos_wt,  color_wt,  0.35);
    simple_violin({hom_data}, pos_hom, color_hom, 0.35);

    scatter(pos_wt  * ones(1,length(wt_data)),  wt_data, 100, 'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor',color_wt, ...
        'LineWidth',1.5,'MarkerFaceAlpha',0.8);

    scatter(pos_hom * ones(1,length(hom_data)), hom_data, 100, 'o', ...
        'MarkerEdgeColor','k','MarkerFaceColor',color_hom, ...
        'LineWidth',1.5,'MarkerFaceAlpha',0.8);

    % Mean lines
    plot([pos_wt-0.2  pos_wt+0.2],  [mean(wt_data)  mean(wt_data)],  'k-', 'LineWidth', 2);
    plot([pos_hom-0.2 pos_hom+0.2], [mean(hom_data) mean(hom_data)], 'k-', 'LineWidth', 2);

    % P value text
    p = pvals(panel);
    if p < 0.001
        pstr = 'p < 0.001';
    else
        pstr = sprintf('p = %.3f', p);
    end

    if panel == 1
        ylim([1.2 4]);
        text_y = 3.8;
    else
        ylim([0.6 1.2]);
        text_y = 1.15;
    end

    text(mean([pos_wt pos_hom]), text_y, pstr, ...
        'HorizontalAlignment', 'center', 'FontSize', 15);

    xlim([0.4 2.6]);
    xticks([pos_wt pos_hom]);
    xticklabels({'WT', 'Q331K'});
    ylabel(ylabels{panel}, 'FontSize', 18, 'FontWeight', 'bold');
    title(titles{panel}, 'FontSize', 20, 'FontWeight', 'bold');
    text(-0.15, 1.05, char('a'+panel-1), 'Units','normalized','FontWeight','bold','FontSize',20)
    box off;
    set(gca, 'FontSize', 17, 'LineWidth', 1.5);
end

% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_5.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_5.tiff'), 'Resolution', 300);

