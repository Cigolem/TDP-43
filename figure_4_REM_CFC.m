

 %% Compute Theta-Gamma PAC — REM — All Channels — Early/Mid/Late

clear all; close all;

fs = 1000;
file_path = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM\EEG_arrays';
output_path = 'D:\TDP\submission_sleep\PAC_scripts\final'; 
files = dir(fullfile(file_path, 'EEG_*.mat'));

% files = files(9);  % test on first file only


channels   = {'PFC_LFP', 'PFC_EcOG', 'V1_EcOG'};
phas_freqs = 3:0.5:10;
ampl_freqs = 15:5:100;

nrem_samples = 30 * fs;   % first 30s is NREM
epoch_len    = 10 * fs;   % 10s epochs
epoch_select = [1:6]; % all epochs

PAC_WT  = struct();
PAC_HOM = struct();
for ch = 1:length(channels)
    PAC_WT.(channels{ch})  = {};
    PAC_HOM.(channels{ch}) = {};
end

wt_count = 0; hom_count = 0;

for f = 1:length(files)
    load(fullfile(file_path, files(f).name));
    fprintf('Processing %s\n', files(f).name);
    
    if contains(files(f).name, 'WT')
        group = 'WT'; wt_count = wt_count + 1; mouse_idx = wt_count;
    else
        group = 'HOM'; hom_count = hom_count + 1; mouse_idx = hom_count;
    end
    
    for ch = 1:length(channels)
        chan = channels{ch};
        rem_bouts = EEG.(chan).NREM_REM;
        n_bouts   = length(rem_bouts);
        
        pac_bouts = [];
        valid_count = 0;
        
        for b = 1:n_bouts
            bout = rem_bouts{b};
            
            % Extract REM portion only (skip first 30s NREM)
            if length(bout) < nrem_samples + 6*epoch_len
                continue;
            end
            rem_data = bout(nrem_samples+1 : nrem_samples + 6*epoch_len);
            
            valid_count = valid_count + 1;
            
            % Compute PAC on epochs 1, 4, 6
            for t = 1:length(epoch_select)
                ep = epoch_select(t);
                start_idx = (ep-1)*epoch_len + 1;
                end_idx   = ep*epoch_len;
                segment   = rem_data(start_idx:end_idx);
                pac_bouts(:,:,t,valid_count) = CFC_TDP_core(segment, phas_freqs, ampl_freqs);
            end
        end
        
        fprintf('  %s: %d valid bouts\n', chan, valid_count);
        
        if valid_count > 0
            mouse_pac = squeeze(nanmean(pac_bouts, 4));
        else
            mouse_pac = NaN(length(phas_freqs), length(ampl_freqs), length(epoch_select));
        end
        
        if strcmp(group, 'WT')
            PAC_WT.(chan){mouse_idx}  = mouse_pac;
        else
            PAC_HOM.(chan){mouse_idx} = mouse_pac;
        end
    end
end

fprintf('\nWT mice: %d | HOM mice: %d\n', wt_count, hom_count);
save(fullfile(output_path, 'PAC_REM_epochs_six_epochs.mat'), 'PAC_WT', 'PAC_HOM', 'phas_freqs', 'ampl_freqs');


%% Plot V1 ECoG - Early (ep1-2) and Late (ep5-6) averaged (figure 4 a-d)

load('D:\TDP\submission_sleep\PAC_scripts\final\PAC_REM_epochs_six_epochs.mat')
clearvars -except PAC_HOM PAC_WT phas_freqs ampl_freqs

phas_freqs = 3:0.5:10;
ampl_freqs = 15:5:100;

chan       = 'V1_EcOG';
chan_label = 'V1 ECoG';
clim_range = [0 2];

n_wt  = length(PAC_WT.(chan));
n_hom = length(PAC_HOM.(chan));

wt_sum  = zeros(length(phas_freqs), length(ampl_freqs), 6);
hom_sum = zeros(length(phas_freqs), length(ampl_freqs), 6);

for m = 1:n_wt
    wt_sum = wt_sum + PAC_WT.(chan){m};
end
for m = 1:n_hom
    hom_sum = hom_sum + PAC_HOM.(chan){m};
end

wt_mean  = wt_sum  / n_wt;
hom_mean = hom_sum / n_hom;

% Average first 2 epochs (epoch 1 and 2)
wt_early  = mean(wt_mean(:,:,1:2),  3);
hom_early = mean(hom_mean(:,:,1:2), 3);

% Average last 2 epochs (epoch 5 and 6)
wt_late  = mean(wt_mean(:,:,5:6),  3);
hom_late = mean(hom_mean(:,:,5:6), 3);


% Plot - 1x4 (Early WT/Q331K)

figure('Position', [100 100 1500 600]);

subplot(1,4,1)
contourf(phas_freqs, ampl_freqs, wt_early', 50, 'LineColor', 'none');
clim(clim_range); set(gca, 'FontSize', 18);
ylabel('Amplitude (Hz)', 'FontSize', 20, 'FontWeight', 'bold');
% text(-0.1, 1.05, 'a', 'Units','normalized','FontWeight','bold','FontSize',20)
title('WT - Early REM', 'FontSize', 20, 'FontWeight', 'bold');


subplot(1,4,2)
contourf(phas_freqs, ampl_freqs, hom_early', 50, 'LineColor', 'none');
clim(clim_range); set(gca, 'FontSize', 18);
% text(-0.1, 1.05, 'b', 'Units','normalized','FontWeight','bold','FontSize',20)
title('Q331K - Early REM', 'FontSize', 20, 'FontWeight', 'bold');


subplot(1,4,3)
contourf(phas_freqs, ampl_freqs, wt_late', 50, 'LineColor', 'none');
clim(clim_range); set(gca, 'FontSize', 18);
% text(-0.1, 1.05, 'c', 'Units','normalized','FontWeight','bold','FontSize',20)
title('WT - Late REM', 'FontSize', 20, 'FontWeight', 'bold');


subplot(1,4,4)
contourf(phas_freqs, ampl_freqs, hom_late', 50, 'LineColor', 'none');
clim(clim_range); set(gca, 'FontSize', 18);
% text(-0.1, 1.05, 'd', 'Units','normalized','FontWeight','bold','FontSize',20)
title('Q331K - Late REM', 'FontSize', 20, 'FontWeight', 'bold');


% Shared x-label between the subplots
    han = axes(gcf, 'Visible', 'off');
    han.XLabel.Visible = 'on';
    xlabel(han, 'Phase (Hz)', 'FontSize', 20, 'FontWeight', 'bold');

% Shared colourbar
colormap jet; 
clim(han, clim_range);
h = colorbar;
h.FontSize = 16;
h.Position = [0.92 0.15 0.015 0.7];   % [left bottom width height]
ylabel(h, 'Modulation Index', 'FontSize', 20, 'FontWeight', 'bold');


% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_4_a_d.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_4_a_d.tiff'), 'Resolution', 300);

%% Stats t-tests and Violin Plots of Early and Late REM between groups (figure 4 e-f)

clearvars -except PAC_HOM PAC_WT phas_freqs ampl_freqs

chan       = 'V1_EcOG';
n_wt  = length(PAC_WT.(chan));
n_hom = length(PAC_HOM.(chan));

theta_idx = phas_freqs >= 6 & phas_freqs <= 8;
gamma_idx = ampl_freqs >= 65 & ampl_freqs <= 70;

% Extract MI per mouse for early (ep1-3) and late (ep4-6)
MI_WT_early  = zeros(n_wt,  1);
MI_WT_late   = zeros(n_wt,  1);
MI_HOM_early = zeros(n_hom, 1);
MI_HOM_late  = zeros(n_hom, 1);

for m = 1:n_wt
    MI_WT_early(m) = mean(PAC_WT.(chan){m}(theta_idx, gamma_idx, 1:2), 'all');
    MI_WT_late(m)  = mean(PAC_WT.(chan){m}(theta_idx, gamma_idx, 5:6), 'all');
end
for m = 1:n_hom
    MI_HOM_early(m) = mean(PAC_HOM.(chan){m}(theta_idx, gamma_idx, 1:2), 'all');
    MI_HOM_late(m)  = mean(PAC_HOM.(chan){m}(theta_idx, gamma_idx, 5:6), 'all');
end

% Early: WT vs Q331K
[~, p_early, ~, stats_early] = ttest2(MI_WT_early, MI_HOM_early, 'Vartype', 'unequal');
fprintf('\n=== Early REM (Ep1-3): WT vs Q331K ===\n');
fprintf('WT=%.4f±%.4f  Q331K=%.4f±%.4f  t(%g)=%.3f  p=%.4f\n', ...
    mean(MI_WT_early), std(MI_WT_early), ...
    mean(MI_HOM_early), std(MI_HOM_early), ...
    stats_early.df, stats_early.tstat, p_early);

% Late: WT vs Q331K
[~, p_late, ~, stats_late] = ttest2(MI_WT_late, MI_HOM_late, 'Vartype', 'unequal');
fprintf('\n=== Late REM (Ep4-6): WT vs Q331K ===\n');
fprintf('WT=%.4f±%.4f  Q331K=%.4f±%.4f  t(%g)=%.3f  p=%.4f\n', ...
    mean(MI_WT_late), std(MI_WT_late), ...
    mean(MI_HOM_late), std(MI_HOM_late), ...
    stats_late.df, stats_late.tstat, p_late);


% Violin plot - Early vs Late REM
colors_wt  = [0.25 0.70 0.55];   % green
colors_hom = [0.85 0.25 0.55];   % magenta

% titles  = {'Early REM (Ep1-2)', 'Late REM (Ep5-6)'};
titles  = {'Early REM', 'Late REM'};
MI_WT_list  = {MI_WT_early,  MI_WT_late};
MI_HOM_list = {MI_HOM_early, MI_HOM_late};
p_list      = {p_early,      p_late};

figure('Position', [100 100 700 400]);

for i = 1:2
    subplot(1, 2, i);
    hold on;

    wt_data  = MI_WT_list{i};
    hom_data = MI_HOM_list{i};
    p_val    = p_list{i};

    simple_violin({wt_data},  1, colors_wt,  0.3);
    simple_violin({hom_data}, 2, colors_hom, 0.3);

    % Individual points
    scatter(ones(size(wt_data))*1,  wt_data,  40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors_wt,  'LineWidth', 1);
    scatter(ones(size(hom_data))*2, hom_data, 40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors_hom, 'LineWidth', 1);

    % Mean line
    plot([0.8 1.2], [mean(wt_data)  mean(wt_data)],  'k-', 'LineWidth', 2);
    plot([1.8 2.2], [mean(hom_data) mean(hom_data)], 'k-', 'LineWidth', 2);

    % p-value bar
    y_top = max([wt_data; hom_data]) + 0.05 * max([wt_data; hom_data]);
    text(1.5, y_top*1.01, sprintf('p = %.3f', p_val), 'HorizontalAlignment', 'center', 'FontSize', 16);

    xlim([0.5 2.5]); xticks([1 2]); xticklabels({'WT', 'Q331K'});
    if i == 1
        ylabel('Modulation Index', 'FontSize', 20, 'FontWeight', 'bold');
    end
    title(titles{i}, 'FontSize', 22, 'FontWeight', 'bold');
    % text(0.03, 1.05, char('e'+i-1), 'Units','normalized','FontWeight','bold','FontSize',24) % add panels
    set(gca, 'FontSize', 18); box off;
end

% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_4_e_f.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_4_e_f.tiff'), 'Resolution', 300);


%% Phase amplitude histograms and violin plots of modulation depth (figure 4 g-i)

clear all; close all

% Step 1: Load and process data
fs           = 1000;
file_path    = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM\EEG_arrays';
files        = dir(fullfile(file_path, 'EEG_*.mat'));
nrem_samples = 30 * fs;
epoch_len    = 10 * fs;

freq_phase  = 6;
freq_amp    = 65;
n_hist_bins = 50;
phase_edges = linspace(-180, 180, n_hist_bins+1);
bin_centers = phase_edges(1:end-1) + diff(phase_edges)/2;

% Storage
amp_by_phases_WT  = [];
amp_by_phases_HOM = [];

wt_count = 0; hom_count = 0;

for f = 1:length(files)
    load(fullfile(file_path, files(f).name));
    fprintf('Processing %s\n', files(f).name);

    if contains(files(f).name, 'WT')
        group = 'WT'; wt_count = wt_count + 1;
    else
        group = 'HOM'; hom_count = hom_count + 1;
    end

    rem_bouts   = EEG.V1_EcOG.NREM_REM;
    phase_mouse = [];
    amp_mouse   = [];

    for b = 1:length(rem_bouts)
        bout = rem_bouts{b};
        if length(bout) < nrem_samples + 6*epoch_len; continue; end
        for ep = 1:6
            seg = bout(nrem_samples+(ep-1)*epoch_len+1 : nrem_samples+ep*epoch_len);
            phasefilt   = filterFGx(double(seg), fs, freq_phase, freq_phase*0.4);
            phase       = angle(hilbert(phasefilt));
            ampfilt     = filterFGx(double(seg), fs, freq_amp, freq_amp*0.78);
            amp         = abs(hilbert(ampfilt)).^2;
            phase_mouse = [phase_mouse; phase(:)];
            amp_mouse   = [amp_mouse;   amp(:)];
        end
    end

    if isempty(phase_mouse); continue; end

    % Histogram per mouse
    phase_deg     = rad2deg(phase_mouse);
    amp_by_phases = zeros(1, n_hist_bins);
    for i = 1:n_hist_bins-1
        amp_by_phases(i) = mean(amp_mouse(phase_deg > phase_edges(i) & phase_deg < phase_edges(i+1)));
    end

    if strcmp(group, 'WT')
        amp_by_phases_WT = [amp_by_phases_WT; amp_by_phases];
    else
        amp_by_phases_HOM = [amp_by_phases_HOM; amp_by_phases];
    end
end

% Group averages
mean_amp_WT  = mean(amp_by_phases_WT,  1);
mean_amp_HOM = mean(amp_by_phases_HOM, 1);

% Modulation depth per mouse
mod_depth_WT  = (max(amp_by_phases_WT,  [], 2) - min(amp_by_phases_WT,  [], 2)) ./ mean(amp_by_phases_WT,  2);
mod_depth_HOM = (max(amp_by_phases_HOM, [], 2) - min(amp_by_phases_HOM, [], 2)) ./ mean(amp_by_phases_HOM, 2);

[~, p, ~, stats] = ttest2(mod_depth_WT, mod_depth_HOM, 'Vartype', 'unequal');
fprintf('\n=== Modulation Depth ===\n');
fprintf('WT=%.3f±%.3f  Q331K=%.3f±%.3f  t(%g)=%.3f  p=%.4f\n', ...
    mean(mod_depth_WT), std(mod_depth_WT), ...
    mean(mod_depth_HOM), std(mod_depth_HOM), ...
    stats.df, stats.tstat, p);

% Step 2: Three-panel figure

colors_wt  = [0.25 0.70 0.55];
colors_hom = [0.85 0.25 0.55];

figure('Position', [100 100 1400 400]);

% Panel g: WT histogram
subplot(1, 3, 1);
bar(bin_centers, mean_amp_WT, 'FaceColor', colors_wt, 'EdgeColor', 'none');
xlabel('Phase at 6 Hz (deg.)', 'FontSize', 30, 'FontWeight', 'bold');
ylabel(['Power at ' num2str(freq_amp) ' Hz'], 'FontSize', 30, 'FontWeight', 'bold');
title('WT', 'FontSize', 50,'FontWeight','bold');
xlim([-180 180]); xticks(-180:60:180); ylim([0 450]);
% text(-0.15, 1.05, 'g', 'Units','normalized','FontWeight','bold','FontSize',26)
set(gca, 'FontSize', 20); box off;

% Panel h: Q331K histogram
subplot(1, 3, 2);
bar(bin_centers, mean_amp_HOM, 'FaceColor', colors_hom, 'EdgeColor', 'none');
xlabel('Phase at 6 Hz (deg.)', 'FontSize', 30, 'FontWeight', 'bold');
% ylabel(['Power at ' num2str(freq_amp) ' Hz'], 'FontSize', 20,'FontWeight','bold');
title('Q331K', 'FontSize', 50, 'FontWeight', 'bold');
xlim([-180 180]); xticks(-180:60:180); ylim([0 450]);
% text(-0.15, 1.05, 'h', 'Units','normalized','FontWeight','bold','FontSize',26)
set(gca, 'FontSize', 20); box off;

% Panel i: Modulation Depth violin
subplot(1, 3, 3);
hold on;
simple_violin({mod_depth_WT},  1, colors_wt,  0.3);
simple_violin({mod_depth_HOM}, 2, colors_hom, 0.3);
scatter(ones(size(mod_depth_WT))*1,  mod_depth_WT,  40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors_wt,  'LineWidth', 1);
scatter(ones(size(mod_depth_HOM))*2, mod_depth_HOM, 40, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors_hom, 'LineWidth', 1);
plot([0.8 1.2], [mean(mod_depth_WT)  mean(mod_depth_WT)],  'k-', 'LineWidth', 2);
plot([1.8 2.2], [mean(mod_depth_HOM) mean(mod_depth_HOM)], 'k-', 'LineWidth', 2);
y_top = max([mod_depth_WT; mod_depth_HOM]) + 0.05;
text(1.5, y_top, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center', 'FontSize', 20);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'WT', 'Q331K'});
ylabel('Modulation Depth', 'FontSize', 40, 'FontWeight', 'bold');
%title('Modulation Depth', 'FontSize', 22, 'FontWeight', 'bold');
% text(-0.15, 1.05, 'i', 'Units','normalized','FontWeight','bold','FontSize',26)
set(gca, 'FontSize', 20, 'FontWeight', 'bold'); box off;


% Shared x-label under panels g and h only
% han = axes(gcf, 'Visible', 'off');
% han.Position = [0.15 0.07 0.44 0.02];   % spans panels 1-2
% han.XLabel.Visible = 'on';
% xlabel(han, ['Phase at ' num2str(freq_phase) ' Hz (deg.)'], 'FontSize', 22, 'FontWeight', 'bold');



% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_4_g_i.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_4_g_i.tiff'), 'Resolution', 300);


%% end
