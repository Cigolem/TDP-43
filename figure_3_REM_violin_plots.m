

clear all; close all;

%% Plot FC EcOG theta and high beta - light on

output_path = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM\EEG_arrays';
eeg_files   = dir(fullfile(output_path, 'EEG_*.mat'));


fs       = 1000;
winsize  = 4 * fs;
noverlap = winsize / 2;
nfft     = 2^(nextpow2(winsize));
rem_window  = [30001, 90000];
total_band  = [1 50];
channel = 'PFC_EcOG';


bands.delta   = [1  4];
bands.theta   = [4  8];
bands.sigma   = [8  12];
bands.low_beta    = [12 20];
bands.high_beta = [20 30];
bands.gamma    = [30 45];
band_names = fieldnames(bands);

color_wt  = [0.25 0.65 0.55];
color_hom = [0.9  0.2  0.7];

%% Loop over bands

figure('Position', [100 100 400 700]);

for band_idx = 1:length(band_names)
    band_name  = band_names{band_idx};
    band_range = bands.(band_name);
    
    WT_REM = []; HOM_REM = [];
    WT_count = 0; HOM_count = 0;
    
    for file_idx = 1:length(eeg_files)
        load(fullfile(output_path, eeg_files(file_idx).name));
        mouse_id = eeg_files(file_idx).name;
        
        if contains(mouse_id, 'WT')
            WT_count  = WT_count + 1;
            mouse_idx = WT_count;
            genotype  = 'WT';
        elseif contains(mouse_id, 'HOM')
            HOM_count = HOM_count + 1;
            mouse_idx = HOM_count;
            genotype  = 'HOM';
        else
            continue;
        end
        
        mouse_power_REM = [];
        
        for bout = 1:length(EEG.transitions)
            signal     = EEG.(channel).NREM_REM{1, bout};
            rem_signal = signal(rem_window(1):rem_window(2));
            
            [psd, freqs] = pwelch(rem_signal, winsize, noverlap, nfft, fs);
            band_power   = trapz(freqs(freqs >= band_range(1) & freqs <= band_range(2)), ...
                                 psd(freqs  >= band_range(1) & freqs <= band_range(2)));
            total_power  = trapz(freqs(freqs >= total_band(1) & freqs <= total_band(2)), ...
                                 psd(freqs  >= total_band(1) & freqs <= total_band(2)));
            rel_power    = (band_power / total_power) * 100;
            mouse_power_REM = [mouse_power_REM; rel_power];
        end
        
        if strcmp(genotype, 'WT')
            WT_REM(mouse_idx)  = mean(mouse_power_REM);
        else
            HOM_REM(mouse_idx) = mean(mouse_power_REM);
        end
    end
    
    %% t-test
    [~, p, ~, stats] = ttest2(WT_REM, HOM_REM, 'Vartype', 'unequal');
    fprintf('\n=== %s | %s (%d-%d Hz) ===\n', channel, band_name, band_range(1), band_range(2));
    fprintf('WT=%.2f±%.2f  Q331K=%.2f±%.2f  t(%g)=%.3f  p=%.4f\n', ...
        mean(WT_REM), std(WT_REM), mean(HOM_REM), std(HOM_REM), stats.df, stats.tstat, p);
    

    %% Violin plot
    
     plot_slot = find(strcmp(band_name, {'theta','high_beta'}));
          if isempty(plot_slot), continue; end

    subplot(2, 1, plot_slot); hold on;
   
     % Plot

    simple_violin({WT_REM},  1, color_wt,  0.35);
    simple_violin({HOM_REM}, 2, color_hom, 0.35);
    
    scatter(ones(1,WT_count),  WT_REM,  120, 'o', 'MarkerEdgeColor','k','MarkerFaceColor',color_wt,  'LineWidth',1.5,'MarkerFaceAlpha',0.7);
    scatter(2*ones(1,HOM_count), HOM_REM, 120, 'o', 'MarkerEdgeColor','k','MarkerFaceColor',color_hom, 'LineWidth',1.5,'MarkerFaceAlpha',0.7);
    
    % p-value text
    y_max = max([WT_REM HOM_REM]) * 1.1;
    text(0.6, 0.90, sprintf('p = %.3f', p), 'Units','normalized','HorizontalAlignment','right','FontSize',14);

    
    xlim([0.5 2.5]);
    xticks([1 2]); xticklabels({'WT', 'Q331K'});
    band_names_label = {'Theta', 'High Beta'};
    title(sprintf(' REM %s (%d-%d Hz)', band_names_label{plot_slot}, band_range(1), band_range(2)), 'FontSize', 18);
    text(0.003, 1.10, char('e'+plot_slot-1), 'Units','normalized','FontWeight','bold','FontSize',20) % add panels
    box off; set(gca, 'FontSize', 16, 'LineWidth', 1.5);
    

end


% Shared y-label between the two subplots
    han = axes(gcf, 'Visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Band Power (% total)', 'FontSize', 16, 'FontWeight', 'bold');
    

 % Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_3_e_f.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_3_e_f.tiff'), 'Resolution', 300);


clear all; close all;


%% Plot V1 EcOG theta and high beta - light on

output_path = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM\EEG_arrays';
eeg_files   = dir(fullfile(output_path, 'EEG_*.mat'));


fs       = 1000;
winsize  = 4 * fs;
noverlap = winsize / 2;
nfft     = 2^(nextpow2(winsize));
rem_window  = [30001, 90000];
total_band  = [1 50];
channel = 'V1_EcOG';


bands.delta   = [1  4];
bands.theta   = [4  8];
bands.sigma   = [8  12];
bands.low_beta    = [12 20];
bands.high_beta = [20 30];
bands.gamma    = [30 45];
band_names = fieldnames(bands);

color_wt  = [0.25 0.65 0.55];
color_hom = [0.9  0.2  0.7];

%% Loop over bands

figure('Position', [100 100 400 700]);

for band_idx = 1:length(band_names)
    band_name  = band_names{band_idx};
    band_range = bands.(band_name);
    
    WT_REM = []; HOM_REM = [];
    WT_count = 0; HOM_count = 0;
    
    for file_idx = 1:length(eeg_files)
        load(fullfile(output_path, eeg_files(file_idx).name));
        mouse_id = eeg_files(file_idx).name;
        
        if contains(mouse_id, 'WT')
            WT_count  = WT_count + 1;
            mouse_idx = WT_count;
            genotype  = 'WT';
        elseif contains(mouse_id, 'HOM')
            HOM_count = HOM_count + 1;
            mouse_idx = HOM_count;
            genotype  = 'HOM';
        else
            continue;
        end
        
        mouse_power_REM = [];
        
        for bout = 1:length(EEG.transitions)
            signal     = EEG.(channel).NREM_REM{1, bout};
            rem_signal = signal(rem_window(1):rem_window(2));
            
            [psd, freqs] = pwelch(rem_signal, winsize, noverlap, nfft, fs);
            band_power   = trapz(freqs(freqs >= band_range(1) & freqs <= band_range(2)), ...
                                 psd(freqs  >= band_range(1) & freqs <= band_range(2)));
            total_power  = trapz(freqs(freqs >= total_band(1) & freqs <= total_band(2)), ...
                                 psd(freqs  >= total_band(1) & freqs <= total_band(2)));
            rel_power    = (band_power / total_power) * 100;
            mouse_power_REM = [mouse_power_REM; rel_power];
        end
        
        if strcmp(genotype, 'WT')
            WT_REM(mouse_idx)  = mean(mouse_power_REM);
        else
            HOM_REM(mouse_idx) = mean(mouse_power_REM);
        end
    end
    
    %% t-test
    [~, p, ~, stats] = ttest2(WT_REM, HOM_REM, 'Vartype', 'unequal');
    fprintf('\n=== %s | %s (%d-%d Hz) ===\n', channel, band_name, band_range(1), band_range(2));
    fprintf('WT=%.2f±%.2f  Q331K=%.2f±%.2f  t(%g)=%.3f  p=%.4f\n', ...
        mean(WT_REM), std(WT_REM), mean(HOM_REM), std(HOM_REM), stats.df, stats.tstat, p);
    

    %% Violin plot
    
     plot_slot = find(strcmp(band_name, {'theta','low_beta'}));
          if isempty(plot_slot), continue; end

    subplot(2, 1, plot_slot); hold on;
   
     % Plot

    simple_violin({WT_REM},  1, color_wt,  0.35);
    simple_violin({HOM_REM}, 2, color_hom, 0.35);
    
    scatter(ones(1,WT_count),  WT_REM,  120, 'o', 'MarkerEdgeColor','k','MarkerFaceColor',color_wt,  'LineWidth',1.5,'MarkerFaceAlpha',0.7);
    scatter(2*ones(1,HOM_count), HOM_REM, 120, 'o', 'MarkerEdgeColor','k','MarkerFaceColor',color_hom, 'LineWidth',1.5,'MarkerFaceAlpha',0.7);
    
    % p-value text
    y_max = max([WT_REM HOM_REM]) * 1.1;
    text(0.6, 0.90, sprintf('p = %.3f', p), 'Units','normalized','HorizontalAlignment','right','FontSize',14);

    
    xlim([0.5 2.5]);
    xticks([1 2]); xticklabels({'WT', 'Q331K'});
    band_names_label = {'Theta', 'Low Beta'};
    title(sprintf(' REM %s (%d-%d Hz)', band_names_label{plot_slot}, band_range(1), band_range(2)), 'FontSize', 18);
    text(0.003, 1.10, char('k'+plot_slot-1), 'Units','normalized','FontWeight','bold','FontSize',20) % add panels
    box off; set(gca, 'FontSize', 16, 'LineWidth', 1.5);
    

end


% Shared y-label between the two subplots
    han = axes(gcf, 'Visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Band Power (% total)', 'FontSize', 16, 'FontWeight', 'bold');
    

 % Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_3_k_l.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_3_k_l.tiff'), 'Resolution', 300);

%% end