
%% Wavelet analysis - RELATIVE POWER (Instantaneous Normalization)
% This script computes relative power at each time point
% Relative power = (power at freq) / (total 1-30Hz power at that time) * 100
% 


clear all
close all

output_path = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM\EEG_arrays';
eeg_files = dir(fullfile(output_path, 'EEG_*.mat'));

% Initialize storage - each mouse will be stored separately
WT_PFC_LFP_mice = {};
WT_PFC_ECoG_mice = {};
WT_V1_ECoG_mice = {};
HOM_PFC_LFP_mice = {};
HOM_PFC_ECoG_mice = {};
HOM_V1_ECoG_mice = {};

WT_count = 0;
HOM_count = 0;

fprintf('Processing files with RELATIVE POWER calculation...\n');

for file_idx = 1:length(eeg_files)
    % Load EEG structure
    load(fullfile(output_path, eeg_files(file_idx).name));
    mouse_id = eeg_files(file_idx).name;
    num_transitions = length(EEG.transitions);
    fprintf('Processing %s: %d transitions\n', mouse_id, num_transitions);
    
    % Determine genotype
    if contains(mouse_id, 'WT')
        genotype = 'WT';
        WT_count = WT_count + 1;
        mouse_idx = WT_count;
    elseif contains(mouse_id, 'HOM')
        genotype = 'HOM';
        HOM_count = HOM_count + 1;
        mouse_idx = HOM_count;
    else
        continue;
    end
    
    % Process PFC LFP
    Power_PFC_LFP = {};
    for bout = 1:num_transitions
        signal = EEG.PFC_LFP.NREM_REM{1, bout};
        [Power_rel, ~] = Wavelet_CG_relative(signal);
        Power_PFC_LFP{bout} = Power_rel;  % Store relative power (non-log)
    end
    % Average across bouts for THIS mouse
    Power_PFC_LFP_mean = mean(cat(3, Power_PFC_LFP{:}), 3);
    
    % Process PFC ECoG
    Power_PFC_ECoG = {};
    for bout = 1:num_transitions
        signal = EEG.PFC_EcOG.NREM_REM{1, bout};
        [Power_rel, ~] = Wavelet_CG_relative(signal);
        Power_PFC_ECoG{bout} = Power_rel;
    end
    Power_PFC_ECoG_mean = mean(cat(3, Power_PFC_ECoG{:}), 3);
    
    % Process V1 ECoG
    Power_V1_ECoG = {};
    for bout = 1:num_transitions
        signal = EEG.V1_EcOG.NREM_REM{1, bout};
        [Power_rel, ~] = Wavelet_CG_relative(signal);
        Power_V1_ECoG{bout} = Power_rel;
    end
    Power_V1_ECoG_mean = mean(cat(3, Power_V1_ECoG{:}), 3);
    
    % Store this mouse's data by genotype
    if strcmp(genotype, 'WT')
        WT_PFC_LFP_mice{mouse_idx} = Power_PFC_LFP_mean;
        WT_PFC_ECoG_mice{mouse_idx} = Power_PFC_ECoG_mean;
        WT_V1_ECoG_mice{mouse_idx} = Power_V1_ECoG_mean;
    else
        HOM_PFC_LFP_mice{mouse_idx} = Power_PFC_LFP_mean;
        HOM_PFC_ECoG_mice{mouse_idx} = Power_PFC_ECoG_mean;
        HOM_V1_ECoG_mice{mouse_idx} = Power_V1_ECoG_mean;
    end
end


fprintf('\nProcessing complete: %d WT mice, %d HOM mice\n', WT_count, HOM_count);

% Save individual mouse data (relative power version)
save_path_output = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM';

% Save the data
save(fullfile(save_path_output, 'Morlet_relative_averaged_transitions_lights_on.mat'), ...
    'WT_PFC_LFP_mice', 'WT_PFC_ECoG_mice', 'WT_V1_ECoG_mice', ...
    'HOM_PFC_LFP_mice', 'HOM_PFC_ECoG_mice', 'HOM_V1_ECoG_mice', ...
    'WT_count', 'HOM_count', '-v7.3');

fprintf('Saved RELATIVE POWER data to: %s\n', fullfile(save_path_output, 'Morlet_relative_averaged_transitions_lights_on.mat'));

%% plot
clear

load('D:\TDP\lights_on\shorter_extracts\ALL_NREM_REM\Morlet_relative_averaged_transitions_lights_on.mat')

clearvars -except HOM_count WT_count HOM_V1_ECoG_mice HOM_PFC_ECoG_mice HOM_PFC_LFP_mice WT_V1_ECoG_mice...
    WT_PFC_ECoG_mice WT_PFC_LFP_mice

% Get frequency and time vectors
maxfrequency = 30;
numoctaves = 6;
numsuboctaves = 4;
dj = 1/numsuboctaves;
s0 = 1/maxfrequency;
j1 = numoctaves/dj;
Freq = 2.^(log2(s0) + (0:j1)*dj);
Freq = 1./Freq;
Freq = fliplr(Freq);

Fs = 1000;
[~, ntime] = size(WT_PFC_LFP_mice{1});
time = (0:ntime-1)/Fs;

%% Plot settings
% Define channels and data
channels = {'PFC_LFP', 'PFC_ECoG', 'V1_ECoG'};
channel_data = struct();
channel_data.PFC_LFP  = struct('WT', {WT_PFC_LFP_mice},  'HOM', {HOM_PFC_LFP_mice});
channel_data.PFC_ECoG = struct('WT', {WT_PFC_ECoG_mice}, 'HOM', {HOM_PFC_ECoG_mice});
channel_data.V1_ECoG  = struct('WT', {WT_V1_ECoG_mice},  'HOM', {HOM_V1_ECoG_mice});

% Set color bar limits
clims_high_all = struct('PFC_LFP',[0 2],  'PFC_ECoG',[0 3],  'V1_ECoG',[0 3]);
clims_low_all  = struct('PFC_LFP',[0 15], 'PFC_ECoG',[0 10], 'V1_ECoG',[0 18]);

% Create plots for all three channels

% for ch_idx = 2:length(channels)
    for ch_idx = 2 % plot one channel at a time )
    channel = channels{ch_idx};
    channel_label = strrep(strrep(channel, '_', ' '), 'PFC', 'FC');
    fprintf('\nPlotting %s...\n', channel_label);

    WT_data  = channel_data.(channel).WT;
    HOM_data = channel_data.(channel).HOM;

    nWT_mice  = length(WT_data);
    nHOM_mice = length(HOM_data);

    [nfreqs, ntime] = size(WT_data{1});

    WT_power  = NaN(nfreqs, ntime, nWT_mice);
    HOM_power = NaN(nfreqs, ntime, nHOM_mice);

    for i = 1:nWT_mice,  WT_power(:,:,i)  = WT_data{i};  end
    for i = 1:nHOM_mice, HOM_power(:,:,i) = HOM_data{i}; end

    WT_mn  = mean(WT_power,  3);
    HOM_mn = mean(HOM_power, 3);

    % Separate color limits for each channel
    clims_high = clims_high_all.(channel);
    clims_low  = clims_low_all.(channel);

    figure('Position', [100 100 1400 700])
    colormap('jet')

    genotype_labels = {'WT', 'Q331K'};
    data_list = {WT_mn, HOM_mn};

    for g = 1:2
        dat = data_list{g};

        % --- Top panel: low frequency 1-12 Hz ---
        subplot(2, 2, g)
        contourf(time, Freq, dat, 50, 'linestyle', 'none')
        ylim([1 12])
        xlim([0 90]); xticks(0:10:90)
        caxis(clims_low)
        colorbar
        % xlabel('Time (sec)', 'FontSize', 18)
        % title(genotype_labels{g}, 'FontSize', 20)
        text(40, 13, genotype_labels{g}, 'HorizontalAlignment','center','FontWeight','bold','FontSize',20)
        % NREM/REM labels centred in each half (top panels only)
        text(17, 13, 'NREM', 'HorizontalAlignment','center','FontWeight','bold','FontSize',16)
        text(62, 13, 'REM',  'HorizontalAlignment','center','FontWeight','bold','FontSize',16)
        set(gca, 'FontSize', 18, 'TickDir', 'out', 'box', 'off')
        text(0.03, 1.10, char('a'+2*(g-1)), 'Units','normalized','FontWeight','bold','FontSize',20) % add panels
        hold on; xline(30, 'k--', 'LineWidth', 2); hold off

        % --- Bottom panel: high frequency 12-30 Hz ---
        subplot(2, 2, g+2)
        contourf(time, Freq, dat, 50, 'linestyle', 'none')
        ylim([12 30])
        yticks([12 15 20 25 30])
        xlim([0 90]); xticks(0:10:90)
        caxis(clims_high)
        colorbar
        xlabel('Time (sec)', 'FontSize', 18)
        set(gca, 'FontSize', 18, 'TickDir', 'out', 'box', 'off')
        text(0.03, 1.10, char('b'+2*(g-1)), 'Units','normalized','FontWeight','bold','FontSize',20) % add panels
        hold on; xline(30, 'k--', 'LineWidth', 2); hold off

    end  % end g loop

    % Shared labels
    annotation('textbox', [0.99 0.17 0.3 0.2], 'String', 'Relative Power (%)', ...
        'FontSize', 16, 'Rotation', 90, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    annotation('textbox', [0.12 0.17 0.3 0.2], 'String', 'Frequency (Hz)', ...
        'FontSize', 18, 'Rotation', 90, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    sgtitle(channel_label, 'FontSize', 22)

    % Save figure — inside ch_idx loop so all three channels are saved

    save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';


   % Save as png + tiff, tight margins
   exportgraphics(gcf, fullfile(save_path_plot, sprintf('figure_3_a_d.png', channel)),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, sprintf('figure_3_a_d.tiff', channel)), 'Resolution', 300);
   

end  % end ch_idx loop

%%%%%%%%%%%%%  other channel %%%%%%%%%%%%

clearvars -except HOM_count WT_count HOM_V1_ECoG_mice HOM_PFC_ECoG_mice HOM_PFC_LFP_mice WT_V1_ECoG_mice...
    WT_PFC_ECoG_mice WT_PFC_LFP_mice

% Get frequency and time vectors
maxfrequency = 30;
numoctaves = 6;
numsuboctaves = 4;
dj = 1/numsuboctaves;
s0 = 1/maxfrequency;
j1 = numoctaves/dj;
Freq = 2.^(log2(s0) + (0:j1)*dj);
Freq = 1./Freq;
Freq = fliplr(Freq);

Fs = 1000;
[~, ntime] = size(WT_PFC_LFP_mice{1});
time = (0:ntime-1)/Fs;

%% Plot settings
% Define channels and data
channels = {'PFC_LFP', 'PFC_ECoG', 'V1_ECoG'};
channel_data = struct();
channel_data.PFC_LFP  = struct('WT', {WT_PFC_LFP_mice},  'HOM', {HOM_PFC_LFP_mice});
channel_data.PFC_ECoG = struct('WT', {WT_PFC_ECoG_mice}, 'HOM', {HOM_PFC_ECoG_mice});
channel_data.V1_ECoG  = struct('WT', {WT_V1_ECoG_mice},  'HOM', {HOM_V1_ECoG_mice});

% Set color bar limits
clims_high_all = struct('PFC_LFP',[0 2],  'PFC_ECoG',[0 3],  'V1_ECoG',[0 3]);
clims_low_all  = struct('PFC_LFP',[0 15], 'PFC_ECoG',[0 10], 'V1_ECoG',[0 18]);

% Create plots for all three channels

% for ch_idx = 2:length(channels)
  for ch_idx = 3 % plot one channel at a time )
    channel = channels{ch_idx};
    channel_label = strrep(strrep(channel, '_', ' '), 'PFC', 'FC');
    fprintf('\nPlotting %s...\n', channel_label);

    WT_data  = channel_data.(channel).WT;
    HOM_data = channel_data.(channel).HOM;

    nWT_mice  = length(WT_data);
    nHOM_mice = length(HOM_data);

    [nfreqs, ntime] = size(WT_data{1});

    WT_power  = NaN(nfreqs, ntime, nWT_mice);
    HOM_power = NaN(nfreqs, ntime, nHOM_mice);

    for i = 1:nWT_mice,  WT_power(:,:,i)  = WT_data{i};  end
    for i = 1:nHOM_mice, HOM_power(:,:,i) = HOM_data{i}; end

    WT_mn  = mean(WT_power,  3);
    HOM_mn = mean(HOM_power, 3);

    % Separate color limits for each channel
    clims_high = clims_high_all.(channel);
    clims_low  = clims_low_all.(channel);

    figure('Position', [100 100 1400 700])
    colormap('jet')

    genotype_labels = {'WT', 'Q331K'};
    data_list = {WT_mn, HOM_mn};

    for g = 1:2
        dat = data_list{g};

        % --- Top panel: low frequency 1-12 Hz ---
        subplot(2, 2, g)
        contourf(time, Freq, dat, 50, 'linestyle', 'none')
        ylim([1 12])
        xlim([0 90]); xticks(0:10:90)
        caxis(clims_low)
        colorbar
        % xlabel('Time (sec)', 'FontSize', 18)
        % title(genotype_labels{g}, 'FontSize', 20)
        text(40, 13, genotype_labels{g}, 'HorizontalAlignment','center','FontWeight','bold','FontSize',20)
        % NREM/REM labels centred in each half (top panels only)
        text(17, 13, 'NREM', 'HorizontalAlignment','center','FontWeight','bold','FontSize',16)
        text(62, 13, 'REM',  'HorizontalAlignment','center','FontWeight','bold','FontSize',16)
        set(gca, 'FontSize', 18, 'TickDir', 'out', 'box', 'off')
        text(0.03, 1.10, char('g'+2*(g-1)), 'Units','normalized','FontWeight','bold','FontSize',20) % add panels
        hold on; xline(30, 'k--', 'LineWidth', 2); hold off

        % --- Bottom panel: high frequency 12-30 Hz ---
        subplot(2, 2, g+2)
        contourf(time, Freq, dat, 50, 'linestyle', 'none')
        ylim([12 30])
        yticks([12 15 20 25 30])
        xlim([0 90]); xticks(0:10:90)
        caxis(clims_high)
        colorbar
        xlabel('Time (sec)', 'FontSize', 18)
        set(gca, 'FontSize', 18, 'TickDir', 'out', 'box', 'off')
        text(0.03, 1.10, char('h'+2*(g-1)), 'Units','normalized','FontWeight','bold','FontSize',20) % add panels
        hold on; xline(30, 'k--', 'LineWidth', 2); hold off

    end  % end g loop

    % Shared labels
    annotation('textbox', [0.99 0.17 0.3 0.2], 'String', 'Relative Power (%)', ...
        'FontSize', 16, 'Rotation', 90, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    annotation('textbox', [0.12 0.17 0.3 0.2], 'String', 'Frequency (Hz)', ...
        'FontSize', 18, 'Rotation', 90, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    sgtitle(channel_label, 'FontSize', 22)

    % Save figure 
    save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';


   % Save as png + tiff, tight margins
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_3_g_j.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'figure_3_g_j.tiff'), 'Resolution', 300);


end  % end ch_idx loop
