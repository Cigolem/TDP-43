
close all; clear all;

%% Set parameters
file_path = 'D:\TDP\lights_on\shorter_extracts\ALL_NREM\EEG_arrays';
file_name = 'EEG_WT3.mat';  % Pick a mouse

channel_name = 'PFC_EcOG';   % Pick a channel
%channel_name = 'PFC_EcOG';   % Pick a channel

% initiate result str array to store all outputs
results = struct;
results.mouse       = file_name;
results.channel     = channel_name;

%% Step 1: Load EEG file
load(fullfile(file_path, file_name));
fs = 1000;
n_bouts = length(EEG.(channel_name).NREM_eeg);

%% Step 2: Design Butterworth filter (0.5-4 Hz, zero-phase, 4th order)
order = 3; % 3 used as trade-off between filter sharpness and stability
[b, a] = butter(order, [0.5 4]/(fs/2), 'bandpass');

%% Step 3: Apply filter to all NREM bouts

filtered_data = cell(n_bouts, 1);

for i = 1:n_bouts
    eeg_signal = EEG.(channel_name).NREM_eeg{i};
    filtered_data{i} = filtfilt(b,a, eeg_signal);
end

%% Step 4: Find positive-to-negative zero-crossings (Niethard method)

for bout = 1:n_bouts
    signal{bout} = filtered_data{bout};
    
    % Find where signal crosses from positive to negative
    % diff(sign(signal)) == -2 means: was +1, now -1 (positive to negative)
    pos_to_neg_crossings{bout} = find(diff(sign(signal{bout})) == -2);
end

%% Step 5: Find slow wave candidates with duration and amplitude criteria

min_duration = 0.4;  % seconds - full cycle
max_duration = 2.0;  % seconds

% Initialize storage for all bouts
all_slow_waves = cell(n_bouts, 1);

% Collect ALL potential SWAs to calculate amplitude thresholds
all_min_amps = [];
all_peak_to_peak = [];

for bout = 1:n_bouts
    crossings = pos_to_neg_crossings{bout};
    
    % Look at intervals between consecutive crossings
    for i = 1:(length(crossings)-1)
        onset_idx = crossings(i);      % First pos→neg crossing
        end_idx = crossings(i+1);      % Next pos→neg crossing
        duration = (end_idx - onset_idx) / fs;
        
        % Check duration criterion
        if duration >= min_duration && duration <= max_duration
            segment = signal{bout}(onset_idx:end_idx);
            
            % Find minimum (trough) and maximum (positive peak)
            [min_val, ~] = min(segment);
            [max_val, ~] = max(segment);
            
            % Store for threshold calculation
            all_min_amps = [all_min_amps; abs(min_val)];
            all_peak_to_peak = [all_peak_to_peak; (max_val - min_val)];
        end
    end
end

% Calculate 66.6% thresholds
min_amp_threshold = 0.666 * mean(all_min_amps);
peak_to_peak_threshold = 0.666 * mean(all_peak_to_peak);

%% Step 6: Store valid SWAs with full event info

for bout = 1:n_bouts
    crossings = pos_to_neg_crossings{bout};
    swa_events = [];
    
    for i = 1:(length(crossings)-1)
        onset_idx = crossings(i);
        end_idx = crossings(i+1);
        duration = (end_idx - onset_idx) / fs;
        
        if duration >= min_duration && duration <= max_duration
            segment = signal{bout}(onset_idx:end_idx);
            
            % Find trough and peak LOCATIONS within this segment
            [min_val, min_idx_rel] = min(segment);
            [max_val, max_idx_rel] = max(segment);
            
            % Convert to absolute indices
            trough_idx = onset_idx + min_idx_rel - 1;
            peak_idx = onset_idx + max_idx_rel - 1;
            
            % Check amplitude criteria
            peak_to_peak = max_val - min_val;
            if abs(min_val) >= min_amp_threshold && peak_to_peak >= peak_to_peak_threshold
                % Store: [onset, end, trough_idx, peak_idx, trough_val, peak_val]
                swa_events = [swa_events; onset_idx, end_idx, trough_idx, peak_idx, min_val, max_val];
            end
        end
    end
    
    all_slow_waves{bout} = swa_events;
    % fprintf('Bout %d: %d slow waves detected\n', bout, size(swa_events, 1));
end


%% Step 7: Spindle detection 

clearvars -except results EEG all_slow_waves fs channel_name n_bouts...
          pos_to_neg_crossings signal file_name;


% Step 7.1: Design Butterworth filter (10-15 Hz, zero-phase, 4th order)

order = 4;
[b, a] = butter(order, [10 15]/(fs/2), 'bandpass');



% Step 7.2: Apply filter to all NREM bouts

filtered_data = cell(n_bouts, 1);

for i = 1:n_bouts
    eeg_signal = EEG.(channel_name).NREM_eeg{i};
    filtered_data{i} = filtfilt(b,a, eeg_signal);
end


% Step 7.3: Compute Hilbert amplitude envelope
hilbert_env = cell(n_bouts, 1);
for i = 1:n_bouts
    analytic_signal = hilbert(filtered_data{i});
    hilbert_env{i} = abs(analytic_signal);
end

smooth_window = round(0.1 * fs);  % 100ms
for i = 1:n_bouts
    hilbert_env{i} = movmean(hilbert_env{i}, smooth_window);
end



% Step 7.4: Calculate amplitude thresholds (percentile of all envelope)
all_envelope = [];
for i = 1:n_bouts
    all_envelope = [all_envelope, hilbert_env{i}];
end

mean_env = mean(all_envelope);
std_env = std(all_envelope);

upper_threshold = prctile(all_envelope, 97);  
lower_threshold = prctile(all_envelope, 85);  



% Step 7.5: Detect spindle using amplitude threshold 

spindles = cell(n_bouts, 1);
for bout = 1:n_bouts
    env_sig = hilbert_env{bout};       
    spindle_events = [];
    
    above_lower = env_sig > lower_threshold;
    
    diff_lower = diff([0, above_lower, 0]);
    starts = find(diff_lower == 1);
    ends   = find(diff_lower == -1) - 1;
    ends   = min(ends, length(env_sig));  % <-- safety cap
    
    for j = 1:length(starts)
        segment = env_sig(starts(j):ends(j));   % <-- uses env_sig
        if any(segment > upper_threshold)
            duration = (ends(j) - starts(j) + 1) / fs;
            spindle_events = [spindle_events; starts(j), ends(j), duration];
        end
    end
    spindles{bout} = spindle_events;
end

% Step 7.6: Merge spindles with ISI between 0.2 seconds

max_ISI = 0.2;   % seconds

for bout = 1:n_bouts
    if isempty(spindles{bout})
        continue;
    end

    events = spindles{bout};
    merged = [];
    i = 1;

    while i <= size(events, 1)
        current_start = events(i, 1);
        current_end = events(i, 2);

        % Check ISI with next spindle
        while i < size(events, 1)
            ISI = (events(i+1, 1) - current_end) / fs;
            if ISI <= max_ISI
                current_end = events(i+1, 2);  % Extend to next spindle
                i = i + 1;
            else
                break;
            end
        end

        duration = (current_end - current_start + 1) / fs;
        merged = [merged; current_start, current_end, duration];
        i = i + 1;
    end

    spindles{bout} = merged;
end


% Step 7.7: Apply duration criteria (0.5-3 seconds)
for bout = 1:n_bouts
    if isempty(spindles{bout})
        continue;
    end
    events = spindles{bout};
    valid = events(:, 3) >= 0.5 & events(:, 3) <= 3;
    spindles{bout} = events(valid, :);
end


%% Step 8: Classify events into 3 categories Coupling window [trough, peak]

clearvars -except results EEG all_slow_waves spindles fs channel_name n_bouts signal file_name

solitary_swa = cell(n_bouts, 1);
solitary_spindles = cell(n_bouts, 1);
coupled_events = cell(n_bouts, 1);

for bout = 1:n_bouts
    
    swa_events = all_slow_waves{bout};
    sp_events = spindles{bout};
    
    swa_is_coupled = false(size(swa_events, 1), 1);
    sp_is_coupled = false(size(sp_events, 1), 1);
    
    coupled_pairs = [];
    
    for swa_idx = 1:size(swa_events, 1)
        swa_onset = swa_events(swa_idx, 1);
        swa_trough = swa_events(swa_idx, 3);
        swa_peak = swa_events(swa_idx, 4);
        
        % Coupling window: [trough, peak + 100ms]
        window_start = swa_trough;
        % window_end = swa_peak + round(0.1 * fs);  % +100ms
        window_end = swa_peak;  % no extension
        
        for sp_idx = 1:size(sp_events, 1)
            sp_onset = sp_events(sp_idx, 1);
            
            if sp_onset >= window_start && sp_onset <= window_end
                swa_is_coupled(swa_idx) = true;
                sp_is_coupled(sp_idx) = true;
                coupled_pairs = [coupled_pairs; swa_idx, sp_idx];
            end
        end
    end
    
    solitary_swa{bout} = swa_events(~swa_is_coupled, :);
    solitary_spindles{bout} = sp_events(~sp_is_coupled, :);
    coupled_events{bout} = coupled_pairs;
    
end



%% Step 9: Visualize example coupled event


% Collect all coupled events and rank by spindle duration
all_candidates = [];
for b = 1:n_bouts
    for e = 1:size(coupled_events{b}, 1)
        sp_idx = coupled_events{b}(e, 2);
        sp_dur = spindles{b}(sp_idx, 3);
        all_candidates = [all_candidates; b, e, sp_dur];
    end
end

% Sort by duration
all_candidates = sortrows(all_candidates, 3);

% Print all
fprintf('Rank | Bout | Event | Duration\n');
for k = 1:size(all_candidates, 1)
    fprintf('%4d | %4d | %5d | %.3f s\n', k, all_candidates(k,1), all_candidates(k,2), all_candidates(k,3));
end

% Select by rank - change this number
rank_select = 48; % this row (or rank) has the bout 7 and event 3 from WT 3
bout      = all_candidates(rank_select, 1);
event_num = all_candidates(rank_select, 2);
fprintf('\nSelected rank %d: Bout %d, Event %d, Duration %.3f s\n', rank_select, bout, event_num, all_candidates(rank_select,3));

% selected: WT 1 bout 7 event 3  Duration 0.902 sec

if size(coupled_events{bout}, 1) == 0
    error('No coupled events in bout %d, try different bout', bout);
end

% Filtered signals
order_swa = 2;
[b_swa, a_swa] = butter(order_swa, [0.5 4]/(fs/2), 'bandpass');
swa_filt = filtfilt(b_swa, a_swa, EEG.(channel_name).NREM_eeg{bout});

order_sp = 4;
[b_sp, a_sp] = butter(order_sp, [10 15]/(fs/2), 'bandpass');
sp_filt = filtfilt(b_sp, a_sp, EEG.(channel_name).NREM_eeg{bout});

% Get event indices
swa_idx  = coupled_events{bout}(event_num, 1);
sp_idx   = coupled_events{bout}(event_num, 2);

swa = all_slow_waves{bout}(swa_idx, :);
sp  = spindles{bout}(sp_idx, :);

onset    = swa(1); 
trough   = swa(3); 
peak     = swa(4);
sp_onset = sp(1);  
sp_end   = sp(2);

% Print info
fprintf('Bout %d, Event %d\n', bout, event_num);
fprintf('SWA onset: %.3f s | Trough: %.3f s | Peak: %.3f s\n', onset/fs, trough/fs, peak/fs);
fprintf('Spindle onset: %.3f s | End: %.3f s | Duration: %.3f s\n', sp_onset/fs, sp_end/fs, sp(3));

% Window around event
win_start = max(1, onset - fs);
win_end   = min(length(swa_filt), sp_end + fs);

% Time axis in seconds — referenced to sample indices
t = (win_start:win_end) / fs;

fprintf('t starts: %.3f | t ends: %.3f\n', t(1), t(end));
fprintf('patch from %.3f to %.3f\n', trough/fs, peak/fs);

%% Plot Representative Spindle - SWA from FC EcOG WT 3 bout 7 event number 3 (supplementary figure 3)

figure('Position', [100 100 900 450]);

h1 = plot(t, swa_filt(win_start:win_end), 'k', 'LineWidth', 2);
hold on;
h2 = plot(t, sp_filt(win_start:win_end), 'b', 'LineWidth', 1);

% Markers
xline(onset/fs,    'k--', 'LineWidth', 2.5, 'Label', 'SWA onset', 'FontSize', 14);
xline(sp_onset/fs, 'm--', 'LineWidth', 2.5, 'Label', 'Sp onset',  'FontSize', 14);
xline(sp_end/fs,   'm--', 'LineWidth', 2.5, 'Label', 'Sp end',    'FontSize', 14);

% Set ylim before patch
xlim([win_start/fs, win_end/fs]);
ylim([-200 250]);

% Draw patch using same time reference as t
h3 = patch([trough/fs, peak/fs, peak/fs, trough/fs], ...
           [-200, -200, 250, 250], ...
           [1 1 0.6], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Bring lines to front
uistack(findobj(gca, 'Type', 'line'), 'top');

% Labels
xlabel('Time (s)', 'FontSize', 20, 'FontWeight','Bold');
ylabel('Amplitude (\muV)', 'FontSize',  20, 'FontWeight','Bold');
title('FC EcOG - Representative Coupled SWA-Spindle', ...
      'FontSize', 22, 'FontWeight', 'bold');
legend([h1 h2 h3], ...
       'SWA filtered (0.5-4 Hz)', ...
       'Spindle filtered (10-15 Hz)', ...
       'Coupling window (trough to peak)', ...
       'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 14);
box off;
set(gca, 'FontSize', 20, 'LineWidth', 1.5);


% Save figure 
   save_path_plot = 'C:\Users\cigde\Dropbox\PDC_Projects\Jemeen\Sleep_paper\figures\';
   
   exportgraphics(gcf, fullfile(save_path_plot, 'supplementary_figure_1.png'),  'Resolution', 300);
   exportgraphics(gcf, fullfile(save_path_plot, 'supplementary_figure_1.tiff'), 'Resolution', 300);

%% end
