
function [results] = detect_swa_sp_core_5(file_path, file_name, channel_name)
% prepared follwoing the methodology described in Niethard et al.

%% initiate result str array to store all outputs
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
% We are looking for signal crossing from POSITIVE to NEGATIVE
% This marks the START of each SO cycle (beginning of downstate)

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
          pos_to_neg_crossings signal;


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

clearvars -except results EEG all_slow_waves spindles fs channel_name n_bouts signal

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

%% Step 9: add fields to results struct

% Store raw event arrays for metric extraction later
results.all_slow_waves    = all_slow_waves;
results.spindles          = spindles;
results.solitary_swa      = solitary_swa;
results.solitary_spindles = solitary_spindles;
results.coupled_events    = coupled_events;

% TF matrices (will be filled in Steps 13-15)
results.tf_sol_swa  = [];
results.tf_sol_swa_norm  = [];
results.tf_sol_sp   = [];
results.tf_sol_sp_norm   = [];
results.tf_coupled  = [];
results.tf_coupled_norm  = [];
results.time_axis   = [];
results.freqs       = [];

% Metrics (will be filled in Step 16+)
results.spindle_density = [];
results.pct_solitary    = [];
results.pct_coupled     = [];
results.dur_solitary    = [];
results.dur_coupled     = [];
results.power_solitary_full = [];
results.power_solitary_slow = [];
results.power_solitary_fast = [];
results.power_coupled_full  = [];
results.power_coupled_slow  = [];
results.power_coupled_fast  = [];
results.peak_freq_solitary = [];
results.peak_freq_coupled  = [];
results.spindle_phases     = [];
results.MRVL              = [];
results.MPA              = [];


%% Step 10: Time-Frequency analyses for solitary SWA

clearvars -except results EEG all_slow_waves spindles solitary_swa solitary_spindles ...
    coupled_events fs channel_name n_bouts signal results

% TF parameters
win_sec     = 2.0;
win_samples = round(win_sec * fs);
time_axis   = (-win_samples:win_samples) / fs;
freqs       = 1:0.5:25;
n_freqs = length(freqs);
% n_cycles    = 5;
n_cycles_vec = linspace(3, 10, n_freqs);  % 3 cycles at 1Hz, 10 cycles at 25Hz
n_freqs     = length(freqs);
n_time      = 2 * win_samples + 1;

% Collect epochs locked to SWA trough
epochs_sol_swa = [];
for bout = 1:n_bouts
    raw     = EEG.(channel_name).NREM_eeg{bout};
    sol_swa = solitary_swa{bout};
    for i = 1:size(sol_swa, 1)
        lock      = sol_swa(i, 3);
        seg_start = lock - win_samples;
        seg_end   = lock + win_samples;
        if seg_start >= 1 && seg_end <= length(raw)
            epochs_sol_swa(end+1, :) = raw(seg_start:seg_end);
        end
    end
end

% Compute TF
n_epochs = size(epochs_sol_swa, 1);
tf_sum   = zeros(n_freqs, n_time);
for ev = 1:n_epochs
    sig = epochs_sol_swa(ev, :);
    for fi = 1:n_freqs
        f       = freqs(fi);
        t_wav   = -2:1/fs:2;
        % sigma   = n_cycles / (2 * pi * f);
        sigma   = n_cycles_vec(fi) / (2 * pi * f);  % <-- uses per-frequency cycles
        wavelet = exp(2i*pi*f.*t_wav) .* exp(-t_wav.^2 / (2*sigma^2));
        wavelet = wavelet / sum(abs(wavelet));
        conv_result = conv(sig, wavelet, 'same');
        tf_sum(fi, :) = tf_sum(fi, :) + abs(conv_result).^2;
    end
end
tf_sol_swa = tf_sum / n_epochs;

% Baseline normalise (-2 to -1s)
baseline_start = dsearchn(time_axis', -2);
baseline_end   = dsearchn(time_axis', -1);
baseline_mean  = mean(tf_sol_swa(:, baseline_start:baseline_end), 2);
tf_sol_swa_norm        = tf_sol_swa ./ baseline_mean;


% Store TF plots
results.tf_sol_swa  = tf_sol_swa;
results.tf_sol_swa_norm  = tf_sol_swa_norm;
results.time_axis   = time_axis;
results.freqs       = freqs;

%% Step 11: Time-Frequency plots for solitary Spindles
clearvars -except results EEG all_slow_waves spindles solitary_swa solitary_spindles ...
    coupled_events fs channel_name n_bouts signal

% TF parameters
win_sec     = 2.0;
win_samples = round(win_sec * fs);
time_axis   = (-win_samples:win_samples) / fs;
freqs       = 1:0.5:25;
n_freqs = length(freqs);
% n_cycles    = 5;
n_cycles_vec = linspace(3, 10, n_freqs);
n_freqs     = length(freqs);
n_time      = 2 * win_samples + 1;

% Collect epochs
epochs_sol_sp = [];
for bout = 1:n_bouts
    raw    = EEG.(channel_name).NREM_eeg{bout};
    sol_sp = solitary_spindles{bout};
    for i = 1:size(sol_sp, 1)
        lock      = sol_sp(i, 1);
        seg_start = lock - win_samples;
        seg_end   = lock + win_samples;
        if seg_start >= 1 && seg_end <= length(raw)
            epochs_sol_sp(end+1, :) = raw(seg_start:seg_end);
        end
    end
end
fprintf('Solitary spindle epochs collected: %d\n', size(epochs_sol_sp, 1));

% Compute TF
n_epochs = size(epochs_sol_sp, 1);
tf_sum   = zeros(n_freqs, n_time);
for ev = 1:n_epochs
    sig = epochs_sol_sp(ev, :);
    for fi = 1:n_freqs
        f       = freqs(fi);
        t_wav   = -2:1/fs:2;
        % sigma   = n_cycles / (2 * pi * f);
        sigma   = n_cycles_vec(fi) / (2 * pi * f);  % <-- uses per-frequency cycles
        wavelet = exp(2i*pi*f.*t_wav) .* exp(-t_wav.^2 / (2*sigma^2));
        wavelet = wavelet / sum(abs(wavelet));
        conv_result = conv(sig, wavelet, 'same');
        tf_sum(fi, :) = tf_sum(fi, :) + abs(conv_result).^2;
    end
end
tf_sol_sp = tf_sum / n_epochs;

% Baseline normalise (-2 to -1s)
baseline_start = dsearchn(time_axis', -2);
baseline_end   = dsearchn(time_axis', -1);
baseline_mean  = mean(tf_sol_sp(:, baseline_start:baseline_end), 2);
tf_sol_sp_norm        = tf_sol_sp ./ baseline_mean;


% Store TF results before clearing
results.tf_sol_sp   = tf_sol_sp;
results.tf_sol_sp_norm  = tf_sol_sp_norm ;

%% Step 12: Coupled Events TF

clearvars -except results EEG all_slow_waves spindles solitary_swa solitary_spindles ...
    coupled_events fs channel_name n_bouts signal

% TF parameters
win_sec     = 2.0;
win_samples = round(win_sec * fs);
time_axis   = (-win_samples:win_samples) / fs;
freqs       = 1:0.5:25;
n_cycles_vec = linspace(3, 10, length(freqs));
n_freqs     = length(freqs);
n_time      = 2 * win_samples + 1;

% Collect epochs locked to SWA trough
epochs_coupled = [];
for bout = 1:n_bouts
    raw = EEG.(channel_name).NREM_eeg{bout};
    cp  = coupled_events{bout};
    for i = 1:size(cp, 1)
        swa_idx   = cp(i, 1);
        lock      = all_slow_waves{bout}(swa_idx, 3);  % trough
        seg_start = lock - win_samples;
        seg_end   = lock + win_samples;
        if seg_start >= 1 && seg_end <= length(raw)
            epochs_coupled(end+1, :) = raw(seg_start:seg_end);
        end
    end
end
% fprintf('Coupled epochs collected: %d\n', size(epochs_coupled, 1));

% Compute TF
n_epochs = size(epochs_coupled, 1);
tf_sum   = zeros(n_freqs, n_time);
for ev = 1:n_epochs
    sig = epochs_coupled(ev, :);
    for fi = 1:n_freqs
        f       = freqs(fi);
        t_wav   = -2:1/fs:2;
        sigma   = n_cycles_vec(fi) / (2 * pi * f);
        wavelet = exp(2i*pi*f.*t_wav) .* exp(-t_wav.^2 / (2*sigma^2));
        wavelet = wavelet / sum(abs(wavelet));
        conv_result = conv(sig, wavelet, 'same');
        tf_sum(fi, :) = tf_sum(fi, :) + abs(conv_result).^2;
    end
end
tf_coupled = tf_sum / n_epochs;

% Baseline normalise
baseline_start = dsearchn(time_axis', -2);
baseline_end   = dsearchn(time_axis', -1);
baseline_mean  = mean(tf_coupled(:, baseline_start:baseline_end), 2);
tf_coupled_norm        = tf_coupled ./ baseline_mean;


% Store TF results before clearing
results.tf_coupled  = tf_coupled;
results.tf_coupled_norm  = tf_coupled_norm;

%% Step 13: Spindle counts and density

clearvars -except results EEG all_slow_waves spindles solitary_swa solitary_spindles ...
    coupled_events fs channel_name n_bouts signal 


total_NREM_duration = 0;
for bout = 1:n_bouts
    total_NREM_duration = total_NREM_duration + length(EEG.(channel_name).NREM_eeg{bout}) / fs;
end
total_NREM_min = total_NREM_duration / 60;

n_solitary_sp = sum(cellfun(@(x) size(x,1), solitary_spindles));
n_coupled_sp  = sum(cellfun(@(x) size(x,1), coupled_events));
n_total_sp    = n_solitary_sp + n_coupled_sp;

spindle_density     = n_total_sp / total_NREM_min;
pct_solitary        = 100 * n_solitary_sp / n_total_sp;
pct_coupled         = 100 * n_coupled_sp  / n_total_sp;

% Store in results array
results.spindle_density = spindle_density; % density of all spindles
results.pct_solitary    = pct_solitary;
results.pct_coupled     = pct_coupled;

%% Step 14: Spindle duration

clearvars -except results EEG all_slow_waves spindles solitary_swa solitary_spindles ...
    coupled_events fs channel_name n_bouts signal 

sol_dur = [];
for bout = 1:n_bouts
    sp = solitary_spindles{bout};
    if ~isempty(sp)
        sol_dur = [sol_dur; sp(:, 3)];  % duration is column 3
    end
end

coup_dur = [];
for bout = 1:n_bouts
    cp = coupled_events{bout};
    for i = 1:size(cp, 1)
        sp_idx = cp(i, 2);
        coup_dur = [coup_dur; spindles{bout}(sp_idx, 3)];
    end
end

% Add to results array
if ~isempty(sol_dur),  results.dur_solitary = mean(sol_dur);  end
if ~isempty(coup_dur), results.dur_coupled  = mean(coup_dur); end


%% Step 15: Spindle power and peak frequency (wavelet relative power)
% core script wavelet_relative_cg 
    % Returns relative power (% of 1-30 Hz total at each time point)
    % Output:Power_relative: Relative power (non-log, as percentages)

clearvars -except results EEG all_slow_waves spindles solitary_swa solitary_spindles ...
    coupled_events fs channel_name n_bouts signal 

bands.full = [10 15];
bands.slow = [10 12];
bands.fast = [12 15];

sol_full = []; sol_slow = []; sol_fast = []; sol_peak = [];
coup_full = []; coup_slow = []; coup_fast = []; coup_peak = [];

for bout = 1:n_bouts
    raw = EEG.(channel_name).NREM_eeg{bout};
    
    for i = 1:size(solitary_spindles{bout}, 1)
        seg = raw(solitary_spindles{bout}(i,1):solitary_spindles{bout}(i,2));
        [Prel, ~, Freq] = Wavelet_CG_relative(seg);
        sol_full(end+1) = mean(Prel(Freq>=10  & Freq<=15, :), 'all');
        sol_slow(end+1) = mean(Prel(Freq>=10  & Freq<=12, :), 'all');
        sol_fast(end+1) = mean(Prel(Freq>=12 & Freq<=15, :), 'all');
        idx = Freq>=10 & Freq<=15;
        [~, mi] = max(mean(Prel(idx,:), 2));
        sol_peak(end+1) = Freq(find(idx,1) + mi - 1);
    end
    
    for i = 1:size(coupled_events{bout}, 1)
        sp_idx = coupled_events{bout}(i,2);
        seg = raw(spindles{bout}(sp_idx,1):spindles{bout}(sp_idx,2));
        [Prel, ~, Freq] = Wavelet_CG_relative(seg);
        coup_full(end+1) = mean(Prel(Freq>=10  & Freq<=15, :), 'all');
        coup_slow(end+1) = mean(Prel(Freq>=10  & Freq<=12, :), 'all');
        coup_fast(end+1) = mean(Prel(Freq>=12 & Freq<=15, :), 'all');
        idx = Freq>=10 & Freq<=15;
        [~, mi] = max(mean(Prel(idx,:), 2));
        coup_peak(end+1) = Freq(find(idx,1) + mi - 1);
    end
end

% add to results array
results.power_solitary_full = mean(sol_full);
results.power_solitary_slow = mean(sol_slow);
results.power_solitary_fast = mean(sol_fast);
results.power_coupled_full  = mean(coup_full);
results.power_coupled_slow  = mean(coup_slow);
results.power_coupled_fast  = mean(coup_fast);
results.peak_freq_solitary = mean(sol_peak);
results.peak_freq_coupled  = mean(coup_peak);

%% Step 16: MRL and MPA for coupled events
clearvars -except results coupled_events EEG channel_name all_slow_waves spindles

fs        = 1000;
n_bouts = length(EEG.(channel_name).NREM_eeg);

order = 2;
[b, a] = butter(order, [0.5 4]/(fs/2), 'bandpass');

spindle_phases = [];

for bout = 1:n_bouts
    raw      = EEG.(channel_name).NREM_eeg{bout}; % get raw EEG/LFP from the bout
    ph_filt  = filtfilt(b, a, raw); % filter at SWA range
    swa_phase = angle(hilbert(ph_filt));  % -pi to pi
    
    cp = coupled_events{bout};
    for i = 1:size(cp, 1)
        sp_idx   = cp(i, 2);
        sp_onset = spindles{bout}(sp_idx, 1);
        spindle_phases(end+1) = swa_phase(sp_onset);
    end
end

% MRL
MRL = abs(mean(exp(1i * spindle_phases)));
mean_angle = rad2deg(angle(mean(exp(1i * spindle_phases))));

results.MRVL = MRL;
results.MPA = mean_angle;
results.spindle_phases = spindle_phases;

end

%% end of script
