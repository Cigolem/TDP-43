

function [Power_relative, Power_relative_log10, Freq] = Wavelet_CG_relative(signal)
%% Wavelet analysis with INSTANTANEOUS RELATIVE POWER normalization
% Returns relative power (% of 1-30 Hz total at each time point)
%
% Outputs:
%   Power_relative: Relative power (non-log, as percentages)
%   Power_relative_log10: log10 of relative power (for visualization)

%% Wavelet Parameters
Fs = 1000;
dt = 1/Fs;
pad = 1;
maxfrequency = 30;
numoctaves = 6;
numsuboctaves = 4;
dj = 1/numsuboctaves;
s0 = 1/maxfrequency;
j1 = numoctaves/dj;
mother = 'Morlet';

%% Detrend and normalize signal
signal1 = detrend(signal);
signal2 = signal1 - nanmean(signal1);
n = length(signal2);
time = [0:length(signal2)] * dt;
time1 = time(1, 2:end);

%% Compute Morlet wavelet transform
[Power, Freq] = wavelet(signal2, dt, pad, dj, s0, j1, mother);
Power = flipud((abs(Power)).^2);  % Absolute power
Freq = fliplr(1./Freq);

%% Calculate instantaneous relative power
% Define frequency range for normalization (1-30 Hz)
freq_range = Freq >= 1 & Freq <= 30;

% Calculate total power (1-30 Hz) at EACH time point
total_power_time = sum(Power(freq_range, :), 1);  % Sum across frequencies, for each time point

% Calculate relative power: (power at each freq) / (total power at that time) * 100
Power_relative = zeros(size(Power));
for t = 1:size(Power, 2)
    if total_power_time(t) > 0  % Avoid division by zero
        Power_relative(:, t) = (Power(:, t) / total_power_time(t)) * 100;
    end
end

%% Log transform for visualization
% Add small constant to avoid log(0)
Power_relative_log10 = log10(Power_relative + eps);


end

