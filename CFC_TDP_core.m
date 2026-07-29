function [phaseamp] = CFC_TDP_core(eeg, phas_freqs, ampl_freqs)
%% Croos Frequency Coupling to measure how much observed modulation index is far away from the mean of null hypothesis 
% distribution of MI(boot strapped values of MI) - in the end we plot z-modulations. 
%% Load the signal (to test the script)
% clearvars -except EEG
% eeg = EEG.active_wake_two_sec_epochs.WT_PFC_LFP{1,1} (1,:); % 1 X 2000 ( 2000 time points one trial data)

srate = 1000;
npnts = length(eeg);

%  full-width at half maximum (FWHM) - the width of the peak at half the maximum amplitude

%% Calculate CFC across frequency range (generate "phaseamp" matrix) 

% number of iterations used for permutation testing
n_iter = 200; 

% initialize output phase-amplitude matrix
phaseamp = zeros(length(phas_freqs),length(ampl_freqs));

% loop over frequencies for phase
for lower_fi=1:length(phas_freqs)
    
    % get phase values
    phasefilt = filterFGx(eeg,srate,phas_freqs(lower_fi),phas_freqs(lower_fi)*.25);
    phase = angle(hilbert(phasefilt));
    
    for upper_fi=1:length(ampl_freqs)
        
        % get power (or amp) values 
        ampfilt = filterFGx(eeg,srate,ampl_freqs(upper_fi),ampl_freqs(upper_fi)*.2);
        amplit = abs(hilbert(ampfilt)).^2;
        
        % calculate observed modulation index
        modidx = abs(mean(amplit.*exp(1i*phase)));

        % use permutation testing to get Z-value
        bm = zeros(1,n_iter);
        for bi=1:n_iter
            cutpoint = randsample(round(npnts/10):round(npnts*.9),1);
            bm(bi) = abs(mean(amplit([ cutpoint:end 1:cutpoint-1 ]).*exp(1i*phase)));
        end

        % the value we use is the normalized distance away from the mean of
        % boot-strapped values
        phaseamp(lower_fi,upper_fi) = (modidx-mean(bm))/std(bm);
    end % end upper frequency loop (for amplitude)
end % end lower frequency loop (for phase)

