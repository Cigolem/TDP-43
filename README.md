# TDP-43
All codes used for the manuscript titled "Disrupted sleep oscillatory architecture in a TDP-43 Q331K knock-in mouse model of ALS-FTD"


Analysis of sleep architecture, spectral, aperiodic, cross-frequency coupling, and spindle measures from ECoG and LFP recordings in wild-type (WT) and Q331K knock-in mice.

## Figure scripts

| File | Figure | Language | Description |
|------|--------|----------|-------------|
| `fig_1_perc_stats.ipynb` | Fig 1 | Python | Hourly % time in Wake/NREM/REM across 24 h; circadian statistics |
| `fig_2_sleep_arch_plot_stats.ipynb` | Fig 2 | Python | Sleep architecture (transitions, bout duration, total time); mixed ANOVA |
| `figure_3_morlet_plots.m` | Fig 3 | MATLAB | Time–frequency spectrograms of NREM–REM transitions |
| `figure_3_REM_violin_plots.m` | Fig 3 | MATLAB | REM band-power violin plots (theta, beta) |
| `figure_4_REM_CFC.m` | Fig 4 | MATLAB | REM theta–gamma phase–amplitude coupling |
| `figure_5_spindle.m` | Fig 5 | MATLAB | Spindle density and duration |
| `figure_6_spindle.m` | Fig 6 | MATLAB | Event-locked spindle spectrograms and band power |
| `figure_7_spindle.m` | Fig 7 | MATLAB | Spindle–slow wave coupling phase |
| `figure_8_aperiodic.ipynb` | Fig 8 | Python | Aperiodic (FOOOF/specparam) exponent and offset |
| `suppl_fig_1_spindle.m` | Suppl. Fig 1 | MATLAB | Representative coupled slow wave–spindle event |

## Core and helper functions

| File | Description |
|------|-------------|
| `detect_swa_sp_core_5.m` | Slow-wave and spindle detection, event classification (solitary/coupled), time–frequency analysis, density, duration, band power, and coupling phase (MRVL/MPA). Called by the spindle figure scripts. |
| `Wavelet_CG_relative.m` | Morlet wavelet transform returning relative power (% of total 1–30 Hz power at each time point). |
| `CFC_TDP_core.m` | Theta–gamma modulation index with surrogate z-scoring, for the PAC analysis. |
| `simple_violin.m` | Violin-plot helper (median and interquartile range). |
| `wavelet.m` | Continuous wavelet transform (Torrence & Compo, 1998). |

## Detection summary (`detect_swa_sp_core_5.m`)

- Slow waves: 0.5–4 Hz zero-phase Butterworth; complete cycles between positive-to-negative zero-crossings, 0.4–2.0 s; negative-peak and peak-to-peak amplitudes ≥ 66.6% of their respective means.
- Spindles: 10–15 Hz zero-phase Butterworth; smoothed Hilbert envelope; dual threshold (85th percentile, with ≥1 sample above the 97th); events < 200 ms apart merged; retained if 0.5–3.0 s.
- Coupling: a spindle is coupled if its onset falls within the slow-wave trough-to-peak (upstate) window; otherwise solitary.
- Phase convention: 0° = positive peak, ±180° = trough.

## Requirements

- **Python**: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `fooof` (specparam)
- **MATLAB** (R2019b or later): Signal Processing Toolbox

## Notes on frequency ranges

Ranges differ by analysis and are described in the Methods / Supplementary Methods: Morlet time–frequency and spindle relative power are normalised to total 1–30 Hz power; PSD band power is normalised to total 1–45 Hz power; aperiodic fitting (FOOOF) is performed over 5–45 Hz.

## Data availability

Raw recordings are available from the corresponding author on reasonable request.

## Citation

Gelegen van Eijl C, et al. Disrupted sleep oscillatory architecture in a TDP-43 Q331K knock-in mouse model of ALS-FTD. *[Journal, year — to be updated on publication].*

## License

Released under the MIT License (see `LICENSE`).
