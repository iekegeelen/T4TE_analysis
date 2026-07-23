% =========================================================================
% T4TE Study 1.1 — Pre-stimulus μ-alpha Power Estimation
% =========================================================================
% Computes pre-stimulus μ-alpha band power per trial, matched exactly to
% the valid suprathreshold trials used in T4TE_phase_MEP_analysis_v7.m.
% Trial selection logic is replicated from v7 to ensure perfect alignment
% with the phase and MEP data already saved in the results struct.
%
% Input:
%   BEL_S0X_eeg_final.set   — preprocessed EEG (EEGLAB)
%   BEL_S0X_MEP.mat         — MEP amplitudes + eeg_mask
%   BEL_S0X_emg.mat         — data_emg.trialinfo
%   T4TE_phase_MEP_results_v7_IAF.mat — alpha_band per subject
%
% Output: BEL_S0X_prestim_power.mat per subject, containing:
%   power_IAF       [n_valid x 1]  raw power in IAF band (µV²)
%   power_fixed     [n_valid x 1]  raw power in fixed 8-13 Hz band (µV²)
%   power_IAF_z     [n_valid x 1]  z-scored within subject
%   power_fixed_z   [n_valid x 1]  z-scored within subject
%   phase           [n_valid x 1]  phase angles (copied from v7 results)
%   mep_z           [n_valid x 1]  z-MEP (copied from v7 results)
%   IAF_used        scalar         IAF (Hz)
%   alpha_band      [1x2]          filter band used
%   n_valid         scalar
%   subj            string
%
% Author:  E.W.M. Dresens
% Date:    June 2026
% =========================================================================

clear; close all; clc

% =========================================================================
%% SETTINGS — change subj per run
% =========================================================================
subj      = ;
base_path = ;

% Hjorth-Laplacian (identical to v7)
eeg_channel       = 'C3';
hjorth_neighbours = {'FC1','CP1','FC5','CP5'};

% Pre-stimulus power window (ms) — same epoch used by Phastimate
% Use -1000 to -50 ms to avoid edge effects at both ends
prestim_win = [-1000, -50];

% Fixed band
fixed_band = [8, 13];

% =========================================================================
%% TOOLBOXES
% =========================================================================
ft_path     = '';
eeglab_path = '';
addpath(ft_path); ft_defaults;
addpath(eeglab_path); eeglab nogui;

fprintf('=== T4TE Pre-stimulus Power: %s ===\n', subj);

% =========================================================================
%% LOAD ALPHA BAND FROM V7 RESULTS
% =========================================================================
results_file = fullfile(base_path, 'T4TE_phase_MEP_results_v7_IAF.mat');
tmp = load(results_file, 'results');
results_all = tmp.results;

subj_idx = find(strcmpi({results_all.subject}, subj));
if isempty(subj_idx)
    error('Subject %s not found in results file.', subj);
end

IAF_used   = results_all(subj_idx).IAF_used;
alpha_band = results_all(subj_idx).alpha_band;
phase_v7   = results_all(subj_idx).phase(:);   % [n_valid x 1]
mep_z_v7   = results_all(subj_idx).mep_z(:);

fprintf('IAF used: %.2f Hz | Band: [%.2f %.2f Hz]\n', ...
    IAF_used, alpha_band(1), alpha_band(2));
fprintf('N valid trials from v7: %d\n', length(phase_v7));

% =========================================================================
%% LOAD EEG, MEP, EMG  (replicate v7 trial selection)
% =========================================================================
proc_path = fullfile(base_path, subj, 'processed');

EEG = pop_loadset('filename', [subj '_eeg_final.set'], 'filepath', proc_path);
load(fullfile(proc_path, [subj '_MEP.mat']), 'MEP', 'eeg_mask');
load(fullfile(proc_path, [subj '_emg.mat']), 'data_emg');

fs_eeg   = EEG.srate;
time_eeg = EEG.times / 1000;   % ms -> s

% =========================================================================
%% REPLICATE V7 TRIAL SELECTION
% =========================================================================
% Step 1: get condition per EEG epoch
cond_eeg = zeros(1, EEG.trials);
for ep = 1:EEG.trials
    evtypes = [EEG.epoch(ep).eventtype];
    if iscell(evtypes); evtypes = evtypes{1}; end
    if isnumeric(evtypes)
        cond_eeg(ep) = evtypes;
    else
        cond_eeg(ep) = str2double(evtypes);
    end
end

supra_eeg_idx  = find(cond_eeg == 3);
n_supra_eeg    = numel(supra_eeg_idx);

supra_orig_idx = find(data_emg.trialinfo == 3);
eeg_mask_supra = eeg_mask(supra_orig_idx);
mep_eeg_valid  = MEP.suprathreshold(eeg_mask_supra);
n_eeg_valid    = sum(eeg_mask_supra);

% Handle off-by-one (same as v7)
if n_eeg_valid == n_supra_eeg + 1
    mep_eeg_valid = mep_eeg_valid(1:n_supra_eeg);
    warning('%s: trimming 1 trial.', subj);
elseif n_eeg_valid ~= n_supra_eeg
    error('%s: trial count mismatch (%d vs %d).', subj, n_eeg_valid, n_supra_eeg);
end

% Step 2: remove NaN MEP trials
mep_valid_mask2 = ~isnan(mep_eeg_valid);
valid_ep_idx    = supra_eeg_idx(mep_valid_mask2);
n_valid_matched = numel(valid_ep_idx);

fprintf('Valid suprathreshold epochs (replicated): %d\n', n_valid_matched);

% Sanity check against v7
if n_valid_matched ~= length(phase_v7)
    warning(['Trial count mismatch: replicated=%d, v7=%d. ' ...
        'Power may not align perfectly with phase/MEP from v7.'], ...
        n_valid_matched, length(phase_v7));
end

% =========================================================================
%% HJORTH-LAPLACIAN SPATIAL FILTER
% =========================================================================
c3_idx = find(strcmpi({EEG.chanlocs.labels}, eeg_channel));
hjorth_idx = zeros(1, 4);
for hh = 1:4
    idx_h = find(strcmpi({EEG.chanlocs.labels}, hjorth_neighbours{hh}));
    if isempty(idx_h)
        warning('Hjorth neighbour %s not found, using C3.', hjorth_neighbours{hh});
        hjorth_idx(hh) = c3_idx;
    else
        hjorth_idx(hh) = idx_h;
    end
end

% =========================================================================
%% PRE-STIMULUS WINDOW INDICES
% =========================================================================
t_start = find(time_eeg >= prestim_win(1)/1000, 1, 'first');
t_end   = find(time_eeg <= prestim_win(2)/1000, 1, 'last');
n_samp  = t_end - t_start + 1;
fprintf('Power window: %d to %d ms (%d samples)\n', ...
    prestim_win(1), prestim_win(2), n_samp);

% =========================================================================
%% BANDPASS FILTERS
% =========================================================================
fir_order = 128;
b_IAF     = fir1(fir_order, alpha_band / (fs_eeg/2), 'bandpass');
b_fixed   = fir1(fir_order, fixed_band  / (fs_eeg/2), 'bandpass');

% =========================================================================
%% COMPUTE POWER PER VALID TRIAL
% =========================================================================
n_valid     = numel(valid_ep_idx);
power_IAF   = zeros(n_valid, 1);
power_fixed = zeros(n_valid, 1);

for tt = 1:n_valid
    ep  = valid_ep_idx(tt);
    c3  = double(squeeze(EEG.data(c3_idx, :, ep)));
    nbr = double(squeeze(mean(EEG.data(hjorth_idx, :, ep), 1)));
    hjorth_sig = c3 - 0.25 * nbr;   % same weighting as v7

    % Bandpass filter full trial, then extract prestim window
    seg_IAF   = filtfilt(b_IAF,   1, hjorth_sig);
    seg_fixed = filtfilt(b_fixed, 1, hjorth_sig);

    power_IAF(tt)   = mean(seg_IAF(t_start:t_end)   .^ 2);
    power_fixed(tt) = mean(seg_fixed(t_start:t_end)  .^ 2);
end

% Z-score within subject
power_IAF_z   = (power_IAF   - mean(power_IAF))   / std(power_IAF);
power_fixed_z = (power_fixed - mean(power_fixed)) / std(power_fixed);

fprintf('Power computed. IAF: mean=%.4f µV², SD=%.4f\n', ...
    mean(power_IAF), std(power_IAF));
fprintf('Power computed. Fixed: mean=%.4f µV², SD=%.4f\n', ...
    mean(power_fixed), std(power_fixed));

% =========================================================================
%% SAVE
% =========================================================================
out_file = fullfile(proc_path, sprintf('%s_prestim_power.mat', subj));

save(out_file, ...
    'power_IAF', 'power_IAF_z', ...
    'power_fixed', 'power_fixed_z', ...
    'phase_v7', 'mep_z_v7', ...
    'IAF_used', 'alpha_band', 'fixed_band', ...
    'prestim_win', 'n_valid', 'subj', 'fs_eeg');

fprintf('Saved: %s\n', out_file);
fprintf('=== Done: %s ===\n\n', subj);
