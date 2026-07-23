% =========================================================================
% T4TE Study - Resting State EEG Preprocessing 
% =========================================================================

% =========================================================================

clear; close all; clc

subj_code  = 
raw_file   = 
raw_path   = 
proc_path  = 

emg_chans  = {'APB', 'FDI', 'EDC'};
bad_chans  = {'Fpz', 'Oz', 'T7', 'T8'};

hjorth_required = {'C3', 'FC1', 'CP1', 'FC5', 'CP5'};
target_fs    = 500;
bandpass_hz  = [1 45];
epoch_len_s  = 2;

ft_path     = 
eeglab_path = 

addpath(ft_path); ft_defaults;
addpath(eeglab_path); eeglab nogui;
if ~exist(proc_path, 'dir'); mkdir(proc_path); end

log              = struct();
log.subj_code    = subj_code;
log.date         = datestr(now, 'yyyy-mm-dd HH:MM');
log.raw_file     = raw_file;
log.pipeline_ver = 'RS_v1';

fprintf('\n=== STEP 1: Loading raw data ===\n');
EEG = pop_loadbv(raw_path, [raw_file '.vhdr']);
EEG = eeg_checkset(EEG);
orig_fs = EEG.srate;
fprintf('Channels: %d | Rate: %d Hz | Duration: %.1f s\n', EEG.nbchan, EEG.srate, EEG.xmax);
log.orig_fs = orig_fs; log.duration_s_raw = EEG.xmax; log.n_chans_raw = EEG.nbchan;
data = eeglab2fieldtrip(EEG, 'preprocessing'); clear EEG

fprintf('\n=== STEP 2: Removing EMG channels (%s) ===\n', strjoin(emg_chans, ', '));
emg_present = intersect(emg_chans, data.label);
if ~isempty(emg_present)
    excl_emg = cellfun(@(c) ['-' c], emg_present, 'UniformOutput', false);
    cfg = []; cfg.channel = ft_channelselection(['all'; excl_emg(:)], data.label);
    data = ft_selectdata(cfg, data);
    fprintf('Removed: %s\n', strjoin(emg_present, ', '));
end
log.emg_chans_removed = emg_present;

fprintf('\n=== STEP 3: Removing bad channels ===\n');
fprintf('Removing: %s\n', strjoin(bad_chans, ', '));
excl = cellfun(@(c) ['-' c], bad_chans, 'UniformOutput', false);
cfg = []; cfg.channel = ft_channelselection(['all'; excl(:)], data.label);
data = ft_selectdata(cfg, data);
fprintf('Channels after removal: %d\n', numel(data.label));
log.bad_chans_removed = bad_chans; log.n_chans_after_removal = numel(data.label);
missing_hjorth = setdiff(hjorth_required, data.label);
if ~isempty(missing_hjorth)
    warning('Missing Hjorth channels: %s', strjoin(missing_hjorth, ', '));
else; fprintf('Hjorth channels: all present\n'); end

fprintf('\n=== STEP 4: Downsampling %d -> %d Hz ===\n', orig_fs, target_fs);
cfg = []; cfg.resamplefs = target_fs; cfg.detrend = 'no';
data = ft_resampledata(cfg, data);
log.target_fs = target_fs;

fprintf('\n=== STEP 5: Bandpass filter %d-%d Hz ===\n', bandpass_hz(1), bandpass_hz(2));
cfg = []; cfg.bpfilter = 'yes'; cfg.bpfreq = bandpass_hz;
cfg.bpfilttype = 'firws'; cfg.bpfiltdir = 'twopass';
data = ft_preprocessing(cfg, data);
log.bandpass_hz = bandpass_hz;

% OPTIONAL: 50 Hz notch — uncomment if line noise visible in PSD
% cfg = []; cfg.bsfilter = 'yes'; cfg.bsfreq = [49 51];
% cfg.bsfilttype = 'firws'; cfg.bsfiltdir = 'twopass';
% data = ft_preprocessing(cfg, data);

fprintf('\n=== STEP 6: Average re-reference ===\n');
cfg = []; cfg.reref = 'yes'; cfg.refchannel = 'all';
data = ft_preprocessing(cfg, data);
fprintf('Re-reference applied over %d channels.\n', numel(data.label));

fprintf('\n=== STEP 7: Epoching into %d s segments ===\n', epoch_len_s);
cfg = []; cfg.length = epoch_len_s; cfg.overlap = 0;
data_epoched = ft_redefinetrial(cfg, data); clear data
n_epochs_total = numel(data_epoched.trial);
fprintf('Epochs: %d\n', n_epochs_total);
log.epoch_len_s = epoch_len_s; log.n_epochs_total = n_epochs_total;

fprintf('\n=== STEP 8: Artifact rejection (visual summary) ===\n');
fprintf('Inspect summary plot. Click outlier trials to reject, then press quit.\n');
fprintf('Look for: amplitude spikes, muscle bursts, sustained drift\n\n');
cfg = []; cfg.method = 'summary'; cfg.keepchannel = 'yes'; cfg.metric = 'maxabs';
data_clean = ft_rejectvisual(cfg, data_epoched);
n_epochs_clean = numel(data_clean.trial);
n_epochs_rejected = n_epochs_total - n_epochs_clean;
pct_rejected = 100 * n_epochs_rejected / n_epochs_total;
fprintf('Rejected: %d (%.1f%%) | Retained: %d\n', n_epochs_rejected, pct_rejected, n_epochs_clean);
if pct_rejected > 30; warning('More than 30%% rejected.'); end
log.n_epochs_rejected = n_epochs_rejected; log.n_epochs_clean = n_epochs_clean;
log.pct_rejected = pct_rejected;

fprintf('\n=== STEP 9: Saving ===\n');
save(fullfile(proc_path, [subj_code '_RS_clean.mat']), 'data_clean', '-v7.3');
log.n_chans_final = numel(data_clean.label); log.chans_final = data_clean.label;
save(fullfile(proc_path, [subj_code '_RS_prepro_log.mat']), 'log');
fid = fopen(fullfile(proc_path, [subj_code '_RS_prepro_log.txt']), 'w');
fprintf(fid, 'Resting State Preprocessing Log — T4TE RS_v1\n');
fprintf(fid, 'Subject: %s | Date: %s | File: %s\n', log.subj_code, log.date, log.raw_file);
fprintf(fid, '-----------------------------------------\n');
fprintf(fid, 'INPUT   Channels: %d | Rate: %d Hz | Duration: %.1f s\n', log.n_chans_raw, log.orig_fs, log.duration_s_raw);
fprintf(fid, 'EMG removed:      %s\n', strjoin(log.emg_chans_removed, ', '));
fprintf(fid, 'CHANNELS removed: %s\n', strjoin(log.bad_chans_removed, ', '));
fprintf(fid, 'RESAMPLED to: %d Hz\n', log.target_fs);
fprintf(fid, 'BANDPASS: %d-%d Hz (firws, twopass)\n', log.bandpass_hz(1), log.bandpass_hz(2));
fprintf(fid, 'EPOCHS   Length: %d s | Total: %d\n', log.epoch_len_s, log.n_epochs_total);
fprintf(fid, 'REJECTED %d (%.1f%%) | Retained: %d\n', log.n_epochs_rejected, log.pct_rejected, log.n_epochs_clean);
fprintf(fid, 'FINAL    Channels: %d\n', log.n_chans_final);
fprintf(fid, '=========================================\n');
fclose(fid);
fprintf('Complete — BEL_S01\nNext: T4TE_IAF_estimation_BEL_S01.m\n');
