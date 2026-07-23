% =========================================================================
% T4TE Study - Full Preprocessing Pipeline v9
% =========================================================================
%
% Pipeline:
%   FIELDTRIP:
%   1.  Load raw data (lazy) + visualization
%   2.  Epoch (-1500 to +500 ms) + visualization
%   3.  Separate EMG channels -> save as .mat for separate MEP script
%   4.  Downsample 5000 -> 500 Hz (naive decimation, no anti-aliasing)
%   4b. NaN mark TMS artifact window (-5 to +10 ms) BEFORE visual inspection
%   5.  Visual inspection (on NaN-marked data) + bad trial/channel removal
%       (blocksize=2s = 1 trial per page; both endpoints of each drag
%       selection used for trial mapping; preallocated concatenation)
%   6.  ICA (runica extended Infomax, EEGLAB path conflict resolved)
%       + visualization
%   6b. ICA component rejection + visualization
%   7.  ft_interpolatenan (pchip, 50 ms context) — FieldTrip native
%
%   EEGLAB / TESA:
%   8.  Convert to EEGLAB
%   9.  Channel locations (3D) + condition labels
%   10. SOUND (lambda=0.1, iter=10)
%   11. High-pass filter 0.5 Hz
%   12. Baseline correction (-1500 to -10 ms)
%   13. SSP-SIR: diagnostic check first, then removal  <-- BEFORE notch
%   14. Notch filter: 50, 100, 150, 200 Hz             <-- AFTER SSP-SIR
%   15. Bad channel interpolation + average reference
%   16. TEP visualization
%   17. Save dataset + preprocessing log
% =========================================================================

clear; close all; clc

% =========================================================================
% SETTINGS — edit per participant
% =========================================================================

subj_code  = 
raw_file   = 
raw_path   = 
proc_path  = 

trig_sub   = 'S  1';   % 90%  rMT (subthreshold)
trig_supra = 'S  3';   % 110% rMT (suprathreshold)

prestim    = 1.5;      % seconds before pulse (-1500 ms for alpha power extraction)
poststim   = 0.5;      % seconds after pulse  (+500 ms covers all TEP components)
tms_win    = [-5 10];  % TMS artifact window to remove (ms)
orig_fs    = 5000;
final_fs   = 500;

emg_chans  = {'APB', 'FDI', 'EDC'};

% Fill in after Step 5 visual inspection, then rerun:
%   bad_chans        = all channels that are bad (flat, noisy, drifting)
%   chans_interpolate = central subset of bad_chans to KEEP through ICA
%                       and interpolate later (Step 15)
%                       peripheral bad channels not listed here are removed
% Example: bad_chans = {'POz','Iz','C5'}; chans_interpolate = {'C5'};
bad_chans         = {'Fpz', 'Oz', 'T7', 'T8'};
chans_interpolate = {};

% Fill in after Step 6 visualization, then run Step 6b:
% Look for: eye blinks (frontal), eye movements (lateral frontal),
%           heartbeat (diffuse ~1 Hz), TMS muscle artifact (focal, early)
% Example: ica_reject = [1 3 7];
ica_reject = [1 2 3 9 18];

% Fill in after Step 13 diagnostic, then rerun Step 13:
% Example: sspsir_pc = 3;
sspsir_pc  = 4;

plot_chan  = 'C3';   % channel for intermediate visualizations

% =========================================================================
% SETUP
% =========================================================================

addpath('/Users/e.w.m.dresens/Documents/MATLAB/fieldtrip-20250106/');
ft_defaults;
addpath('/Users/e.w.m.dresens/Documents/MATLAB/eeglab2026.0.0/');
eeglab;

if ~exist(proc_path, 'dir'); mkdir(proc_path); end
full_path = fullfile(raw_path, raw_file);

% Preprocessing log — saved as .mat and .txt at Step 17
log              = struct();
log.subj_code    = subj_code;
log.date         = datestr(now, 'yyyy-mm-dd HH:MM');
log.raw_file     = raw_file;
log.pipeline_ver = 'v9';

% =========================================================================
% PHASE 1: FIELDTRIP
% =========================================================================

%% STEP 1: Load raw data + visualization
% Only header and events are loaded into RAM; ft_databrowser streams from disk
fprintf('\n=== STEP 1: Loading raw data ===\n');

hdr   = ft_read_header(full_path);
event = ft_read_event(full_path);
fprintf('Rate: %d Hz | Channels: %d\n', hdr.Fs, hdr.nChans);

stim_events = event(strcmp({event.type}, 'Stimulus'));
vals        = {stim_events.value};
n_sub       = sum(strcmp(vals, trig_sub));
n_supra     = sum(strcmp(vals, trig_supra));
n_orig      = n_sub + n_supra;
fprintf('Subthreshold: %d | Suprathreshold: %d\n', n_sub, n_supra);

log.n_trials_input_sub   = n_sub;
log.n_trials_input_supra = n_supra;
log.n_trials_input_total = n_orig;

cfg_browser                = [];
cfg_browser.dataset        = full_path;   % streams from disk — no RAM load
cfg_browser.preproc.demean = 'yes';
ft_databrowser(cfg_browser);

%% STEP 2: Epoch + visualization

% trl matrix (col 4 = trigger value) is preserved for condition tracking
fprintf('\n=== STEP 2: Epoching ===\n');
cfg                     = [];
cfg.dataset             = full_path;
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = {trig_sub, trig_supra};
cfg.trialdef.prestim    = prestim;
cfg.trialdef.poststim   = poststim;
cfg                     = ft_definetrial(cfg);
trl                     = cfg.trl;

% -------------------------------------------------------------------------
% BEL_S02 ONLY: remove 11 aborted suprathreshold trials from block 1.
% Delivered before TMS trigger cable was replaced — t=0 unreliable.
% Identified as the first 11 rows in trl where column 4 == 3 (S  3).
% -------------------------------------------------------------------------
if strcmp(subj_code, 'BEL_S02')
    supra_rows = find(trl(:, 4) == 3);
    trl(supra_rows(1:11), :) = [];
    cfg.trl = trl;
    fprintf('BEL_S02: removed 11 aborted suprathreshold trials.\n');
    fprintf('Remaining: %d sub | %d supra\n', ...
        sum(trl(:,4) == 1), sum(trl(:,4) == 3));
end
% -------------------------------------------------------------------------

cfg.reref               = 'no';
cfg.demean              = 'no';
data_epoched            = ft_preprocessing(cfg);
fprintf('%d trials at %d Hz\n', length(data_epoched.trial), data_epoched.fsample);

cfg_browser            = [];
cfg_browser.preproc.demean = 'yes';
ft_databrowser(cfg_browser, data_epoched);

cfg_avg                = [];
cfg_avg.preproc.demean = 'yes';
data_avg               = ft_timelockanalysis(cfg_avg, data_epoched);
ft_databrowser(cfg_avg, data_avg);


%% STEP 3: Separate EMG channels
% EMG (APB, FDI, EDC) retained at 5000 Hz and saved for separate MEP analysis script
fprintf('\n=== STEP 3: Separating EMG channels ===\n');

emg_idx = find(ismember(data_epoched.label, emg_chans));
eeg_idx = find(~ismember(data_epoched.label, emg_chans));
if isempty(emg_idx)
    warning('No EMG channels found: %s', strjoin(emg_chans, ', '));
end

data_emg       = data_epoched;
data_emg.label = data_epoched.label(emg_idx);
for i = 1:length(data_epoched.trial)
    data_emg.trial{i} = data_epoched.trial{i}(emg_idx, :);
end
save(fullfile(proc_path, [subj_code '_emg.mat']), 'data_emg');
fprintf('EMG saved at %d Hz.\n', orig_fs);

data_eeg       = data_epoched;
data_eeg.label = data_epoched.label(eeg_idx);
for i = 1:length(data_epoched.trial)
    data_eeg.trial{i} = data_epoched.trial{i}(eeg_idx, :);
end
fprintf('EEG channels remaining: %d\n', length(data_eeg.label));
clear data_epoched data_emg

%% STEP 4: Downsample 5000 -> 500 Hz
fprintf('\n=== STEP 4: Downsampling %d -> %d Hz ===\n', orig_fs, final_fs);

factor          = orig_fs / final_fs;
data_ds         = data_eeg;
data_ds.fsample = final_fs;

% Update trl sample indices to match downsampled rate
trl_ds      = trl;
trl_ds(:,1) = ceil(trl(:,1) / factor);
trl_ds(:,2) = ceil(trl(:,2) / factor);
trl_ds(:,3) = ceil(trl(:,3) / factor);

for i = 1:length(data_eeg.trial)
    data_ds.trial{i} = data_eeg.trial{i}(:, 1:factor:end);
    data_ds.time{i}  = data_eeg.time{i}(1:factor:end);
end
clear data_eeg

%% STEP 4b: NaN mark TMS artifact window — BEFORE visual inspection
% Replacing the artifact window with NaN before visual inspection removes
% the large TMS spike from view, making it much easier to identify true
% bad trials (muscle bursts, eye movements, drift) and bad channels.
% NaN values are excluded from runica decomposition automatically by
% FieldTrip's runica wrapper, and filled via ft_interpolatenan after ICA.
fprintf('\n=== STEP 4b: NaN marking TMS artifact window (%d to %d ms) ===\n', ...
    tms_win(1), tms_win(2));

t_vec     = data_ds.time{1};
t0_samp   = find(abs(t_vec) < 1e-9);
assert(isscalar(t0_samp), 'Expected exactly one sample at t=0. Check epoch alignment.');
pre_samp  = round(abs(tms_win(1)) * final_fs / 1000);
post_samp = round(tms_win(2)      * final_fs / 1000);
win_start = t0_samp - pre_samp;
win_end   = t0_samp + post_samp;

data_nan = data_ds;
for itrial = 1:length(data_ds.trial)
    data_nan.trial{itrial}(:, win_start:win_end) = NaN;
end
fprintf('NaN marked: samples %d to %d (%.0f to %.0f ms)\n', ...
    win_start, win_end, t_vec(win_start)*1000, t_vec(win_end)*1000);
clear data_ds

%% STEP 5: Visual inspection + bad trial/channel removal
% Data shown WITHOUT the TMS artifact (NaN window) for clearer inspection.
% blocksize = 2 s = exactly 1 trial per page (epoch = -1500 to +500 ms).
% Both start and end sample of each drag selection are used to determine
% which trials to reject, so full-page drags correctly flag both trials.
% Preallocated concatenation used for speed.
% Right-click and drag to mark bad segments. Close window when done.
fprintf('\n=== STEP 5: Visual inspection (TMS artifact removed from view) ===\n');
fprintf('Right-click and drag to mark bad segments. Close window when done.\n');
fprintf('Look for: pre-stimulus muscle bursts, eye movements, jumps, drift\n\n');

nTrials  = length(data_nan.trial);
nSamples = size(data_nan.trial{1}, 2);
n_orig = length(data_nan.trial);

% Preallocate — avoids repeated memory reallocation in growing cat loop
data_pool = zeros(length(data_nan.label), nTrials * nSamples);
for tr = 1:nTrials
    col_start = (tr-1) * nSamples + 1;
    col_end   = tr * nSamples;
    data_pool(:, col_start:col_end) = data_nan.trial{tr};
end

data_temp.label   = data_nan.label;
data_temp.fsample = data_nan.fsample;
data_temp.trial   = {data_pool};
data_temp.time    = {(0:size(data_pool,2)-1) / data_nan.fsample};

cfg_browser                  = [];
cfg_browser.viewmode         = 'vertical';
cfg_browser.preproc.demean   = 'yes';
cfg_browser.preproc.hpfilter = 'yes';
cfg_browser.preproc.hpfreq   = 0.5;
cfg_browser.continuous       = 'yes';
cfg_browser.blocksize        = 2;   % 2 s = 1 trial per page
data_temp = ft_databrowser(cfg_browser, data_temp);

% Extract rejected trial indices from artifact markings
% Both start and end sample used — full-page drags flag both trials
if isfield(data_temp, 'artfctdef') && ...
   isfield(data_temp.artfctdef, 'visual') && ...
   ~isempty(data_temp.artfctdef.visual.artifact)
    trl2rej = [];
    for a = 1:size(data_temp.artfctdef.visual.artifact, 1)
        t_start = ceil(data_temp.artfctdef.visual.artifact(a,1) / nSamples);
        t_end   = ceil(data_temp.artfctdef.visual.artifact(a,2) / nSamples);
        trl2rej = [trl2rej, t_start:t_end];
    end
    trl2rej = unique(trl2rej);
    % clamp to valid trial range (prevents phantom index from drag past end)
trl2rej = trl2rej(trl2rej >= 1 & trl2rej <= size(trl_ds, 1));
fprintf('Trials marked for rejection: %d\n', length(trl2rej));

else
    trl2rej = [];
end
fprintf('Trials marked for rejection: %d\n', length(trl2rej));
clear data_temp data_pool

% Per-condition rejection breakdown (trl_ds still has all rows at this point)
trig_vals   = trl_ds(trl2rej, 4);
n_rej_sub   = sum(trig_vals == 1);
n_rej_supra = sum(trig_vals == 3);
fprintf('Rejected subthreshold: %d | Rejected suprathreshold: %d\n', n_rej_sub, n_rej_supra);

% Remove bad trials from both data_nan and trl_ds
data_clean = data_nan;
if ~isempty(trl2rej)
    trl2rej = sort(trl2rej, 'descend');
    for tr = trl2rej
        if tr <= length(data_clean.trial)
            data_clean.trial(tr) = [];
            data_clean.time(tr)  = [];
            if isfield(data_clean, 'trialinfo');  data_clean.trialinfo(tr,:)  = []; end
            if isfield(data_clean, 'sampleinfo'); data_clean.sampleinfo(tr,:) = []; end
            trl_ds(tr,:) = [];
        end
    end
end
fprintf('Trials remaining: %d\n', length(data_clean.trial));

% Check SOP rejection threshold
pct_visual = (1 - length(data_clean.trial) / n_orig) * 100;
fprintf('Visual rejection rate: %.1f%%\n', pct_visual);
if pct_visual > 30
    warning('More than 30%% rejected — check SOP exclusion criteria.');
end

% Remove bad channels
% NOTE: chans_interpolate are kept through ICA and interpolated at Step 15
if ~isempty(bad_chans)
    chans_remove = setdiff(bad_chans, chans_interpolate);
    if ~isempty(chans_remove)
        cfg_sel         = [];
        cfg_sel.channel = setdiff(data_clean.label, chans_remove);
        data_clean      = ft_selectdata(cfg_sel, data_clean);
        fprintf('Channels removed: %s\n', strjoin(chans_remove, ', '));
    end
    if ~isempty(chans_interpolate)
        fprintf('Channels kept for interpolation: %s\n', strjoin(chans_interpolate, ', '));
    end
end
fprintf('Channels remaining: %d\n', length(data_clean.label));

log.trials_rejected_visual    = trl2rej;
log.n_trials_after_visual     = length(data_clean.trial);
log.n_trials_rej_visual_sub   = n_rej_sub;
log.n_trials_rej_visual_supra = n_rej_supra;

% v9: explicitly name variable to guarantee strjoin() compatibility in Step 17
chans_removed            = setdiff(bad_chans, chans_interpolate);
log.bad_chans_removed    = chans_removed;
log.bad_chans_for_interp = chans_interpolate;

% Save post-visual checkpoint
vis_file = fullfile(proc_path, [subj_code '_postvis.mat']);
save(vis_file, 'data_clean', 'trl_ds', 'log', 'n_orig', '-v7.3');
fprintf('Checkpoint saved: %s\n', vis_file);



%% STEP 6: ICA (runica extended Infomax) + visualization
% runica extended is preferred over FastICA: more stable, models both
% sub- and super-Gaussian sources (TESA recommendation for TMS-EEG).
% ICA runs on NaN-marked data — artifact window excluded from decomposition.
% EEGLAB's sigprocfunc is temporarily removed to avoid runica path conflict.
fprintf('\n=== STEP 6: ICA (runica extended) — this will take 1-2 hours ===\n');

if ~exist('data_clean', 'var')
    vis_file = fullfile(proc_path, [subj_code '_postvis.mat']);
    if exist(vis_file, 'file')
        fprintf('Loading post-visual checkpoint: %s\n', vis_file);
        load(vis_file);   % restores data_clean, trl_ds, log, n_orig
        % v9: guard — reconstruct bad_chans log fields if missing from older checkpoint
        if ~isfield(log, 'bad_chans_removed')
            log.bad_chans_removed    = setdiff(bad_chans, chans_interpolate);
            log.bad_chans_for_interp = chans_interpolate;
            fprintf('Warning: bad_chans fields missing from checkpoint — reconstructed from SETTINGS.\n');
        end
    else
        error('data_clean not in workspace and no checkpoint found.\nExpected: %s\nRun Steps 1-5 first.', vis_file);
    end
end

% Remove EEGLAB's runica to avoid conflict with FieldTrip's version
rmpath(fullfile('/Users/e.w.m.dresens/Documents/MATLAB/eeglab2026.0.0/functions/sigprocfunc'));

cfg_ica                 = [];
cfg_ica.method          = 'runica';
cfg_ica.channel         = 'all';
cfg_ica.runica.extended = 1;
comp                    = ft_componentanalysis(cfg_ica, data_clean);

% Restore EEGLAB path for subsequent steps
addpath(fullfile('/Users/e.w.m.dresens/Documents/MATLAB/eeglab2026.0.0/functions/sigprocfunc'));

% Topographies + time courses for component identification
figure('Name', 'ICA Topographies');
cfg_plot           = [];
cfg_plot.component = 1:min(64, length(comp.label));
cfg_plot.layout    = 'EEG1010.lay';   % 2D FieldTrip layout, matched by channel label
cfg_plot.comment   = 'no';
ft_topoplotIC(cfg_plot, comp);

cfg_browser2          = [];
cfg_browser2.viewmode = 'component';
cfg_browser2.layout   = 'EEG1010.lay';
ft_databrowser(cfg_browser2, comp);

save('comp_ICA.mat', 'comp', '-v7.3');

fprintf('\n>>> Fill in ica_reject in SETTINGS, then run Step 6b <<<\n');

%% power spectrum 

% Tijdelijk NaN -> 0 voor power spectrum berekening
comp_freq = comp;
for tr = 1:length(comp_freq.trial)
    nan_idx = isnan(comp_freq.trial{tr});
    comp_freq.trial{tr}(nan_idx) = 0;
end

% Power spectrum
cfg_freq            = [];
cfg_freq.method     = 'mtmfft';
cfg_freq.taper      = 'hanning';
cfg_freq.foi        = 1:0.5:80;
cfg_freq.pad        = 'nextpow2';
cfg_freq.keeptrials = 'no';
freq_comp           = ft_freqanalysis(cfg_freq, comp_freq);
clear comp_freq

% Plot
n_plot = min(20, length(comp.label));
figure('Name', 'IC Power Spectra', 'Position', [100 100 1400 900]);
for ic = 1:n_plot
    subplot(4, 5, ic);
    pwr = squeeze(freq_comp.powspctrm(ic, :));
    plot(freq_comp.freq, 10*log10(pwr), 'k', 'LineWidth', 1);
    title(sprintf('IC%d', ic), 'FontSize', 8);
    xlabel('Hz', 'FontSize', 7);
    ylabel('dB', 'FontSize', 7);
    xlim([1 80]); grid on;
    set(gca, 'FontSize', 7);
end
sgtitle('IC Power Spectra (1–80 Hz)');


%% STEP 6b: ICA component rejection
if isempty(ica_reject)
    warning('ica_reject is empty — fill in and rerun this section.');
else
    cfg_rej           = [];
    cfg_rej.component = ica_reject;
    cfg_rej.demean    = 'no';
    data_ica          = ft_rejectcomponent(cfg_rej, comp);
    fprintf('ICA components removed: [%s]\n', num2str(ica_reject));
    log.ica_components_rejected = ica_reject;

    cfg_browser            = [];
    cfg_browser.preproc.demean = 'yes';
    ft_databrowser(cfg_browser, data_ica);
end

save(fullfile(proc_path, [subj_code '_comp_ICA.mat']), 'comp', '-v7.3');

%% STEP 7: ft_interpolatenan — pchip cubic interpolation of TMS artifact window
% Deferred to after ICA so synthetic interpolated samples do not influence
% the decomposition. Uses FieldTrip-native ft_interpolatenan with 50 ms
% pre/post context windows on each side of the NaN artifact gap.
fprintf('\n=== STEP 7: ft_interpolatenan (pchip, 50 ms context) ===\n');

cfg_interp            = [];
cfg_interp.method     = 'pchip';
cfg_interp.prewindow  = 0.05;   % 50 ms of real data before artifact window
cfg_interp.postwindow = 0.05;   % 50 ms of real data after artifact window
data_interp           = ft_interpolatenan(cfg_interp, data_ica);
fprintf('Interpolation complete.\n');

% Visualization: check interpolation at stimulation site
chan_idx = find(strcmp(data_interp.label, plot_chan));
if ~isempty(chan_idx)
    t_vec = data_interp.time{1};   % time axis in seconds
    trial_avg = mean(cat(3, data_interp.trial{:}), 3);
    figure('Name', ['Interpolation check - ' plot_chan]);
    plot(t_vec * 1000, trial_avg(chan_idx, :));
    xlabel('Time (ms)'); ylabel('uV');
    title(['Trial average at ' plot_chan ' after interpolation']);
    xlim([-100 100]);
    xline(0, 'r--', 'TMS');
    xline(tms_win(1), 'g--'); xline(tms_win(2), 'g--');
end

clear data_nan data_ica data_clean

% =========================================================================
% PHASE 2: EEGLAB / TESA
% =========================================================================

%% STEP 8: Convert FieldTrip -> EEGLAB
fprintf('\n=== STEP 8: Converting to EEGLAB ===\n');

EEG         = eeg_emptyset();
EEG.data    = cat(3, data_interp.trial{:});   % [nChan x nSamples x nTrials]
EEG.nbchan  = size(EEG.data, 1);
EEG.pnts    = size(EEG.data, 2);
EEG.trials  = size(EEG.data, 3);
EEG.srate   = data_interp.fsample;
EEG.xmin    = data_interp.time{1}(1);
EEG.xmax    = data_interp.time{1}(end);
EEG.times   = data_interp.time{1} * 1000;
EEG.setname = [subj_code '_postICA'];
for i = 1:length(data_interp.label)
    EEG.chanlocs(i).labels = data_interp.label{i};
end
EEG = eeg_checkset(EEG);
fprintf('Channels: %d | Trials: %d | Rate: %d Hz\n', EEG.nbchan, EEG.trials, EEG.srate);
clear data_interp

%% STEP 9: Channel locations (3D) + condition labels
% standard-10-5-cap385.elp provides 3D (x,y,z) coordinates — required for
% SOUND's forward head model, spherical interpolation (Step 15), and
% accurate EEGLAB topoplot rendering.
% Different from EEG1010.lay used in Steps 6-6b (2D, FieldTrip only).
fprintf('\n=== STEP 9: Channel locations + condition labels ===\n');

EEG          = pop_chanedit(EEG, 'lookup', 'standard-10-5-cap385.elp');
trial_labels = trl_ds(:, 4);   % trigger values from preserved trl matrix

for i = 1:EEG.trials
    lbl = trial_labels(min(i, length(trial_labels)));
    EEG.event(i).latency   = 1 + (i-1) * EEG.pnts;
    EEG.event(i).epoch     = i;
    EEG.event(i).type      = lbl;
    if lbl == 1;      EEG.event(i).condition = 'subthreshold';
    elseif lbl == 3;  EEG.event(i).condition = 'suprathreshold';
    else;             EEG.event(i).condition = 'unknown';
    end
end
EEG = eeg_checkset(EEG, 'eventconsistency');

%% STEP 10: SOUND + visualization
% Attenuates spatially uncorrelated channel noise using a forward head model.
% Channel locations (Step 9) must be loaded before this step.
% Applied before high-pass filter and SSP-SIR.
fprintf('\n=== STEP 10: SOUND (lambda=0.1, iter=10) ===\n');

EEG = tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 10);
pop_eegplot(EEG, 1, 1, 1);

%% STEP 11: High-pass filter 0.5 Hz
fprintf('\n=== STEP 11: High-pass filter (0.5 Hz) ===\n');

EEG = pop_eegfiltnew(EEG, 0.5, []);
pop_eegplot(EEG, 1, 1, 1);

%% STEP 12: Baseline correction
% Full pre-stimulus window used hevere for drift removal.
% For TEP analysis apply the narrower -100 to -10 ms window (T4TE manuscript).
fprintf('\n=== STEP 12: Baseline correction (-1500 to -10 ms) ===\n');

EEG = pop_rmbase(EEG, [-1500 -10]);
pop_eegplot(EEG, 1, 1, 1);


%% STEP 13: SSP-SIR
fprintf('\n=== STEP 13: SSP-SIR ===\n');

EEG_preSSP = EEG;

pc_vals = [3, 4, 5];

EEG_compare = cell(1, 3);
for k = 1:3
    fprintf('Running pc=%d...\n', pc_vals(k));
    EEG_compare{k} = tesa_sspsir(EEG_preSSP, 'artScale', 'automatic', 'PC', pc_vals(k));
end

c3_idx  = find(strcmp({EEG.chanlocs.labels}, 'C3'));
fc5_idx = find(strcmp({EEG.chanlocs.labels}, 'FC5'));
colors  = {'b', 'k', 'r'};
labels  = {'pc=3', 'pc=4', 'pc=5'};

figure; hold on;
for k = 1:3
    plot(EEG.times, mean(EEG_compare{k}.data(c3_idx,:,:), 3), colors{k}, 'LineWidth', 1.5);
end
xline(0,'--k','TMS'); xlim([-200 500]);
legend(labels); title('SSP-SIR — C3'); grid on;

figure; hold on;
for k = 1:3
    plot(EEG.times, mean(EEG_compare{k}.data(fc5_idx,:,:), 3), colors{k}, 'LineWidth', 1.5);
end
xline(0,'--k','TMS'); xlim([-50 150]);
legend(labels); title('SSP-SIR — FC5'); grid on;

EEG = tesa_sspsir(EEG_preSSP, 'artScale', 'automatic', 'PC', sspsir_pc);
fprintf('SSP-SIR applied with pc=%d.\n', sspsir_pc);
log.sspsir_pc = sspsir_pc;
clear EEG_preSSP EEG_compare

%% STEP 14: Notch filter — 50 Hz and harmonics (100, 150, 200 Hz)
% Applied AFTER SSP-SIR to preserve spatial covariance for source modelling.
% Upper limit 200 Hz — 250 Hz is Nyquist for 500 Hz data.
fprintf('\n=== STEP 14: Notch filters (50, 100, 150, 200 Hz) ===\n');

for f = 50:50:200
    EEG = pop_eegfiltnew(EEG, f-2, f+2, [], 1);   % revfilt=1 -> bandstop
    fprintf('  %d Hz\n', f);
end

%% STEP 15: Bad channel interpolation + average reference
% Only channels in chans_interpolate (central) are interpolated.
% Peripheral bad channels were removed at Step 5.
fprintf('\n=== STEP 15: Channel interpolation + average reference ===\n');

if ~isempty(chans_interpolate)
    EEG = pop_interp(EEG, EEG.chaninfo.nodatchans, 'spherical');
    fprintf('Interpolated: %s\n', strjoin(chans_interpolate, ', '));
end
EEG = pop_reref(EEG, []);

%% STEP 16: TEP visualization
fprintf('\n=== STEP 16: TEP visualization ===\n');

TEP      = mean(EEG.data, 3);
t        = EEG.times;
chanlabs = {EEG.chanlocs.labels};

c_idx = find(strcmp(chanlabs, plot_chan));
if ~isempty(c_idx)
    figure('Name', ['TEP - ' plot_chan]);
    plot(t, TEP(c_idx, :), 'k', 'LineWidth', 1.2);
    xlabel('Time (ms)'); ylabel('\muV');
    title(['TEP — ' plot_chan ' (n=' num2str(EEG.trials) ')']);
    xlim([-200 500]); xline(0, 'r--', 'TMS'); grid on;
end

% Hjorth-Laplacian C3: center minus mean of 4 neighbours (T4TE manuscript)
c3_idx = find(strcmp(chanlabs, 'C3'));
nb_idx = find(ismember(chanlabs, {'FC1','CP1','FC5','CP5'}));

if ~isempty(c3_idx) && length(nb_idx) == 4
    Hjorth_avg = mean(squeeze(EEG.data(c3_idx,:,:)) - ...
                 0.25 * squeeze(sum(EEG.data(nb_idx,:,:), 1)), 2);
    figure('Name', 'TEP - Hjorth C3');
    plot(t, Hjorth_avg, 'k', 'LineWidth', 1.2);
    xlabel('Time (ms)'); ylabel('\muV');
    title('TEP — Hjorth-Laplacian C3 (FC1, CP1, FC5, CP5)');
    xlim([-200 500]); xline(0, 'r--', 'TMS'); grid on;
else
    fprintf('Warning: Hjorth neighbours not all found.\n');
end

%TEP per condition — sub vs supra apart
sub_idx   = find([EEG.event.type] == 1);
supra_idx = find([EEG.event.type] == 3);

c3_idx = find(strcmp({EEG.chanlocs.labels}, 'C3'));

figure('Name', 'TEP per conditie');
plot(EEG.times, mean(EEG.data(c3_idx,:,sub_idx), 3), 'b', 'LineWidth', 1.5); hold on;
plot(EEG.times, mean(EEG.data(c3_idx,:,supra_idx), 3), 'r', 'LineWidth', 1.5);
xline(0, '--k', 'TMS');
xlim([-200 500]); xlabel('ms'); ylabel('µV');
legend('Subthreshold', 'Suprathreshold');
title([subj_code ' — TEP C3 per condition']); grid on;

% Hjorth per conditie
matlabnb_idx = find(ismember({EEG.chanlocs.labels}, {'FC1','CP1','FC5','CP5'}));

Hjorth_sub   = mean(squeeze(EEG.data(c3_idx,:,sub_idx))   - 0.25*squeeze(sum(EEG.data(nb_idx,:,sub_idx),1)),   2);
Hjorth_supra = mean(squeeze(EEG.data(c3_idx,:,supra_idx)) - 0.25*squeeze(sum(EEG.data(nb_idx,:,supra_idx),1)), 2);

figure('Name', 'Hjorth per conditie');
plot(EEG.times, Hjorth_sub,   'b', 'LineWidth', 1.5); hold on;
plot(EEG.times, Hjorth_supra, 'r', 'LineWidth', 1.5);
xline(0, '--k', 'TMS');
xlim([-200 500]); xlabel('ms'); ylabel('µV');
legend('Subthreshold', 'Suprathreshold');
title([subj_code ' — Hjorth C3 per condition']); grid on;

% Butterfly plot — alle kanalen
figure('Name', 'Butterfly plot');
plot(EEG.times, mean(EEG.data, 3), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold on;
plot(EEG.times, mean(EEG.data(c3_idx,:,:), 3), 'k', 'LineWidth', 2);
xline(0, '--r', 'TMS');
xlim([-200 500]); xlabel('ms'); ylabel('µV');
title([subj_code ' — Butterfly (grey = all channels, black=C3)']); grid on;

% Save figures
figs = get(0, 'Children');
for i = 1:length(figs)
    exportgraphics(figs(i), fullfile(proc_path, ...
        sprintf('%s_fig%d.png', subj_code, figs(i).Number)), 'Resolution', 150);
end
fprintf('figures saved.\n');


%% STEP 17: Save dataset + preprocessing log
fprintf('\n=== STEP 17: Saving ===\n');

EEG.setname = [subj_code '_eeg_final'];
EEG = pop_saveset(EEG, 'filename', [subj_code '_eeg_final.set'], 'filepath', proc_path);

n_sub_final   = sum([EEG.event.type] == 1);
n_supra_final = sum([EEG.event.type] == 3);
log.n_sub_final        = n_sub_final;
log.n_supra_final      = n_supra_final;
log.n_trials_final     = EEG.trials;
log.pct_rejected_total = (1 - EEG.trials / n_orig) * 100;

save(fullfile(proc_path, [subj_code '_prepro_log.mat']), 'log');

fid = fopen(fullfile(proc_path, [subj_code '_prepro_log.txt']), 'w');
fprintf(fid, 'Preprocessing Log — T4TE v9\n');
fprintf(fid, 'Subject: %s | Date: %s | File: %s\n', log.subj_code, log.date, log.raw_file);
fprintf(fid, '-----------------------------------------\n');
fprintf(fid, 'INPUT  Sub: %d | Supra: %d | Total: %d\n', ...
    log.n_trials_input_sub, log.n_trials_input_supra, log.n_trials_input_total);
fprintf(fid, 'REJECTION (visual)\n');
fprintf(fid, '  Total:  %d rejected | %d remaining\n', ...
    length(log.trials_rejected_visual), log.n_trials_after_visual);
fprintf(fid, '  Sub:    %d rejected\n', log.n_trials_rej_visual_sub);
fprintf(fid, '  Supra:  %d rejected\n', log.n_trials_rej_visual_supra);

% v9: safe strjoin — prints 'none' instead of blank when list is empty
if isempty(log.bad_chans_removed)
    fprintf(fid, 'CHANNELS removed: none\n');
else
    fprintf(fid, 'CHANNELS removed: %s\n', strjoin(log.bad_chans_removed, ', '));
end
if isempty(log.bad_chans_for_interp)
    fprintf(fid, 'CHANNELS interp:  none\n');
else
    fprintf(fid, 'CHANNELS interp:  %s\n', strjoin(log.bad_chans_for_interp, ', '));
end
fprintf(fid, 'ICA rejected:     [%s]\n', num2str(log.ica_components_rejected));
fprintf(fid, 'SSP-SIR PC:       %d\n',   log.sspsir_pc);
fprintf(fid, 'FINAL  Sub: %d | Supra: %d | Total: %d (%.1f%% rejected)\n', ...
    log.n_sub_final, log.n_supra_final, log.n_trials_final, log.pct_rejected_total);
fprintf(fid, '=========================================\n');
fclose(fid);

fprintf('Saved: %s_eeg_final.set | %s_prepro_log.txt\n', subj_code, subj_code);
