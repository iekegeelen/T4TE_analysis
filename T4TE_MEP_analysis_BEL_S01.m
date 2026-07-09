%% T4TE_MEP_analysis_BEL_S01.m
% MEP computation and interactive trial-by-trial validation
% Adapted from Elena Mongiardini's MEP_validation.m for the T4TE study
%
% v2 changes vs v1:
%   - EEG rejection mask: trials removed during visual inspection (Step 5)
%     in the preprocessing pipeline are automatically set to NaN before
%     interactive validation, so MEP and TEP analyses are always based on
%     the same trial set.
%   - BEL_S02 flag: option to additionally exclude aborted suprathreshold
%     trials (cable failure) from the EMG trial numbering.
%   - Summary statistics split by rejection source (EEG vs MEP-quality).
%
% INPUT:  <subj_code>_emg.mat         — FieldTrip EMG struct (all original trials,
%                                        saved at Step 3 of preprocessing pipeline)
%         <subj_code>_prepro_log.mat  — preprocessing log (saved at Step 5/17),
%                                        must contain log.trials_rejected_visual
%
% OUTPUT: MEP       — struct with fields:
%                     .subthreshold   : peak-to-peak amplitude per trial (µV)
%                     .suprathreshold : peak-to-peak amplitude per trial (µV)
%                     NaN = rejected trial (EEG mask OR MEP quality)
%         saved to: <proc_path>/<subj_code>_MEP.mat
%
% Trial indexing:
%   All vectors are indexed over the ORIGINAL trial numbering from data_emg
%   (i.e., before any rejection). This preserves direct correspondence with
%   trialinfo and with the preprocessing log indices.
%
% NOTE: MEP detection is only physiologically meaningful in the
%       SUPRATHRESHOLD condition. Subthreshold trials will not produce
%       a visible MEP, but amplitudes are computed for completeness.
%       During validation, reject subthreshold trials only if the signal
%       is clearly artefactual (e.g., large pre-stimulus burst).
%       EEG-based trial exclusion is most critical for suprathreshold
%       trials used in TEP analyses.
%
% Ieke | T4TE study | CIMeC, Trento | 2025
% =========================================================================

clc; clear all; close all;

%% -------------------------------------------------------------------------
%  PARAMETERS — edit these for each subject
% -------------------------------------------------------------------------
subj_code  = 'BEL_S01';
proc_path  = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/BEL_S01/processed/';

emg_channel = 'FDI';   % primary MEP channel; alternatives: 'APB', 'EDC'

% MEP search window (ms post-TMS pulse, i.e., post t=0)
mep_win_ms = [20, 60];

% Baseline window for DC correction (ms pre-TMS, must be negative)
baseline_win_ms = -10;   % use everything before -10 ms

% Condition codes in trialinfo
cond_codes  = [1, 3];
cond_labels = {'subthreshold', 'suprathreshold'};

% -------------------------------------------------------------------------
% BEL_S02 ONLY: aborted suprathreshold trials due to cable failure.
% These were removed AFTER data_emg was saved, so they must be excluded
% manually here. Set to true for BEL_S02, leave false for all others.
% n_aborted: number of aborted trials at the START of the suprathreshold
% block (they appear first in the trial sequence for that condition).
% -------------------------------------------------------------------------
remove_aborted_supra = false;   % BEL_S01: no aborted trials
n_aborted            = 0;

% =========================================================================

%% -------------------------------------------------------------------------
%  Load EMG data
% -------------------------------------------------------------------------
fprintf('Loading EMG data for %s...\n', subj_code);
load(fullfile(proc_path, [subj_code '_emg.mat']));  % loads data_emg

fs        = data_emg.fsample;   % 5000 Hz
ch_labels = data_emg.label;     % {'APB','FDI','EDC'}
n_trials  = length(data_emg.trial);

% Find the MEP channel index
ch_idx = find(strcmp(ch_labels, emg_channel));
if isempty(ch_idx)
    error('Channel %s not found. Available: %s', emg_channel, strjoin(ch_labels, ', '));
end

% Time vector (same for all trials)
time_vec = data_emg.time{1};

% Baseline indices (pre-stimulus, for DC correction)
baseline_idx = time_vec < (baseline_win_ms / 1000);

% Convert MEP window from ms to sample indices
[~, mep_idx(1)] = min(abs(time_vec - mep_win_ms(1)/1000));
[~, mep_idx(2)] = min(abs(time_vec - mep_win_ms(2)/1000));

fprintf('MEP window:  %.0f ms to %.0f ms (samples %d to %d)\n', ...
    mep_win_ms(1), mep_win_ms(2), mep_idx(1), mep_idx(2));
fprintf('Baseline:    before %.0f ms (%d samples)\n', ...
    baseline_win_ms, sum(baseline_idx));
fprintf('Total trials in EMG struct: %d\n', n_trials);

%% -------------------------------------------------------------------------
%  Automatic MEP computation (peak-to-peak on DC-corrected signal)
% -------------------------------------------------------------------------
mep_auto = nan(1, n_trials);
mep_max  = nan(1, n_trials);
mep_min  = nan(1, n_trials);

for tr = 1:n_trials
    sig_full = data_emg.trial{tr}(ch_idx, :);
    sig_full = sig_full - mean(sig_full(baseline_idx));
    sig      = sig_full(mep_idx(1):mep_idx(2));
    mep_max(tr)  = max(sig);
    mep_min(tr)  = min(sig);
    mep_auto(tr) = mep_max(tr) - mep_min(tr);
end

fprintf('\nAutomatic MEP computation complete.\n');
fprintf('Suprathreshold median MEP (all trials, pre-mask): %.1f uV\n', ...
    median(mep_auto(data_emg.trialinfo == 3), 'omitnan'));

%% -------------------------------------------------------------------------
%  EEG rejection mask
%  Trials removed during visual inspection in the preprocessing pipeline
%  (Step 5, log.trials_rejected_visual) are set to NaN here so that MEP
%  and TEP analyses are always computed on the same trial set.
% -------------------------------------------------------------------------
fprintf('\n=== EEG rejection mask ===\n');

mep_validated = mep_auto;   % start from automatic values
eeg_mask      = true(1, n_trials);   % true = kept, false = excluded by EEG

log_file = fullfile(proc_path, [subj_code '_prepro_log.mat']);

if exist(log_file, 'file')
    load(log_file, 'log');

    if isfield(log, 'trials_rejected_visual') && ~isempty(log.trials_rejected_visual)
        eeg_rejected = log.trials_rejected_visual;

        % Validate indices are within the EMG trial count
        out_of_range = eeg_rejected(eeg_rejected < 1 | eeg_rejected > n_trials);
        if ~isempty(out_of_range)
            warning('%d rejected trial indices are out of range [1 %d] and will be skipped: %s', ...
                numel(out_of_range), n_trials, num2str(out_of_range));
            eeg_rejected = eeg_rejected(eeg_rejected >= 1 & eeg_rejected <= n_trials);
        end

        eeg_mask(eeg_rejected) = false;
        mep_validated(eeg_rejected) = NaN;

        fprintf('Log loaded: %s\n', log_file);
        fprintf('EEG mask applied: %d trials excluded (%.1f%% of %d)\n', ...
            numel(eeg_rejected), 100 * numel(eeg_rejected) / n_trials, n_trials);

        % Per-condition breakdown of EEG-excluded trials
        for c = 1:numel(cond_codes)
            cond_all = find(data_emg.trialinfo == cond_codes(c));
            n_eeg_rej_cond = sum(~eeg_mask(cond_all));
            fprintf('  %s: %d EEG-rejected\n', cond_labels{c}, n_eeg_rej_cond);
        end
    else
        fprintf('log.trials_rejected_visual is empty — no EEG trials were rejected.\n');
    end
else
    warning('Preprocessing log not found: %s\nNo EEG mask applied. MEP and TEP trial sets may not match.', log_file);
end

%% -------------------------------------------------------------------------
%  BEL_S02: additional exclusion of aborted suprathreshold trials
% -------------------------------------------------------------------------
if remove_aborted_supra && n_aborted > 0
    fprintf('\n=== BEL_S02: excluding %d aborted suprathreshold trials ===\n', n_aborted);

    supra_all = find(data_emg.trialinfo == 3);
    if n_aborted > numel(supra_all)
        error('n_aborted (%d) exceeds number of suprathreshold trials (%d).', ...
            n_aborted, numel(supra_all));
    end
    aborted_idx = supra_all(1:n_aborted);
    eeg_mask(aborted_idx)      = false;
    mep_validated(aborted_idx) = NaN;
    fprintf('Aborted trial indices: %s\n', num2str(aborted_idx));
end

%% -------------------------------------------------------------------------
%  Interactive trial-by-trial validation
%  Only trials NOT already set to NaN by the EEG mask are shown.
% -------------------------------------------------------------------------
fprintf('\n=== Interactive MEP validation ===\n');
fprintf('Keys: SPACE = accept | R = reject (set NaN) | C = manually reselect max/min\n');
fprintf('EEG-masked trials are skipped automatically.\n\n');

for c = 1:numel(cond_codes)

    cond_trials = find(data_emg.trialinfo == cond_codes(c));
    valid_trials = cond_trials(eeg_mask(cond_trials));   % skip EEG-excluded

    fprintf('\n--- Condition: %s (%d trials after EEG mask, %d total) ---\n', ...
        cond_labels{c}, numel(valid_trials), numel(cond_trials));

    for tt = 1:numel(valid_trials)
        tr = valid_trials(tt);

        % DC-corrected full signal for plotting
        sig = data_emg.trial{tr}(ch_idx, :);
        sig = sig - mean(sig(baseline_idx));

        fig = figure('Name', sprintf('%s | %s | Trial %d/%d (orig #%d)', ...
            subj_code, cond_labels{c}, tt, numel(valid_trials), tr), ...
            'NumberTitle', 'off');

        plot(time_vec * 1000, sig, 'b', 'LineWidth', 0.8);
        hold on;

        yl = [-500 500];

        % Shade MEP search window
        patch([mep_win_ms(1) mep_win_ms(2) mep_win_ms(2) mep_win_ms(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha', 0.4);
        uistack(findobj(gca,'Type','patch'), 'bottom');

        % Mark automatic max/min within MEP window
        [~, maxLoc_rel] = max(sig(mep_idx(1):mep_idx(2)));
        [~, minLoc_rel] = min(sig(mep_idx(1):mep_idx(2)));
        maxIdx = mep_idx(1) + maxLoc_rel - 1;
        minIdx = mep_idx(1) + minLoc_rel - 1;

        hMax = plot(time_vec(maxIdx)*1000, sig(maxIdx), 'ro', ...
            'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','Max');
        hMin = plot(time_vec(minIdx)*1000, sig(minIdx), 'go', ...
            'MarkerSize', 10, 'LineWidth', 2, 'DisplayName','Min');

        ylim(yl);
        xlim([-100 200]);
        xlabel('Time (ms)');
        ylabel([emg_channel ' amplitude (\muV)']);
        title(sprintf('%s | %s | Trial %d (orig #%d) | MEP = %.1f \muV', ...
            subj_code, cond_labels{c}, tt, tr, mep_validated(tr)));
        legend([hMax hMin], {'Max','Min'}, 'Location','best');
        grid on;
        xline(0, '--k', 'TMS', 'LineWidth', 1);

        % Wait for keypress
        keyPressed = false;
        while ~keyPressed
            waitforbuttonpress;
            k = get(fig, 'CurrentCharacter');

            if k == 'r' || k == 'R'
                mep_validated(tr) = NaN;
                fprintf('Trial %d (orig #%d, %s) rejected → NaN\n', tt, tr, cond_labels{c});
                keyPressed = true;

            elseif k == 'c' || k == 'C'
                fprintf('Click: (1) new max position, then (2) new min position\n');
                [x_click, ~] = ginput(2);
                newMaxIdx = round(x_click(1)/1000 * fs) + find(time_vec >= 0, 1);
                newMinIdx = round(x_click(2)/1000 * fs) + find(time_vec >= 0, 1);
                newMaxIdx = max(1, min(newMaxIdx, length(time_vec)));
                newMinIdx = max(1, min(newMinIdx, length(time_vec)));
                set(hMax, 'XData', time_vec(newMaxIdx)*1000, 'YData', sig(newMaxIdx));
                set(hMin, 'XData', time_vec(newMinIdx)*1000, 'YData', sig(newMinIdx));
                mep_validated(tr) = sig(newMaxIdx) - sig(newMinIdx);
                title(sprintf('%s | %s | Trial %d (orig #%d) | MEP = %.1f \muV [updated]', ...
                    subj_code, cond_labels{c}, tt, tr, mep_validated(tr)));
                fprintf('Trial %d updated → MEP = %.1f uV\n', tr, mep_validated(tr));

            elseif k == ' '
                fprintf('Trial %d (orig #%d, %s) accepted → MEP = %.1f uV\n', ...
                    tt, tr, cond_labels{c}, mep_validated(tr));
                keyPressed = true;
            end
        end

        close(fig);
    end
end

%% -------------------------------------------------------------------------
%  Package MEP struct by condition
% -------------------------------------------------------------------------
MEP = struct();
for c = 1:numel(cond_codes)
    cond_trials = find(data_emg.trialinfo == cond_codes(c));
    MEP.(cond_labels{c}) = mep_validated(cond_trials);
end

fprintf('\n=== Validation complete ===\n');
fprintf('\nSummary per condition:\n');
for c = 1:numel(cond_codes)
    cond_trials  = find(data_emg.trialinfo == cond_codes(c));
    n_total      = numel(cond_trials);
    n_eeg_rej    = sum(~eeg_mask(cond_trials));
    n_mep_rej    = sum(isnan(MEP.(cond_labels{c}))) - n_eeg_rej;
    n_valid      = sum(~isnan(MEP.(cond_labels{c})));

    fprintf('  %s:\n', upper(cond_labels{c}));
    fprintf('    Total trials:          %d\n', n_total);
    fprintf('    EEG-masked (pipeline): %d\n', n_eeg_rej);
    fprintf('    MEP-quality rejected:  %d\n', n_mep_rej);
    fprintf('    Valid (used in both):  %d\n', n_valid);
    if n_valid > 0
        fprintf('    Median MEP amplitude:  %.1f uV\n', ...
            median(MEP.(cond_labels{c}), 'omitnan'));
    end
end

%% -------------------------------------------------------------------------
%  Save
% -------------------------------------------------------------------------
save_path = fullfile(proc_path, [subj_code '_MEP.mat']);
save(save_path, 'MEP', 'eeg_mask', 'mep_auto', 'mep_validated');
fprintf('\nSaved: %s\n', save_path);
fprintf('Variables saved: MEP, eeg_mask, mep_auto, mep_validated\n');

%% -------------------------------------------------------------------------
%  Summary figure
% -------------------------------------------------------------------------
figure('Name', [subj_code ' — MEP summary'], 'Position', [100 100 800 500]);

for c = 1:numel(cond_codes)
    subplot(2, 2, c);
    vals = MEP.(cond_labels{c});
    histogram(vals(~isnan(vals)), 20, ...
        'FaceColor', [0.4 0.6 0.8] + (c-1)*[0.4 -0.2 -0.4]);
    xlabel('Peak-to-peak amplitude (\muV)');
    ylabel('Trial count');
    n_valid = sum(~isnan(vals));
    title(sprintf('%s (n=%d valid)', cond_labels{c}, n_valid));
    grid on;
end

% Suprathreshold trial-by-trial trace (valid only)
subplot(2, 2, [3 4]);
supra_vals = MEP.suprathreshold;
supra_valid_idx = find(~isnan(supra_vals));
plot(supra_valid_idx, supra_vals(supra_valid_idx), 'o-', ...
    'Color', [0.8 0.4 0.4], 'MarkerSize', 4, 'LineWidth', 0.8);
xlabel('Trial index (within condition)');
ylabel('Peak-to-peak MEP (\muV)');
title('Suprathreshold: trial-by-trial MEP (valid trials only)');
grid on;
yline(median(supra_vals, 'omitnan'), '--k', 'Median', 'LineWidth', 1.2);

sgtitle(sprintf('%s — %s MEP amplitudes (DC-corrected, EEG-synced)', subj_code, emg_channel));

exportgraphics(gcf, fullfile(proc_path, [subj_code '_MEP_summary.png']), 'Resolution', 150);
fprintf('Figure saved: %s_MEP_summary.png\n', subj_code);
