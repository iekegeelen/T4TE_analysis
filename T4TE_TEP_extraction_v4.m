% =========================================================================
% T4TE Study 1.1 — TEP Extraction  (v4)
% =========================================================================
%
% Verified against BEL_S01_eeg_final.set (inspected 2026-06-19):
%   - Linked EEGLAB set (.set + .fdt), trial-level data [chan x time x trials]
%   - srate = 500 Hz, time axis = -1500 to 498 ms (1000 samples)
%   - EEG.epoch(i).eventtype = numeric 1 (sub, 90% rMT) or 3 (supra, 110% rMT)
%
% Channel handling:
%   Bad peripheral channels were hard-removed per subject in pipeline v9
%   (Step 5). Central channels were interpolated within the pipeline (Step 15).
%   Channel counts therefore vary across subjects (56-60 ch).
%   Solution: compute channel intersection across all subjects on first pass,
%   then extract TEPs only on those 50 common channels. No new interpolation.
%   All key channels confirmed present in all subjects:
%   C3, FC1, CP1, FC5, CP5.
%
% Steps:
%   1. First pass: load each subject, collect channel labels, find intersection
%   2. Second pass: load each subject, select common channels, baseline correct,
%      split by condition, accumulate subject-average TEPs
%   3. Grand average + SEM
%   4. Hjorth-Laplacian TEP at C3
%   5. TEP component extraction via TEP_peak.m (interactive)
%   6. Topoplots via custum_topoplot.m
%   7. Save TEP_1results.mat
%
% Dependencies: shadedErrorBar.m, TEP_peak.m, custum_topoplot.m,
%               dcc_customized_acticap64.mat, EEGLAB, FieldTrip
%
% Author: Ieke Dresens — T4TE internship, CIMeC Trento, 2025/2026
% =========================================================================




% =========================================================================
% PATHS — edit once
% =========================================================================

proc_path   = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/';
eeglab_path = '/Users/e.w.m.dresens/Documents/MATLAB/eeglab2026.0.0/';
ft_path     = '/Users/e.w.m.dresens/Documents/MATLAB/fieldtrip-20250106/';
tools_path  = '/Users/e.w.m.dresens/Documents/MATLAB/';

addpath(eeglab_path);
addpath(ft_path);
addpath(tools_path);
ft_defaults;
eeglab('nogui');

group_path = fullfile(proc_path, 'group');
if ~exist(group_path, 'dir'); mkdir(group_path); end

% =========================================================================
% SETTINGS
% =========================================================================

subj_list = {'BEL_S01','BEL_S02','BEL_S03','BEL_S04','BEL_S05','BEL_S06', ...
             'BEL_S07','BEL_S08','BEL_S09','BEL_S10','BEL_S11','BEL_S12'};
n_sub = numel(subj_list);

trig_vals   = [1, 3];
cond_labels = {'subthreshold', 'suprathreshold'};
n_cond      = numel(cond_labels);

baseline_win = [-100, -10];   % ms — narrow re-baseline for TEP
artifact_win = [-5,   10];    % ms — blanked in plots only
plot_chan    = 'C3';
hjorth_nb    = {'FC1','CP1','FC5','CP5'};

p_peaks     = [25, 60];
n_peaks     = [45, 100];
peak_window = [10, 200];

yl_fixed = [-5, 5];
col_cond = {[0 0 0], [0 0.45 0.74]};   % black=sub, blue=supra

% =========================================================================
% FIRST PASS — collect channel labels per subject, compute intersection
% =========================================================================

fprintf('\n=== Pass 1: Computing common channel set ===\n');

chan_sets = cell(1, n_sub);   % channel labels per subject

for s = 1:n_sub

    subj      = subj_list{s};
    setfile   = [subj '_eeg_final.set'];
    subdir    = fullfile(proc_path, subj, 'processed');
    full_path = fullfile(subdir, setfile);

    if ~exist(full_path, 'file')
        warning('  %s not found — will be skipped.', setfile);
        chan_sets{s} = {};
        continue
    end

    % Load header only (no data) to read chanlocs
    EEG_hdr = pop_loadset('filename', setfile, 'filepath', subdir, 'loadmode', 'info');

    chan_sets{s} = {EEG_hdr.chanlocs.labels};
    fprintf('  %s: %d channels\n', subj, numel(chan_sets{s}));

end

% Intersection: channels present in every subject that was found
common_chans = chan_sets{find(~cellfun(@isempty, chan_sets), 1)};   % start from first valid
for s = 1:n_sub
    if ~isempty(chan_sets{s})
        common_chans = intersect(common_chans, chan_sets{s}, 'stable');
    end
end

fprintf('\nCommon channel set: %d channels\n', numel(common_chans));

% Confirm key channels present
key_chans = [plot_chan, hjorth_nb];
for k = 1:numel(key_chans)
    if ~ismember(key_chans{k}, common_chans)
        error('Key channel %s missing from common set — cannot proceed.', key_chans{k});
    end
end
fprintf('C3 and all Hjorth neighbours confirmed in common set.\n');

% =========================================================================
% SECOND PASS — load, select common channels, accumulate TEPs
% =========================================================================

fprintf('\n=== Pass 2: Extracting TEPs ===\n\n');

TEP_data   = cell(1, n_cond);   % {c}: [n_common x n_time x n_sub]
trial_data = cell(1, n_cond);   % {c}.(subj): [n_common x n_time x n_trials]
for c = 1:n_cond
    TEP_data{c}   = [];
    trial_data{c} = struct();
end

t        = [];
chanlabs = {};   % will be set to common_chans on first subject
n_time   = 0;

for s = 1:n_sub

    subj      = subj_list{s};
    subdir    = fullfile(proc_path, subj, 'processed');
    setfile   = [subj '_eeg_final.set'];
    full_path = fullfile(subdir, setfile);

    fprintf('--- Subject %d/%d: %s ---\n', s, n_sub, subj);

    if ~exist(full_path, 'file')
        warning('  Skipping — file not found.');
        continue
    end

    EEG = pop_loadset('filename', setfile, 'filepath', subdir);
    EEG = eeg_checkset(EEG);

    fprintf('  Loaded: %d ch × %d samples × %d epochs\n', ...
        EEG.nbchan, EEG.pnts, EEG.trials);

    % Store time axis on first subject
    if isempty(t)
        t      = EEG.times;   % ms [1 x n_time]
        n_time = EEG.pnts;
        chanlabs = common_chans;
    end

    % -----------------------------------------------------------------
    % Select common channels
    % -----------------------------------------------------------------
    [~, chan_idx] = ismember(common_chans, {EEG.chanlocs.labels});

    if any(chan_idx == 0)
        missing = common_chans(chan_idx == 0);
        error('  %s: common channel(s) not found: %s', subj, strjoin(missing, ', '));
    end

    EEG.data     = EEG.data(chan_idx, :, :);
    EEG.chanlocs = EEG.chanlocs(chan_idx);
    EEG.nbchan   = numel(common_chans);

    % -----------------------------------------------------------------
    % Narrow baseline correction: subtract mean of [-100 -10] ms
    % -----------------------------------------------------------------
    [~, bl_s] = min(abs(t - baseline_win(1)));
    [~, bl_e] = min(abs(t - baseline_win(2)));

    bl_mean  = mean(EEG.data(:, bl_s:bl_e, :), 2);   % [ch x 1 x trial]
    EEG.data = EEG.data - bl_mean;

    % -----------------------------------------------------------------
    % Get condition trigger per epoch via EEG.epoch(i).eventtype
    % -----------------------------------------------------------------
    epoch_trig = zeros(1, EEG.trials);
    for ep = 1:EEG.trials
        raw = EEG.epoch(ep).eventtype;
        if iscell(raw); raw = raw{1}; end
        epoch_trig(ep) = double(raw);
    end

    fprintf('  Triggers — type 1: %d   type 3: %d\n', ...
        sum(epoch_trig == 1), sum(epoch_trig == 3));

    % -----------------------------------------------------------------
    % Split by condition and accumulate
    % -----------------------------------------------------------------
    for c = 1:n_cond

        idx_trials = find(epoch_trig == trig_vals(c));

        if isempty(idx_trials)
            warning('  %s: no trials for %s (trigger %d).', ...
                subj, cond_labels{c}, trig_vals(c));
            continue
        end

        cond_data = EEG.data(:, :, idx_trials);   % [n_common x n_time x n_trials]
        fprintf('  %s: %d trials\n', cond_labels{c}, size(cond_data, 3));

        sub_tep = mean(cond_data, 3);              % [n_common x n_time]

        TEP_data{c}          = cat(3, TEP_data{c}, sub_tep);
        trial_data{c}.(subj) = cond_data;

    end

    clear EEG bl_mean cond_data sub_tep

end

fprintf('\nAll subjects processed.\n');
fprintf('Grand average matrix: %d channels × %d samples × %d subjects\n', ...
    size(TEP_data{1}, 1), size(TEP_data{1}, 2), size(TEP_data{1}, 3));

% =========================================================================
% CHANNEL INDICES (on common channel set)
% =========================================================================

c3_idx = find(strcmp(chanlabs, plot_chan));
nb_idx = find(ismember(chanlabs, hjorth_nb));

idx_pre  = t <= artifact_win(1);
idx_post = t >= artifact_win(2);

% =========================================================================
% FIGURE 1 — GRAND AVERAGE TEP AT C3, SEPARATE SUBPLOTS
% =========================================================================

fprintf('\n--- Figure 1: Grand average TEP ---\n');

fig_ga = figure('Name', 'Grand Average TEP', 'Position', [50 50 1000 700]);

for c = 1:n_cond

    n_s    = size(TEP_data{c}, 3);
    ga     = mean(TEP_data{c}, 3);
    ga_sem = std(TEP_data{c}, 0, 3) / sqrt(n_s);
    m      = ga(c3_idx, :);
    e      = ga_sem(c3_idx, :);
    lp     = {'-', 'Color', col_cond{c}};

    subplot(n_cond, 1, c); hold on;

    patch([artifact_win(1) artifact_win(2) artifact_win(2) artifact_win(1)], ...
          [yl_fixed(1) yl_fixed(1) yl_fixed(2) yl_fixed(2)], ...
          [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

    shadedErrorBar(t(idx_pre),  m(idx_pre),  e(idx_pre),  'lineprops', lp, 'transparent', true);
    shadedErrorBar(t(idx_post), m(idx_post), e(idx_post), 'lineprops', lp, 'transparent', true);

    xline(0, 'r--', 'LineWidth', 1.2);
    plot(xlim, [0 0], 'k', 'LineWidth', 0.7);
    for lat = [25 45 60 100]
        xline(lat, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.9);
    end

    xlim([-100 300]); ylim(yl_fixed);
    xlabel('Time (ms)'); ylabel('Amplitude (µV)');
    title(sprintf('Grand Average TEP — %s — %s  (n = %d)', plot_chan, cond_labels{c}, n_s));
    grid on; box off; set(gca, 'FontSize', 12);

end

saveas(fig_ga, fullfile(group_path, 'TEP_grandaverage_C3.png'));

% =========================================================================
% FIGURE 2 — BOTH CONDITIONS OVERLAID
% =========================================================================

fprintf('--- Figure 2: Overlay ---\n');

fig_ov = figure('Name', 'Grand Average TEP — Overlay', 'Position', [80 80 900 480]);
hold on;

patch([artifact_win(1) artifact_win(2) artifact_win(2) artifact_win(1)], ...
      [yl_fixed(1) yl_fixed(1) yl_fixed(2) yl_fixed(2)], ...
      [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

leg_h = gobjects(1, n_cond);

for c = 1:n_cond
    n_s = size(TEP_data{c}, 3);
    m   = squeeze(mean(TEP_data{c}(c3_idx, :, :), 3));
    e   = squeeze(std(TEP_data{c}(c3_idx, :, :), 0, 3)) / sqrt(n_s);
    lp  = {'-', 'Color', col_cond{c}};

    shadedErrorBar(t(idx_pre),  m(idx_pre),  e(idx_pre),  'lineprops', lp, 'transparent', true);
    H = shadedErrorBar(t(idx_post), m(idx_post), e(idx_post), 'lineprops', lp, 'transparent', true);
    leg_h(c) = H.mainLine;
end

xline(0, 'r--', 'LineWidth', 1.2);
plot(xlim, [0 0], 'k', 'LineWidth', 0.7);
for lat = [25 45 60 100]
    xline(lat, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.9);
end

xlim([-100 300]); ylim(yl_fixed);
xlabel('Time (ms)'); ylabel('Amplitude (µV)');
title(sprintf('Grand Average TEP — %s — both conditions', plot_chan));
legend(leg_h, cond_labels, 'Location', 'best', 'Box', 'off');
grid on; box off; set(gca, 'FontSize', 12);

saveas(fig_ov, fullfile(group_path, 'TEP_grandaverage_overlay.png'));

% =========================================================================
% FIGURE 3 — HJORTH-LAPLACIAN TEP AT C3
% =========================================================================

fprintf('--- Figure 3: Hjorth-Laplacian TEP ---\n');

fig_hj = figure('Name', 'Hjorth-Laplacian TEP', 'Position', [100 100 1000 700]);

for c = 1:n_cond

    n_s         = size(TEP_data{c}, 3);
    hjorth_subs = zeros(n_s, n_time);

    for s = 1:n_s
        c3_sig  = squeeze(TEP_data{c}(c3_idx, :, s));
        nb_mean = squeeze(mean(TEP_data{c}(nb_idx, :, s), 1));
        hjorth_subs(s, :) = c3_sig - nb_mean;
    end

    m_hj = mean(hjorth_subs, 1);
    e_hj = std(hjorth_subs, 0, 1) / sqrt(n_s);
    lp   = {'-', 'Color', col_cond{c}};

    subplot(n_cond, 1, c); hold on;

    patch([artifact_win(1) artifact_win(2) artifact_win(2) artifact_win(1)], ...
          [yl_fixed(1) yl_fixed(1) yl_fixed(2) yl_fixed(2)], ...
          [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

    shadedErrorBar(t(idx_pre),  m_hj(idx_pre),  e_hj(idx_pre),  'lineprops', lp, 'transparent', true);
    shadedErrorBar(t(idx_post), m_hj(idx_post), e_hj(idx_post), 'lineprops', lp, 'transparent', true);

    xline(0, 'r--', 'LineWidth', 1.2);
    plot(xlim, [0 0], 'k', 'LineWidth', 0.7);
    for lat = [25 45 60 100]
        xline(lat, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.9);
    end

    xlim([-100 300]); ylim(yl_fixed);
    xlabel('Time (ms)'); ylabel('Amplitude (µV)');
    title(sprintf('Hjorth-Laplacian TEP — %s — %s  (n = %d)', plot_chan, cond_labels{c}, n_s));
    grid on; box off; set(gca, 'FontSize', 12);

end

saveas(fig_hj, fullfile(group_path, 'TEP_Hjorth_C3.png'));

% =========================================================================
% TEP COMPONENT EXTRACTION — TEP_peak.m (interactive, once per condition)
% =========================================================================
%
% Opens a figure of the C3 grand average and prompts in the command window:
%   "Modify peaks? (1=yes / 0=no)"
%   0 = accept automatic findpeaks result
%   1 = click peak positions and intervals manually

fprintf('\n--- TEP component extraction ---\n');
fprintf('TEP_peak.m will prompt "Modify peaks?" per condition.\n\n');

PEAK = struct();

for c = 1:n_cond

    fprintf('Condition %d/%d: %s\n', c, n_cond, cond_labels{c});

    m_c3 = squeeze(mean(TEP_data{c}(c3_idx, :, :), 3));   % [1 x n_time]

    pp = TEP_peak(p_peaks, n_peaks, peak_window, t, m_c3);

    PEAK(c).condition = cond_labels{c};
    PEAK(c).peak      = pp;

    fprintf('  Components:\n');
    for k = 1:numel(pp)
        fprintf('    %-6s  lat=%5.1f ms   amp=%+.3f µV   interval=[%5.1f  %5.1f]\n', ...
            char(pp(k).label), pp(k).latency, pp(k).amplitude, ...
            pp(k).interval(1), pp(k).interval(2));
    end
    fprintf('\n');

end

% =========================================================================
% TOPOPLOTS PER COMPONENT
% =========================================================================

fprintf('--- Topoplots ---\n');

n_comp = numel(PEAK(1).peak);

for p = 1:n_comp

    comp_label = char(PEAK(1).peak(p).label);
    fprintf('  %s\n', comp_label);

    fig_topo = figure('Name', ['Topo — ' comp_label], 'Position', [100 100 1000 500]);
    k = 1;

    for c = 1:n_cond

        ga_all = TEP_data{c};
        ga     = mean(ga_all, 3);
        m_c3   = ga(c3_idx, :);
        iv     = PEAK(c).peak(p).interval;

        % Waveform panel
        subplot(n_cond, 2, k); k = k + 1;
        idx_win = (t >= peak_window(1)) & (t <= peak_window(2));
        plot(t(idx_win), m_c3(idx_win), 'k', 'LineWidth', 1.2);
        hold on; yl = ylim;
        patch([iv(1) iv(2) iv(2) iv(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.4 0.1 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(PEAK(c).peak(p).latency, PEAK(c).peak(p).amplitude, 'm.', 'MarkerSize', 14);
        xlim(peak_window); xlabel('Time (ms)'); ylabel('µV');
        title(['Grand avg — ' cond_labels{c}]);
        subtitle(['Interval: [' num2str(round(iv(1))) '  ' num2str(round(iv(2))) '] ms']);
        set(gca, 'FontSize', 11);

        % Topoplot panel
        % custum_topoplot uses ft_topoplotER with dcc_customized_acticap64.mat
        % Only common channels are passed — peripheral removed channels are absent
        subplot(n_cond, 2, k); k = k + 1;
        [~, l_idx] = min(abs(t - iv(1)));
        [~, r_idx] = min(abs(t - iv(2)));
        topo_data  = mean(mean(ga_all(:, l_idx:r_idx, :), 3), 2);

        custum_topoplot(topo_data, cellstr(chanlabs)');
        cb = colorbar; cb.Label.String = 'Amplitude (µV)'; cb.Label.FontSize = 9;
        title(cond_labels{c});

    end

    sgtitle(comp_label, 'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig_topo, fullfile(group_path, ['TEP_topo_' comp_label '.png']));

end

% =========================================================================
% SAVE TEP_results.mat
% =========================================================================

fprintf('\n--- Saving TEP_results.mat ---\n');

TEP_results.subj_list    = subj_list;
TEP_results.cond_labels  = cond_labels;
TEP_results.trig_vals    = trig_vals;
TEP_results.t            = t;              % ms [1 x 1000]
TEP_results.chanlabs     = chanlabs;       % {1 x 50} common channel set
TEP_results.c3_idx       = c3_idx;        % index of C3 in chanlabs
TEP_results.nb_idx       = nb_idx;        % Hjorth neighbour indices
TEP_results.baseline_win = baseline_win;
TEP_results.artifact_win = artifact_win;
TEP_results.TEP_data     = TEP_data;      % {1x2} each [50 x 1000 x 12]
TEP_results.PEAK         = PEAK;          % from TEP_peak.m
TEP_results.trial_data   = trial_data;    % {1x2}.(subj) [50 x 1000 x n_trials]

save(fullfile(group_path, 'TEP_results.mat'), 'TEP_results', '-v7.3');
fprintf('Saved: %s\n', fullfile(group_path, 'TEP_results.mat'));
fprintf('\n=== Done ===\n');
