% =========================================================================
% T4TE Study 1.1 — Phase-Power-MEP Visualization
% =========================================================================
% Produces two figures inspired by Kirchhoff et al. (2026, Clin Neurophysiol):
%
% Figure 1 — Scatter plot: pre-stimulus phase vs. z-MEP amplitude
%   - All trials pooled across 12 subjects
%   - Dots colour-coded by z-scored IAF power (blue=low, red=high)
%   - Sinusoidal fit lines for high power (>1 SD) and low power (<-1 SD)
%   - Dotted line shows the reference mu oscillation shape
%
% Figure 2 — Boxplot: z-MEP at trough vs. peak, split by power
%   - Trough: phase within +/-30 deg of -pi/pi
%   - Peak:   phase within +/-30 deg of 0
%   - Two power groups: high (>1 SD above mean) and low (<1 SD below mean)
%   - Shows whether phase effect is power-dependent
%
% Input: BEL_S0X_prestim_power.mat per subject (from T4TE_prestim_power.m)
%
% Author:  E.W.M. Dresens
% Date:    June 2026
% =========================================================================

clear; close all; clc

% =========================================================================
%% SETTINGS
% =========================================================================
base_path = 
subjects = 

% Phase bin definitions (radians)
trough_centre = -pi;
peak_centre   =  0;
bin_hw        = 30 * pi/180;   % +/- 30 degrees

% Power split threshold (z-score)
high_mask = all_power_z >= median(all_power_z);
low_mask  = all_power_z <  median(all_power_z);

% Colours
col_high = [0.75 0.15 0.10];   % dark red  — high power
col_low  = [0.12 0.35 0.75];   % dark blue — low power

% =========================================================================
%% LOAD AND POOL DATA
% =========================================================================
all_phase   = [];
all_mep_z   = [];
all_power_z = [];   % z-scored IAF power (within subject, then pooled)

for s = 1:numel(subjects)
    subj     = subjects{s};
    pwr_file = fullfile(base_path, subj, 'processed', ...
                        sprintf('%s_prestim_power.mat', subj));
    if ~exist(pwr_file, 'file')
        warning('%s: power file not found, skipping.', subj);
        continue
    end
    load(pwr_file, 'phase_v7', 'mep_z_v7', 'power_IAF_z');

    all_phase   = [all_phase;   phase_v7(:)];
    all_mep_z   = [all_mep_z;   mep_z_v7(:)];
    all_power_z = [all_power_z; power_IAF_z(:)];
end

n_total = numel(all_phase);
fprintf('Total trials pooled: %d\n', n_total);

% =========================================================================
%% SINUSOIDAL FIT FUNCTION
% =========================================================================
% Fit: y = A*sin(x + phi) + C  using least squares
sin_fit = @(p, x) p(1)*sin(x + p(2)) + p(3);

% =========================================================================
%% FIGURE 1 — SCATTER PLOT WITH SINUSOIDAL FITS
% =========================================================================
fig1 = figure('Name', 'T4TE Phase-Power-MEP Scatter', ...
    'Position', [50 50 820 500], 'NumberTitle', 'off');

ax = axes('Parent', fig1);
hold(ax, 'on');

% --- Scatter: colour by power_z ---
% Use a diverging colormap (blue->white->red)
scatter(ax, all_phase, all_mep_z, 8, all_power_z, ...
    'filled', 'MarkerFaceAlpha', 0.35);
colormap(ax, brewermap(256, '*RdBu'));   % requires brewermap toolbox
clim(ax, [-2.5 2.5]);
cb = colorbar(ax);
cb.Label.String = 'Pre-stimulus \mu-power (z-score)';
cb.Label.FontSize = 10;

% --- Reference oscillation (dotted) ---
x_ref  = linspace(-pi, pi, 300);
y_ref  = -0.4 * cos(x_ref) - 1.5;   % scaled cosine, shifted below data
plot(ax, x_ref, y_ref, ':k', 'LineWidth', 1.5);
text(ax, pi - 0.1, y_ref(end) + 0.08, '\mu-osc.', ...
    'HorizontalAlignment','right', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

% --- Sinusoidal fits for high and low power ---
high_mask = all_power_z >= median(all_power_z);
low_mask  = all_power_z <  median(all_power_z);

phase_fit = linspace(-pi, pi, 300);

for grp = 1:2
    if grp == 1
        mask = high_mask;
        col  = col_high;
        lbl  = sprintf('High power (> +%d SD, n=%d)', power_thresh, sum(mask));
    else
        mask = low_mask;
        col  = col_low;
        lbl  = sprintf('Low power (< -%d SD, n=%d)', power_thresh, sum(mask));
    end

    if sum(mask) < 10
        warning('Too few trials for sinusoidal fit (group %d).', grp);
        continue
    end

    ph_grp  = all_phase(mask);
    mep_grp = all_mep_z(mask);

    % Initial parameter guess: [amplitude, phase_offset, offset]
    p0 = [0.1, 0, mean(mep_grp)];
    opts = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000);
    try
        p_fit = fminsearch(@(p) sum((sin_fit(p, ph_grp) - mep_grp).^2), p0, opts);
        y_fit = sin_fit(p_fit, phase_fit);
        plot(ax, phase_fit, y_fit, '-', 'Color', col, 'LineWidth', 2.5, ...
            'DisplayName', lbl);
    catch
        warning('Sinusoidal fit failed for group %d.', grp);
    end
end

% --- Formatting ---
xline(ax, trough_centre, ':', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5);
xline(ax, peak_centre,   ':', 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
text(ax, trough_centre, ax.YLim(2)*0.9, 'Trough', ...
    'Color', [0.2 0.2 0.8], 'FontSize', 9, 'HorizontalAlignment', 'center');
text(ax, peak_centre, ax.YLim(2)*0.9, 'Peak', ...
    'Color', [0.1 0.6 0.1], 'FontSize', 9, 'HorizontalAlignment', 'center');

set(ax, 'XTick', [-pi -pi/2 0 pi/2 pi], ...
    'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'}, 'FontSize', 10);
xlabel(ax, 'Pre-stimulus \mu-alpha phase (rad)', 'FontSize', 11);
ylabel(ax, 'z-scored MEP amplitude', 'FontSize', 11);
title(ax, 'T4TE — Pre-stimulus \mu-alpha phase vs. MEP amplitude', 'FontSize', 11);
subtitle(ax, sprintf('All trials pooled (n=%d trials, n=12 subjects) | colour = pre-stimulus \mu-power', ...
    n_total), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
legend(ax, 'Location', 'northeast', 'FontSize', 9);
xlim(ax, [-pi pi]);
ylim(ax, [-3 3]);
grid(ax, 'on'); box(ax, 'on');

exportgraphics(fig1, fullfile(base_path, 'T4TE_phase_power_MEP_scatter.png'), ...
    'Resolution', 250);
fprintf('Saved: T4TE_phase_power_MEP_scatter.png\n');

% =========================================================================
%% FIGURE 2 — BOXPLOT: TROUGH vs. PEAK × POWER
% =========================================================================
fig2 = figure('Name', 'T4TE Phase-Power Boxplot', ...
    'Position', [50 600 600 420], 'NumberTitle', 'off');

ax2 = axes('Parent', fig2);
hold(ax2, 'on');

% Define phase bins
circ_diff   = @(a, b) angle(exp(1i * (a - b)));
trough_mask = abs(circ_diff(all_phase, trough_centre)) <= bin_hw;
peak_mask   = abs(circ_diff(all_phase, peak_centre))   <= bin_hw;

% Four groups: [trough-low, trough-high, peak-low, peak-high]
groups = {
    all_mep_z(trough_mask & low_mask),  ...   % trough, low power
    all_mep_z(trough_mask & high_mask), ...   % trough, high power
    all_mep_z(peak_mask   & low_mask),  ...   % peak, low power
    all_mep_z(peak_mask   & high_mask)  ...   % peak, high power
};

ns = cellfun(@numel, groups);
fprintf('\nGroup sizes:\n');
fprintf('  Trough-Low: %d  |  Trough-High: %d\n', ns(1), ns(2));
fprintf('  Peak-Low:   %d  |  Peak-High:   %d\n', ns(3), ns(4));

% X positions
x_pos  = [1, 1.5, 2.5, 3];
cols   = {col_low, col_high, col_low, col_high};
labels = {'Low','High','Low','High'};

for g = 1:4
    data_g = groups{g};
    if isempty(data_g); continue; end

    % Box
    q1  = prctile(data_g, 25);
    med = median(data_g);
    q3  = prctile(data_g, 75);
    iqr_g = q3 - q1;
    w_low  = max(data_g(data_g >= q1 - 1.5*iqr_g));
    w_high = min(data_g(data_g <= q3 + 1.5*iqr_g));

    % Draw box
    x  = x_pos(g);
    bw = 0.18;   % box half-width
    fill(ax2, [x-bw x+bw x+bw x-bw x-bw], [q1 q1 q3 q3 q1], ...
        cols{g}, 'FaceAlpha', 0.45, 'EdgeColor', cols{g}, 'LineWidth', 1.2);
    plot(ax2, [x-bw x+bw], [med med], '-', 'Color', cols{g}*0.6, 'LineWidth', 2);

    % Whiskers
    plot(ax2, [x x], [w_low q1],  '-', 'Color', cols{g}, 'LineWidth', 1);
    plot(ax2, [x x], [q3 w_high], '-', 'Color', cols{g}, 'LineWidth', 1);
    plot(ax2, [x-bw/2 x+bw/2], [w_low  w_low],  '-', 'Color', cols{g}, 'LineWidth', 1);
    plot(ax2, [x-bw/2 x+bw/2], [w_high w_high], '-', 'Color', cols{g}, 'LineWidth', 1);

    % Outliers
    outliers = data_g(data_g < w_low | data_g > w_high);
    if ~isempty(outliers)
        scatter(ax2, x*ones(size(outliers)), outliers, 15, '+', ...
            'MarkerEdgeColor', cols{g}, 'LineWidth', 0.8);
    end

    % n label
    text(ax2, x, w_high + 0.08, sprintf('n=%d', ns(g)), ...
        'HorizontalAlignment','center', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
end

% Reference line
yline(ax2, 0, '--k', 'LineWidth', 0.8, 'Alpha', 0.5);

% Group divider
xline(ax2, 2, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

% X labels
set(ax2, 'XTick', [1 1.5 2.5 3], ...
    'XTickLabel', {'Low','High','Low','High'}, 'FontSize', 10);

% Group labels below x-axis
text(ax2, 1.25, ax2.YLim(1) - 0.25, 'Trough (±30°)', ...
    'HorizontalAlignment','center', 'FontSize', 10, 'FontWeight','bold', ...
    'Color', [0.2 0.2 0.7]);
text(ax2, 2.75, ax2.YLim(1) - 0.25, 'Peak (±30°)', ...
    'HorizontalAlignment','center', 'FontSize', 10, 'FontWeight','bold', ...
    'Color', [0.1 0.6 0.1]);

% Legend patches
patch(ax2, NaN, NaN, col_low,  'FaceAlpha', 0.5, 'DisplayName', 'Low \mu-power (< -1 SD)');
patch(ax2, NaN, NaN, col_high, 'FaceAlpha', 0.5, 'DisplayName', 'High \mu-power (> +1 SD)');
legend(ax2, 'Location', 'northeast', 'FontSize', 9);

ylabel(ax2, 'z-scored MEP amplitude', 'FontSize', 11);
title(ax2, 'T4TE — MEP by phase bin and pre-stimulus \mu-power', 'FontSize', 11);
subtitle(ax2, 'All trials pooled (n=12 subjects) | trough = ±\pi ±30°, peak = 0 ±30°', ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);
grid(ax2, 'on'); box(ax2, 'on');

exportgraphics(fig2, fullfile(base_path, 'T4TE_phase_power_MEP_boxplot.png'), ...
    'Resolution', 250);
fprintf('Saved: T4TE_phase_power_MEP_boxplot.png\n');
fprintf('\n=== Visualization complete ===\n');
