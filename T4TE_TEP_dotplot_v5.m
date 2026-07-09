% =========================================================================
% T4TE — TEP Dot Plot | Trough vs. Falling | Hjorth C3
% =========================================================================
% One figure per intensity. Layout: 1 row x 4 TOIs (N45, P60, N100, P180).
% Bars = group mean per bin (trough, falling). Dots = individual subjects,
% connected across bins. Significance marker: trough vs. falling
% (reused directly from T4TE_TEP_main_analysis_v5.m — NOT recomputed here,
% so the parametric/non-parametric test choice stays consistent).
%
% Input (already computed):
%   T4TE_TEP_phase_results_v4_IAF_supra.mat
%   T4TE_TEP_phase_results_v4_IAF_sub.mat
%   T4TE_TEP_main_results_v5.mat   (for main_ttest — pre-computed stats)
%
% Output:
%   T4TE_TEP_dotplot_supra_v5.png
%   T4TE_TEP_dotplot_sub_v5.png
%
% Author:  E.W.M. Dresens | June 2026
% =========================================================================

clear; close all; clc

base_path = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/';
fig_dir   = fullfile(base_path, 'TEP_figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% -------------------------------------------------------------------------
%  SETTINGS
% -------------------------------------------------------------------------
toi_names  = {'N45','P60','N100','P180'};
rname      = 'hjorth_C3';

conditions  = {'supra','sub'};
cond_labels = {'Suprathreshold (110% rMT)','Subthreshold (90% rMT)'};
cond_files  = { ...
    'T4TE_TEP_phase_results_v4_IAF_supra.mat', ...
    'T4TE_TEP_phase_results_v4_IAF_sub.mat'};

col_trough  = [0.20 0.35 0.75];   % blue
col_falling = [0.80 0.20 0.15];   % red
col_line    = [0.65 0.65 0.65];   % grey connecting lines

fs_ax    = 13;
fs_title = 14;
fs_sub   = 10;
fs_tick  = 11;

%% -------------------------------------------------------------------------
%  LOAD SAVED RESULTS AND PRE-COMPUTED STATISTICS
% -------------------------------------------------------------------------
fprintf('Loading saved results and statistics...\n');

main_stats_path = fullfile(base_path, 'T4TE_TEP_main_results_v5.mat');
if ~exist(main_stats_path, 'file')
    error(['T4TE_TEP_main_results_v5.mat not found. Run ' ...
        'T4TE_TEP_main_analysis_v5.m first to generate main_ttest.']);
end
main_stats = load(main_stats_path, 'main_ttest');

loaded = cell(1, numel(conditions));
for ci = 1:numel(conditions)
    loaded{ci} = load(fullfile(base_path, cond_files{ci}));
    fprintf('  %s: n=%d subjects\n', cond_labels{ci}, numel(loaded{ci}.valid_subj));
end

%% =========================================================================
%  ONE FIGURE PER INTENSITY
% =========================================================================
for ci = 1:numel(conditions)
    cond = conditions{ci};
    rd   = loaded{ci};
    vs   = rd.valid_subj;
    n    = numel(vs);

    fig = figure('Name', sprintf('T4TE TEP dotplot — %s', cond), ...
        'NumberTitle','off', 'Color','w', 'Position', [50 50 1300 380]);

    for ti = 1:numel(toi_names)
        tname = toi_names{ti};
        ax = subplot(1, numel(toi_names), ti);
        hold on;

        % Extract per-subject trough and falling means at Hjorth C3
        tr_vals = nan(1, n);
        fl_vals = nan(1, n);
        for si = 1:n
            sv = vs(si);
            tr_vals(si) = rd.results(sv).toi_stats.(rname).(tname).trough;
            fl_vals(si) = rd.results(sv).toi_stats.(rname).(tname).falling;
        end

        valid = ~isnan(tr_vals) & ~isnan(fl_vals);
        tr_v  = tr_vals(valid);
        fl_v  = fl_vals(valid);
        n_v   = sum(valid);

        if n_v < 2
            axis off;
            title(tname, 'FontSize', fs_ax);
            continue
        end

        tr_mean = mean(tr_v); tr_sem = std(tr_v)/sqrt(n_v);
        fl_mean = mean(fl_v); fl_sem = std(fl_v)/sqrt(n_v);

        % Bars — one call per bin so FaceColor sets individually
        bar(1, tr_mean, 0.5, 'FaceColor', col_trough, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.55, 'HandleVisibility', 'off');
        bar(2, fl_mean, 0.5, 'FaceColor', col_falling, 'EdgeColor', 'none', ...
            'FaceAlpha', 0.55, 'HandleVisibility', 'off');

        errorbar([1 2], [tr_mean fl_mean], [tr_sem fl_sem], ...
            'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 7);

        % Individual subject dots + connecting lines
        jitter = (rand(1, n_v) - 0.5) * 0.12;
        for si = 1:n_v
            plot([1+jitter(si) 2+jitter(si)], [tr_v(si) fl_v(si)], '-', ...
                'Color', [col_line 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            plot(1+jitter(si), tr_v(si), 'o', 'Color', col_trough*0.7, ...
                'MarkerFaceColor', col_trough, 'MarkerSize', 5, 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
            plot(2+jitter(si), fl_v(si), 'o', 'Color', col_falling*0.7, ...
                'MarkerFaceColor', col_falling, 'MarkerSize', 5, 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end

        % Significance marker — reused from main_ttest (NOT recomputed)
        tt_entry = main_stats.main_ttest.(cond).(tname);
        p_val = tt_entry.p;
        d_val = tt_entry.cohens_d;

        if ~isnan(p_val)
            if p_val < 0.001;    sig = '***';
            elseif p_val < 0.01; sig = '**';
            elseif p_val < 0.05; sig = '*';
            elseif p_val < 0.10; sig = '†';
            else;                sig = '';
            end

            if ~isempty(sig)
                y_sig = max([tr_v fl_v]) + range([tr_v fl_v])*0.20;
                plot([1 2], [y_sig y_sig], '-k', 'LineWidth', 1);
                text(1.5, y_sig + range([tr_v fl_v])*0.08, sig, ...
                    'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
            end

            text(0.97, 0.05, sprintf('p=%.3f\nd=%.2f', p_val, d_val), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'FontSize', 8, 'Color', [0.35 0.35 0.35]);
        end

        yline(0, '--k', 'LineWidth', 0.7, 'Alpha', 0.5, 'HandleVisibility', 'off');
        set(ax, 'XTick', [1 2], 'XTickLabel', {'Trough','Falling'}, 'FontSize', fs_tick);
        ylabel('Amplitude (\muV)', 'FontSize', fs_ax);
        xlim([0.4 2.6]); grid on; box on;

        title(tname, 'FontSize', fs_ax, 'FontWeight', 'bold');
        subtitle(sprintf('n=%d', n_v), 'FontSize', fs_sub, 'Color', [0.4 0.4 0.4]);
    end

    sgtitle(sprintf('T4TE — TEP Phase Modulation: Trough vs. Falling | %s | Hjorth C3', ...
        cond_labels{ci}), 'FontSize', fs_title, 'FontWeight', 'bold');

    annotation('textbox', [0.01 0.01 0.98 0.03], ...
        'String', 'Blue = Trough | Red = Falling | Dots = individual subjects | Bar = group mean \pm SEM | * p<.05  ** p<.01  *** p<.001  \dag p<.10', ...
        'FontSize', 8, 'EdgeColor','none','Color',[0.45 0.45 0.45], ...
        'HorizontalAlignment','center');

    fname = sprintf('T4TE_TEP_dotplot_%s_v5.png', cond);
    exportgraphics(fig, fullfile(fig_dir, fname), 'Resolution', 150);
    fprintf('Saved: %s\n', fname);
end

fprintf('\n=== Dot plots complete ===\n');
