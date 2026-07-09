% =========================================================================
% T4TE Study — TEP Main Analysis v6 | Hjorth C3 | One-Way RM-ANOVA (4 bins)
% =========================================================================
% MAIN ANALYSIS — answers the thesis research question using a one-way
% repeated-measures ANOVA across all four phase bins (trough, rising,
% peak, falling), rather than a single pre-specified pairwise contrast.
%
% Statistics, per TOI per intensity:
%
%   (1) PRIMARY — One-way repeated-measures ANOVA (4 levels: trough,
%       rising, peak, falling) on Hjorth C3 TOI amplitude.
%       Mauchly's test of sphericity is run; if violated
%       (p < .05), the Greenhouse-Geisser corrected p-value is used
%       and reported alongside the uncorrected value.
%       If the omnibus ANOVA is significant (p < .05), Bonferroni-
%       corrected pairwise post-hoc comparisons (6 comparisons: all
%       bin pairs) are run and reported.
%
%   (2) SECONDARY / MORE SENSITIVE — Circular-linear correlation
%       (circ_corrcl, CircStat toolbox) between single-trial phase and
%       single-trial TOI amplitude, per subject per TOI per intensity.
%       Group-level test: Wilcoxon signed-rank test of the r-values
%       against zero (unchanged from v5).
%
% Inputs (already computed — Phastimate is NOT re-run here):
%   T4TE_TEP_phase_results_v4_IAF_supra.mat
%   T4TE_TEP_phase_results_v4_IAF_sub.mat
%
% Outputs:
%   T4TE_TEP_main_results_v6.mat
%   T4TE_TEP_main_waveform_v6.png   (Hjorth C3, 4 bins, both intensities)
%   Printed summary table answering the research question directly
%
% Requires: MATLAB Statistics and Machine Learning Toolbox (fitrm, ranova,
%           multcompare) and CircStat toolbox (circ_corrcl).
%
% Author:  E.W.M. Dresens
% Date:    June 2026
% =========================================================================

clear; close all; clc

%% -------------------------------------------------------------------------
%  TOOLBOX PATHS
% -------------------------------------------------------------------------
circstat_path = '/Users/e.w.m.dresens/Documents/MATLAB/CircStat2012a/';
addpath(circstat_path);

%% -------------------------------------------------------------------------
%  PATHS
% -------------------------------------------------------------------------
base_path = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/';
fig_dir   = fullfile(base_path, 'TEP_figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% -------------------------------------------------------------------------
%  SETTINGS
% -------------------------------------------------------------------------
toi_names  = {'N45','P60','N100','P180'};
alpha_sig  = 0.05;

conditions  = {'supra','sub'};
cond_labels = {'Suprathreshold (110% rMT)','Subthreshold (90% rMT)'};
cond_files  = { ...
    'T4TE_TEP_phase_results_v4_IAF_supra.mat', ...
    'T4TE_TEP_phase_results_v4_IAF_sub.mat'};

% Phase bin colours (consistent with MEP figures)
bin_names  = {'trough','rising','peak','falling'};
bin_labels = {'Trough','Rising','Peak','Falling'};
bin_cols   = {[0.20 0.35 0.75], [0.15 0.60 0.30], [0.90 0.50 0.10], [0.80 0.20 0.15]};

tep_plot_win = [-100 350];   % ms
artifact_win = [-2 15];      % ms

% Font sizes
fs_ax    = 14;
fs_title = 14;
fs_tick  = 13;
fs_sub   = 11;

%% =========================================================================
%  LOAD SAVED RESULTS (Phastimate already run — no re-extraction here)
% =========================================================================
fprintf('Loading saved Phastimate-derived results...\n');

loaded = cell(1, numel(conditions));
for ci = 1:numel(conditions)
    fpath = fullfile(base_path, cond_files{ci});
    if ~exist(fpath, 'file')
        error('File not found: %s', fpath);
    end
    loaded{ci} = load(fpath);
    fprintf('  %s: n=%d subjects\n', cond_labels{ci}, numel(loaded{ci}.valid_subj));
end

%% =========================================================================
%  STATISTIC 1 — ONE-WAY REPEATED-MEASURES ANOVA (4 bins), Hjorth C3
% =========================================================================
fprintf('\n\n=== STATISTIC 1: One-way RM-ANOVA (4 phase bins) | Hjorth C3 ===\n');

main_anova = struct();

% All pairwise bin comparisons for post-hoc (6 total)
bin_pairs = nchoosek(1:4, 2);
n_pairs   = size(bin_pairs, 1);

for ci = 1:numel(conditions)
    cond = conditions{ci};
    rd   = loaded{ci};
    vs   = rd.valid_subj;
    n    = numel(vs);

    fprintf('\n--- %s (n=%d) ---\n', cond_labels{ci}, n);

    for tt = 1:numel(toi_names)
        tname = toi_names{tt};

        % Build [n_subj x 4] data matrix: columns = trough, rising, peak, falling
        data_mat = nan(n, numel(bin_names));
        for bi = 1:numel(bin_names)
            bn = bin_names{bi};
            for si = 1:n
                sv = vs(si);
                data_mat(si, bi) = rd.results(sv).toi_stats.hjorth_C3.(tname).(bn);
            end
        end

        valid_rows = ~any(isnan(data_mat), 2);
        data_v = data_mat(valid_rows, :);
        n_v    = size(data_v, 1);

        if n_v < 4
            fprintf('  %-6s: insufficient data (n=%d)\n', tname, n_v);
            main_anova.(cond).(tname) = struct('test','none','p',NaN,'F',NaN,'n',n_v);
            continue
        end

        % --- Repeated-measures ANOVA via fitrm/ranova ---
        t_tbl = array2table(data_v, 'VariableNames', bin_names);
        within_design = table(bin_names', 'VariableNames', {'Bin'});
        within_design.Bin = categorical(within_design.Bin);

        rm = fitrm(t_tbl, 'trough-falling ~ 1', 'WithinDesign', within_design);
        ranova_tbl = ranova(rm);

        % Mauchly's test of sphericity
        mauchly_tbl = mauchly(rm);
        sphericity_p = mauchly_tbl.pValue(1);
        sphericity_violated = sphericity_p < alpha_sig;

        if sphericity_violated
            % Use Greenhouse-Geisser corrected p-value
            p_anova = ranova_tbl.pValueGG(1);
            f_anova = ranova_tbl.F(1);
            test_used = sprintf('RM-ANOVA (GG-corrected, Mauchly p=%.3f)', sphericity_p);
        else
            p_anova = ranova_tbl.pValue(1);
            f_anova = ranova_tbl.F(1);
            test_used = sprintf('RM-ANOVA (sphericity assumed, Mauchly p=%.3f)', sphericity_p);
        end

        df1 = ranova_tbl.DF(1);
        df2 = ranova_tbl.DF(2);

        if p_anova < 0.001;    sig = '***';
        elseif p_anova < 0.01; sig = '** ';
        elseif p_anova < 0.05; sig = '*  ';
        elseif p_anova < 0.10; sig = '(†)';
        else;                  sig = '   ';
        end

        fprintf('  %-6s: n=%d | %s | F(%d,%d)=%.3f, p=%.4f %s\n', ...
            tname, n_v, test_used, df1, df2, f_anova, p_anova, sig);

        % --- Post-hoc pairwise comparisons (Bonferroni), only if ANOVA significant ---
        posthoc = struct();
        if p_anova < alpha_sig
            mc_tbl = multcompare(rm, 'Bin', 'ComparisonType', 'bonferroni');

            for pp = 1:n_pairs
                b1 = bin_names{bin_pairs(pp,1)};
                b2 = bin_names{bin_pairs(pp,2)};
                pair_label = sprintf('%s_vs_%s', b1, b2);

                % Find matching row in multcompare table (either order)
                match_idx = find( ...
                    (strcmp(string(mc_tbl.Bin_1), b1) & strcmp(string(mc_tbl.Bin_2), b2)) | ...
                    (strcmp(string(mc_tbl.Bin_1), b2) & strcmp(string(mc_tbl.Bin_2), b1)), 1);

                if ~isempty(match_idx)
                    p_ph = mc_tbl.pValue(match_idx);
                    diff_ph = mc_tbl.Difference(match_idx);
                    posthoc.(pair_label) = struct('p', p_ph, 'diff', diff_ph);

                    if p_ph < 0.05
                        fprintf('      post-hoc %s vs %s: p=%.4f (Bonferroni)\n', b1, b2, p_ph);
                    end
                else
                    posthoc.(pair_label) = struct('p', NaN, 'diff', NaN);
                end
            end
        else
            for pp = 1:n_pairs
                b1 = bin_names{bin_pairs(pp,1)};
                b2 = bin_names{bin_pairs(pp,2)};
                pair_label = sprintf('%s_vs_%s', b1, b2);
                posthoc.(pair_label) = struct('p', NaN, 'diff', NaN);
            end
        end

        % Cohen's d for trough vs falling specifically (for figure annotation
        % and continuity with prior reporting), regardless of ANOVA outcome
        diffs_tf = data_v(:,1) - data_v(:,4);   % trough - falling
        cohens_d_tf = mean(diffs_tf) / std(diffs_tf);

        main_anova.(cond).(tname) = struct( ...
            'test', test_used, 'p', p_anova, 'F', f_anova, ...
            'df1', df1, 'df2', df2, 'n', n_v, ...
            'sphericity_p', sphericity_p, 'sphericity_violated', sphericity_violated, ...
            'posthoc', posthoc, 'cohens_d_trough_falling', cohens_d_tf, ...
            'data', data_v);
    end
end

%% =========================================================================
%  STATISTIC 2 — CIRCULAR-LINEAR CORRELATION, per TOI, per intensity
%  (unchanged from v5 — uses continuous phase, not the 4 discrete bins)
% =========================================================================
fprintf('\n\n=== STATISTIC 2: Circular-linear correlation | Hjorth C3 ===\n');

if exist('circ_corrcl', 'file') ~= 2
    warning(['circ_corrcl not found on path. Add CircStat toolbox path ' ...
        'at the top of this script before running.']);
end

main_corr = struct();

for ci = 1:numel(conditions)
    cond = conditions{ci};
    rd   = loaded{ci};
    vs   = rd.valid_subj;
    n    = numel(vs);

    fprintf('\n--- %s (n=%d) ---\n', cond_labels{ci}, n);

    for tt = 1:numel(toi_names)
        tname = toi_names{tt};

        r_vals = nan(1, n);
        p_vals_subj = nan(1, n);

        for si = 1:n
            sv = vs(si);
            phase_trial = rd.results(sv).trial_toi.phase;
            amp_trial   = rd.results(sv).trial_toi.(tname);

            valid_tr = ~isnan(phase_trial) & ~isnan(amp_trial);
            if sum(valid_tr) < 10
                continue
            end

            [r_cc, p_cc] = circ_corrcl(phase_trial(valid_tr)', amp_trial(valid_tr)');
            r_vals(si)      = r_cc;
            p_vals_subj(si) = p_cc;
        end

        valid_r = ~isnan(r_vals);
        r_v     = r_vals(valid_r);
        n_v     = sum(valid_r);

        if n_v < 4
            fprintf('  %-6s: insufficient data (n=%d)\n', tname, n_v);
            main_corr.(cond).(tname) = struct('mean_r',NaN,'p',NaN,'n',n_v);
            continue
        end

        [p_grp, ~, st_grp] = signrank(r_v);
        n_pos = sum(r_v > 0);

        if p_grp < 0.001;    sig = '***';
        elseif p_grp < 0.01; sig = '** ';
        elseif p_grp < 0.05; sig = '*  ';
        elseif p_grp < 0.10; sig = '(†)';
        else;                sig = '   ';
        end

        fprintf('  %-6s: n=%d | mean r=%.4f | %d/%d positive | Wilcoxon p=%.4f %s\n', ...
            tname, n_v, mean(r_v), n_pos, n_v, p_grp, sig);

        main_corr.(cond).(tname) = struct( ...
            'r_vals', r_v, 'mean_r', mean(r_v), 'n_positive', n_pos, ...
            'p', p_grp, 'n', n_v);
    end
end

%% =========================================================================
%  MAIN WAVEFORM FIGURE — Hjorth C3, ALL 4 BINS, both intensities side by side
% =========================================================================
fprintf('\n\n=== Generating main waveform figure (4 bins) ===\n');

vs1    = loaded{1}.valid_subj;
t_ms   = loaded{1}.results(vs1(1)).time;
t_mask = t_ms >= tep_plot_win(1) & t_ms <= tep_plot_win(2);
t_plot = t_ms(t_mask);

toi_centres = [45, 60, 100, 180];

fig = figure('Name', 'T4TE Main TEP — Hjorth C3 (4 bins)', 'NumberTitle','off', ...
    'Color','w', 'Position', [50 50 1100 420]);

for ci = 1:numel(conditions)
    cond = conditions{ci};
    rd   = loaded{ci};
    vs   = rd.valid_subj;
    n_s  = numel(vs);

    ax = subplot(1, 2, ci);
    hold on;

    for bi = 1:numel(bin_names)
        bn  = bin_names{bi};
        col = bin_cols{bi};

        mat = nan(n_s, sum(t_mask));
        for si = 1:n_s
            sv = vs(si);
            wave = rd.results(sv).hjorth_tep.(bn);
            if ~all(isnan(wave))
                mat(si,:) = wave(t_mask);
            end
        end

        gm = mean(mat, 1, 'omitnan');
        gs = std(mat,  0, 1, 'omitnan') / sqrt(n_s);

        fill([t_plot fliplr(t_plot)], [gm+gs fliplr(gm-gs)], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(t_plot, gm, '-', 'Color', col, 'LineWidth', 2.0, ...
            'DisplayName', bin_labels{bi});
    end

    yl = ylim;
    patch([artifact_win(1) artifact_win(2) artifact_win(2) artifact_win(1)], ...
        [yl(1)*5 yl(1)*5 yl(2)*5 yl(2)*5], [0.88 0.88 0.88], ...
        'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');

    for tt = 1:numel(toi_centres)
        xline(toi_centres(tt), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
            'Alpha', 0.6, 'HandleVisibility', 'off');

        % ANOVA significance marker per TOI
        an_entry = main_anova.(cond).(toi_names{tt});
        if ~isnan(an_entry.p) && an_entry.p < 0.05
            mark = '*';
            if an_entry.p < 0.01;  mark = '**';  end
            if an_entry.p < 0.001; mark = '***'; end
            text(toi_centres(tt), max(yl)*0.95, mark, 'FontSize', 12, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
        end
        text(toi_centres(tt), max(yl)*0.85, toi_names{tt}, 'FontSize', fs_tick-5, ...
            'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
    end
    xline(0, '--k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    yline(0, '-k', 'LineWidth', 0.5, 'Alpha', 0.4, 'HandleVisibility', 'off');

    xlabel('Time (ms)', 'FontSize', fs_ax);
    ylabel('Amplitude (\muV)', 'FontSize', fs_ax);
    title(cond_labels{ci}, 'FontSize', fs_title, 'FontWeight', 'bold');
    subtitle(sprintf('Hjorth C3 | n=%d subjects | * RM-ANOVA p<.05', n_s), ...
        'FontSize', fs_sub, 'Color', [0.4 0.4 0.4]);
    xlim(tep_plot_win); grid on; box on;
    set(ax, 'FontSize', fs_tick);

    if ci == 1
        legend('Location', 'northwest', 'FontSize', 9);
    end
end

sgtitle('T4TE — Main TEP Result: Phase-Dependent Modulation (4 Bins) at Hjorth C3', ...
    'FontSize', fs_title, 'FontWeight', 'bold');

exportgraphics(fig, fullfile(fig_dir, 'T4TE_TEP_main_waveform_v6.png'), 'Resolution', 150);
fprintf('Saved: T4TE_TEP_main_waveform_v6.png\n');

%% =========================================================================
%  SUMMARY TABLE — directly answering the research question
% =========================================================================
fprintf('\n\n=================================================================\n');
fprintf('  MAIN RESULT SUMMARY — Hjorth C3 | RM-ANOVA (4 bins) | Circ-Lin Corr\n');
fprintf('=================================================================\n');

for ci = 1:numel(conditions)
    cond = conditions{ci};
    fprintf('\n%s:\n', cond_labels{ci});
    fprintf('  %-6s  %-12s %-10s %-10s |  %-10s %-10s\n', ...
        'TOI', 'ANOVA p', 'F', 'd(T-F)', 'CircCorr p', 'mean r');
    fprintf('  %s\n', repmat('-', 1, 65));
    for tt = 1:numel(toi_names)
        tname = toi_names{tt};
        an_entry = main_anova.(cond).(tname);
        cc_entry = main_corr.(cond).(tname);
        fprintf('  %-6s  %-12.4f %-10.3f %-10.3f |  %-10.4f %-10.4f\n', ...
            tname, an_entry.p, an_entry.F, an_entry.cohens_d_trough_falling, ...
            cc_entry.p, cc_entry.mean_r);
    end
end

%% =========================================================================
%  SAVE
% =========================================================================
save(fullfile(base_path, 'T4TE_TEP_main_results_v6.mat'), ...
    'main_anova', 'main_corr', 'toi_names', 'conditions', 'cond_labels', ...
    'bin_names', 'bin_labels', 'bin_cols');

fprintf('\nSaved: T4TE_TEP_main_results_v6.mat\n');
fprintf('\n=== Main TEP analysis (RM-ANOVA, 4 bins) complete ===\n');
