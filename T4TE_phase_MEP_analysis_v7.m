%% T4TE_phase_MEP_analysis_v6.m
% Pre-stimulus mu-alpha phase vs. MEP amplitude analysis — DUAL BAND
%
% Pipeline per subject (run twice: IAF-based band and fixed 8-13 Hz):
%   1. Load EEG (_eeg_final.set) and MEP (_MEP.mat)
%   2a. IAF run:   Load _IAF_results.mat -> alpha band = IAF_peak (PAF) +/- 2 Hz
%       fixed run: alpha band = [8 13] Hz
%   3. Match valid suprathreshold trials via eeg_mask + data_emg.trialinfo
%   4. Hjorth-Laplacian on C3 (neighbours: FC1, CP1, FC5, CP5)
%   5. Phastimate: FIR bandpass + AR forward prediction + Hilbert -> phase at t=0
%   6. Z-score normalise MEP amplitudes within subject
%   7. Three phase bin analyses:
%        Analysis 1 — narrow trough vs. rest
%        Analysis 2 — narrow peak  vs. rest
%        Analysis 3 — narrow trough vs. narrow peak (direct comparison)
%   8. [NEW v6] Circular-linear correlation per subject (circ_corrcl, CircStat)
%        — trial-level phase x MEP correlation without binning
%        — group-level Wilcoxon signed-rank test on per-subject r-values
%
% Dual-band rationale:
%   The IAF run uses a subject-specific filter centred on the individual
%   alpha frequency (CoG from resting-state EEG), which is physiologically
%   more accurate. The fixed 8-13 Hz run provides a directly comparable
%   conventional reference. Both runs use identical trial sets and phase
%   extraction logic; only the FIR bandpass filter differs.
%
% v6 additions:
%   Circular-linear correlation (circ_corrcl from CircStat toolbox, Berens 2009)
%   is computed per subject on the full trial-by-trial phase x MEP z-score
%   relationship. This avoids the information loss inherent in discretising
%   phase into two bins and is more powerful for small group N. At the group
%   level, per-subject r-values are tested against zero using a one-sided
%   Wilcoxon signed-rank test (expected direction: positive r, i.e. higher
%   MEP at trough phase). Requires CircStat toolbox on the MATLAB path.
%
% Outputs:
%   T4TE_phase_MEP_results_v7_IAF.mat   — results struct, IAF run
%   T4TE_phase_MEP_results_v7_fixed.mat — results struct, fixed run
%   T4TE_phase_MEP_comparison_v6.mat    — both results in one file
%   Figures saved per band in base_path
%
% Ieke | T4TE study | CIMeC, Trento | 2026
% =========================================================================

clc; clear all; close all;

%% -------------------------------------------------------------------------
%  TOOLBOX PATHS — edit once per machine
% -------------------------------------------------------------------------
ft_path         = '/Users/e.w.m.dresens/Documents/MATLAB/fieldtrip-20250106/';
eeglab_path     = '/Users/e.w.m.dresens/Documents/MATLAB/eeglab2026.0.0/';
phastimate_path = '/Users/e.w.m.dresens/Documents/MATLAB/';
circstat_path   ='/Users/e.w.m.dresens/Documents/MATLAB/CircStat2012a/';   % Berens 2009

addpath(ft_path);   ft_defaults;
addpath(eeglab_path);   eeglab nogui;
addpath(phastimate_path);   % Phastimate.m
addpath(circstat_path);     % circ_corrcl.m etc.

% Verify CircStat is available
if ~exist('circ_corrcl', 'file')
    warning(['circ_corrcl not found. CircStat toolbox (Berens 2009) is required ' ...
             'for circular-linear correlation. Download from: ' ...
             'https://github.com/circstat/circstat-matlab\n' ...
             'Circular-linear correlation will be skipped for all subjects.']);
end

%% -------------------------------------------------------------------------
%  PARAMETERS
% -------------------------------------------------------------------------
base_path = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/';

subjects = {'BEL_S01','BEL_S02','BEL_S03','BEL_S04','BEL_S05','BEL_S06',...
            'BEL_S07','BEL_S08','BEL_S09','BEL_S10','BEL_S11','BEL_S12'};

eeg_channel       = 'C3';
hjorth_neighbours = {'FC1','CP1','FC5','CP5'};

% IAF settings
iaf_halfwidth    = 2;           % Hz either side of IAF peak (PAF)
fixed_alpha_band = [8, 13];     % Hz — fixed fallback / fixed reference band

% Phastimate parameters
fir_order      = 128;
ar_order       = 30;
edge_samples   = 64;
hilbert_window = 128;
offset_corr    = 4;

% Phase bin definitions
trough_centre = -pi;
peak_centre   =  0;
narrow_hw     = 30 * pi/180;    % +/-30 degrees half-width

% Band modes to run (both always executed)
band_modes = {'IAF', 'fixed'};

%% =========================================================================
%  MAIN LOOP: run analysis twice — once per band mode
% =========================================================================
all_results = struct();   % will hold: all_results.IAF, all_results.fixed

for bm = 1:numel(band_modes)
    mode = band_modes{bm};   % 'IAF' or 'fixed'
    fprintf('\n\n========================================================\n');
    fprintf('  BAND MODE: %s\n', mode);
    fprintf('========================================================\n');

    results = struct();

    for s = 1:numel(subjects)
        subj      = subjects{s};
        proc_path = fullfile(base_path, subj, 'processed');
        fprintf('\n=== %s [%s] ===\n', subj, mode);

        %% -- Define alpha band for this subject and mode --
        if strcmp(mode, 'IAF')
            iaf_file = fullfile(proc_path, [subj '_IAF_results.mat']);
            if exist(iaf_file, 'file')
                load(iaf_file, 'IAF_results');
                IAF        = IAF_results.IAF_peak;
                alpha_band = [IAF - iaf_halfwidth, IAF + iaf_halfwidth];
                fprintf('  Alpha band (PAF-based): [%.2f – %.2f Hz]  (PAF = %.2f Hz)\n', ...
                    alpha_band(1), alpha_band(2), IAF);
            else
                warning('%s: IAF file not found — falling back to fixed %d-%d Hz.', ...
                    subj, fixed_alpha_band(1), fixed_alpha_band(2));
                alpha_band = fixed_alpha_band;
                IAF        = NaN;
            end
        else
            alpha_band = fixed_alpha_band;
            IAF        = NaN;
            fprintf('  Alpha band (fixed): [%d – %d Hz]\n', alpha_band(1), alpha_band(2));
        end

        %% -- Load EEG --
        eeg_file = fullfile(proc_path, [subj '_eeg_final.set']);
        if ~exist(eeg_file, 'file')
            warning('%s: EEG file not found, skipping.', subj); continue
        end
        EEG = pop_loadset('filename', [subj '_eeg_final.set'], 'filepath', proc_path);

        c3_idx = find(strcmpi({EEG.chanlocs.labels}, eeg_channel));
        if isempty(c3_idx)
            warning('%s: %s not found in EEG, skipping.', subj, eeg_channel); continue
        end

        fs_eeg   = EEG.srate;
        time_eeg = EEG.times / 1000;   % ms -> s

        %% -- Load MEP and EMG --
        mep_file = fullfile(proc_path, [subj '_MEP.mat']);
        emg_file = fullfile(proc_path, [subj '_emg.mat']);
        if ~exist(mep_file, 'file')
            warning('%s: MEP file not found, skipping.', subj); continue
        end
        if ~exist(emg_file, 'file')
            warning('%s: EMG file not found (may be in iCloud), skipping.', subj); continue
        end
        load(mep_file, 'MEP', 'eeg_mask');
        load(emg_file, 'data_emg');

        %% -- Match valid suprathreshold trials --
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

        supra_eeg_idx = find(cond_eeg == 3);
        n_supra_eeg   = numel(supra_eeg_idx);

        supra_orig_idx = find(data_emg.trialinfo == 3);
        eeg_mask_supra = eeg_mask(supra_orig_idx);
        mep_eeg_valid  = MEP.suprathreshold(eeg_mask_supra);
        n_eeg_valid    = sum(eeg_mask_supra);

        if n_eeg_valid == n_supra_eeg + 1
            mep_eeg_valid = mep_eeg_valid(1:n_supra_eeg);
            warning('%s: eeg_mask gives %d supra, EEG has %d — trimming 1 trial.', ...
                subj, n_eeg_valid, n_supra_eeg);
        elseif n_eeg_valid ~= n_supra_eeg
            warning('%s: eeg_mask gives %d supra, EEG has %d epochs. Skipping.', ...
                subj, n_eeg_valid, n_supra_eeg);
            continue
        end

        mep_valid_mask2 = ~isnan(mep_eeg_valid);
        valid_ep_idx    = supra_eeg_idx(mep_valid_mask2);
        mep_vals        = mep_eeg_valid(mep_valid_mask2);
        n_valid         = sum(mep_valid_mask2);

        fprintf('%s: %d valid supra trials\n', subj, n_valid);
        if n_valid < 20
            warning('%s: fewer than 20 valid trials, results unreliable.', subj);
        end

        %% -- Phase extraction via Phastimate --
        hjorth_idx = zeros(1, numel(hjorth_neighbours));
        for hh = 1:numel(hjorth_neighbours)
            idx_h = find(strcmpi({EEG.chanlocs.labels}, hjorth_neighbours{hh}));
            if isempty(idx_h)
                warning('%s: Hjorth neighbour %s not found, using C3.', subj, hjorth_neighbours{hh});
                hjorth_idx(hh) = c3_idx;
            else
                hjorth_idx(hh) = idx_h;
            end
        end

        % Design FIR bandpass for the current band
        b_fir   = fir1(fir_order, alpha_band / (fs_eeg/2), 'bandpass');
        pre_idx = time_eeg < 0;

        phase_at_t0 = nan(1, n_valid);
        for tt = 1:n_valid
            ep  = valid_ep_idx(tt);
            c3  = double(squeeze(EEG.data(c3_idx, :, ep)));
            nbr = double(squeeze(mean(EEG.data(hjorth_idx, :, ep), 1)));
            hjorth_sig = c3 - 0.25 * nbr;
            pre_seg    = hjorth_sig(pre_idx)';

            try
                [ph, ~] = Phastimate(pre_seg, b_fir, edge_samples, ar_order, ...
                    hilbert_window, offset_corr);
                phase_at_t0(tt) = ph;
            catch ME
                warning('%s trial %d: Phastimate failed (%s)', subj, tt, ME.message);
            end
        end

        % Remove failed trials
        failed = isnan(phase_at_t0);
        if any(failed)
            fprintf('  %d trials excluded (Phastimate failed)\n', sum(failed));
            phase_at_t0(failed)  = [];
            mep_vals(failed)     = [];
            valid_ep_idx(failed) = [];
            n_valid = numel(phase_at_t0);
        end

        %% -- Z-score MEP within subject --
        mep_z = (mep_vals - mean(mep_vals,'omitnan')) / std(mep_vals,'omitnan');

        %% -- [v6] Circular-linear correlation (trial-level) --
        % circ_corrcl tests the association between a circular variable
        % (pre-stimulus phase angle in radians) and a linear variable
        % (z-scored MEP amplitude) across all valid trials for this subject,
        % without collapsing to phase bins. This uses the full trial-level
        % information and is more powerful than the binned t-test for small N.
        %
        % Expected direction: positive r (higher MEP near trough, i.e. phase ~-pi).
        % Note: circ_corrcl is unsigned — r is always positive. The group-level
        % Wilcoxon tests whether r is consistently above chance across subjects.
        if exist('circ_corrcl', 'file')
            try
                [r_circ, p_circ] = circ_corrcl(phase_at_t0', mep_z');
                fprintf('  Circ-lin correlation: r=%.3f, p=%.4f\n', r_circ, p_circ);
            catch ME
                warning('%s: circ_corrcl failed — %s', subj, ME.message);
                r_circ = NaN; p_circ = NaN;
            end
        else
            r_circ = NaN; p_circ = NaN;
        end

        %% -- Phase binning --
        ph        = angle(exp(1i * phase_at_t0));
        circ_diff = @(a,b) angle(exp(1i*(a-b)));

        trough_mask    = abs(circ_diff(ph, trough_centre)) <= narrow_hw;
        peak_mask      = abs(circ_diff(ph, peak_centre))   <= narrow_hw;
        nontrough_mask = ~trough_mask;
        nonpeak_mask   = ~peak_mask;

        n_trough    = sum(trough_mask);
        n_nontrough = sum(nontrough_mask);
        n_peak      = sum(peak_mask);
        n_nonpeak   = sum(nonpeak_mask);

        mean_mep_trough    = mean(mep_z(trough_mask),    'omitnan');
        mean_mep_nontrough = mean(mep_z(nontrough_mask), 'omitnan');
        mean_mep_peak      = mean(mep_z(peak_mask),      'omitnan');
        mean_mep_nonpeak   = mean(mep_z(nonpeak_mask),   'omitnan');
        mean_mep_tp_trough = mean_mep_trough;   % same as trough, kept for clarity
        mean_mep_tp_peak   = mean_mep_peak;

        fprintf('  Trough (n=%d): %.3f | Non-trough (n=%d): %.3f\n', ...
            n_trough, mean_mep_trough, n_nontrough, mean_mep_nontrough);
        fprintf('  Peak   (n=%d): %.3f | Non-peak   (n=%d): %.3f\n', ...
            n_peak, mean_mep_peak, n_nonpeak, mean_mep_nonpeak);

        %% -- Store per-subject results --
        results(s).subject        = subj;
        results(s).band_mode      = mode;
        results(s).IAF_peak       = IAF;            % PAF — used for filter band
        results(s).IAF_CoG        = IAF_results.IAF_CoG;  % CoG — reported descriptively
        results(s).IAF_used       = IAF;            % explicit: value used for alpha_band
        results(s).alpha_band     = alpha_band;
        results(s).n_valid        = n_valid;
        results(s).phase          = phase_at_t0;
        results(s).mep_z          = mep_z;
        results(s).trough_mask    = trough_mask;
        results(s).peak_mask      = peak_mask;
        results(s).mean_trough    = mean_mep_trough;
        results(s).mean_nontrough = mean_mep_nontrough;
        results(s).mean_peak      = mean_mep_peak;
        results(s).mean_nonpeak   = mean_mep_nonpeak;
        results(s).mean_tp_trough = mean_mep_tp_trough;
        results(s).mean_tp_peak   = mean_mep_tp_peak;
        results(s).n_trough       = n_trough;
        results(s).n_nontrough    = n_nontrough;
        results(s).n_peak         = n_peak;
        results(s).n_nonpeak      = n_nonpeak;
        results(s).r_circ         = r_circ;   % [v6] circular-linear correlation coefficient
        results(s).p_circ         = p_circ;   % [v6] corresponding p-value

        %% -- Subject-level figure --
        fig = figure('Name', sprintf('%s [%s]', subj, mode), ...
            'Position', [100 100 1200 380], 'NumberTitle', 'off');

        subplot(1,3,1);
        polarhistogram(phase_at_t0, 24, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'w');
        title(sprintf('%s [%s] — Phase (n=%d)', subj, mode, n_valid), 'FontSize', 9);

        subplot(1,3,2);
        bar([1 2], [mean_mep_trough mean_mep_nontrough], 0.5, ...
            'FaceColor', [0.4 0.6 0.8], 'EdgeColor', [0.2 0.4 0.6]);
        hold on;
        scatter(ones(1,n_trough)*1,    mep_z(trough_mask),    20, [0.6 0.6 0.6], ...
            'filled', 'jitter','on','jitterAmount',0.08);
        scatter(ones(1,n_nontrough)*2, mep_z(nontrough_mask), 20, [0.6 0.6 0.6], ...
            'filled', 'jitter','on','jitterAmount',0.08);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Trough','Rest'});
        ylabel('MEP amplitude (z-score)'); title('Narrow trough vs. rest');
        yline(0,'--k','LineWidth',0.8); grid on;

        subplot(1,3,3);
        bar([1 2], [mean_mep_tp_trough mean_mep_tp_peak], 0.5, ...
            'FaceColor', [0.8 0.5 0.4], 'EdgeColor', [0.6 0.3 0.2]);
        hold on;
        scatter(ones(1,n_trough)*1, mep_z(trough_mask), 20, [0.6 0.6 0.6], ...
            'filled', 'jitter','on','jitterAmount',0.08);
        scatter(ones(1,n_peak)*2,   mep_z(peak_mask),   20, [0.6 0.6 0.6], ...
            'filled', 'jitter','on','jitterAmount',0.08);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Trough','Peak'});
        ylabel('MEP amplitude (z-score)'); title('Trough vs. Peak');
        yline(0,'--k','LineWidth',0.8); grid on;

        sgtitle(sprintf('%s — %s band [%.1f–%.1f Hz]', subj, mode, alpha_band(1), alpha_band(2)));
        exportgraphics(fig, fullfile(proc_path, sprintf('%s_phase_MEP_%s.png', subj, mode)), ...
            'Resolution', 150);
        close(fig);
    end   % subject loop

    %% -----------------------------------------------------------------
    %  Group-level analysis for this band mode
    % -----------------------------------------------------------------
    fprintf('\n=== Group analysis [%s band] ===\n', mode);

    valid_subj    = find(~cellfun(@isempty, {results.subject}));
    trough_grp    = [results(valid_subj).mean_trough];
    nontrough_grp = [results(valid_subj).mean_nontrough];
    peak_grp      = [results(valid_subj).mean_peak];
    nonpeak_grp   = [results(valid_subj).mean_nonpeak];
    tp_trough_grp = [results(valid_subj).mean_tp_trough];
    tp_peak_grp   = [results(valid_subj).mean_tp_peak];
    n_subj        = numel(valid_subj);

    [~, p1, ci1, st1] = ttest(trough_grp, nontrough_grp);
    [~, p2, ci2, st2] = ttest(peak_grp,   nonpeak_grp);
    [~, p3, ci3, st3] = ttest(tp_trough_grp, tp_peak_grp);

    fprintf('Analysis 1 (trough vs. rest): t(%d)=%.3f, p=%.4f\n', st1.df, st1.tstat, p1);
    fprintf('Analysis 2 (peak vs. rest):   t(%d)=%.3f, p=%.4f\n', st2.df, st2.tstat, p2);
    fprintf('Analysis 3 (trough vs. peak): t(%d)=%.3f, p=%.4f\n', st3.df, st3.tstat, p3);

    %% -- [v6] Group-level circular-linear correlation --
    % Collect per-subject r and p values. Test whether r-values are
    % consistently positive across subjects using a one-sided Wilcoxon
    % signed-rank test (H1: median r > 0, i.e. phase systematically
    % predicts MEP amplitude in the expected direction).
    % A two-sided test is also reported for completeness.
    r_circ_grp = [results(valid_subj).r_circ];
    p_circ_grp = [results(valid_subj).p_circ];

    % Remove subjects where circ_corrcl failed
    valid_circ = ~isnan(r_circ_grp);
    r_circ_valid = r_circ_grp(valid_circ);

    fprintf('\n--- [v6] Group circular-linear correlation [%s band] ---\n', mode);
    fprintf('  r per subject:  '); fprintf('%.3f  ', r_circ_grp); fprintf('\n');
    fprintf('  p per subject:  '); fprintf('%.4f  ', p_circ_grp); fprintf('\n');
    fprintf('  Subjects with individual p<0.05: %d / %d\n', ...
            sum(p_circ_grp(valid_circ) < 0.05), sum(valid_circ));

    if sum(valid_circ) >= 4
        % One-sided: expected direction is r > 0
        [p_wilcox_1, ~, stats_wilcox] = signrank(r_circ_valid, 0, 'tail', 'right');
        % Two-sided for completeness
        [p_wilcox_2, ~, ~]            = signrank(r_circ_valid, 0, 'tail', 'both');
        fprintf('  Wilcoxon signed-rank (r > 0, one-sided): p=%.4f\n', p_wilcox_1);
        fprintf('  Wilcoxon signed-rank (two-sided):        p=%.4f\n', p_wilcox_2);
        fprintf('  Mean r = %.3f (SD = %.3f)\n', mean(r_circ_valid), std(r_circ_valid));
    else
        warning('Fewer than 4 subjects with valid circ_corrcl — Wilcoxon not computed.');
        p_wilcox_1 = NaN; p_wilcox_2 = NaN; stats_wilcox = [];
    end

    % Store for saving
    all_results.(mode).r_circ_grp    = r_circ_grp;
    all_results.(mode).p_circ_grp    = p_circ_grp;
    all_results.(mode).p_wilcox_1    = p_wilcox_1;
    all_results.(mode).p_wilcox_2    = p_wilcox_2;
    all_results.(mode).stats_wilcox  = stats_wilcox;

    %% -- Group figures --
    % Colour scheme: IAF = blue, fixed = orange
    if strcmp(mode, 'IAF')
        bar_col = [0.3 0.5 0.8];
        subtitle_band = sprintf('IAF +/-%d Hz', iaf_halfwidth);
    else
        bar_col = [0.85 0.55 0.25];
        subtitle_band = 'Fixed 8–13 Hz';
    end

    analyses = {
        trough_grp,    nontrough_grp,  {'Trough','Rest'},   p1, st1, 'Analysis 1: Trough vs. Rest', ...
            sprintf('T4TE_group_phase_MEP_v6_%s_analysis1.png', mode);
        peak_grp,      nonpeak_grp,    {'Peak','Rest'},     p2, st2, 'Analysis 2: Peak vs. Rest', ...
            sprintf('T4TE_group_phase_MEP_v6_%s_analysis2.png', mode);
        tp_trough_grp, tp_peak_grp,    {'Trough','Peak'},   p3, st3, 'Analysis 3: Trough vs. Peak', ...
            sprintf('T4TE_group_phase_MEP_v6_%s_analysis3.png', mode);
    };

    for panel = 1:size(analyses, 1)
        g1      = analyses{panel,1};
        g2      = analyses{panel,2};
        labs    = analyses{panel,3};
        p_pan   = analyses{panel,4};
        st_pan  = analyses{panel,5};
        ttl     = analyses{panel,6};
        fname   = analyses{panel,7};

        if p_pan < 0.001;    sig_p = '***';
        elseif p_pan < 0.01; sig_p = '**';
        elseif p_pan < 0.05; sig_p = '*';
        else;                sig_p = sprintf('p=%.3f', p_pan);
        end

        fig_grp = figure('Name', ttl, 'Position', [100 100 500 500], 'NumberTitle','off');
        hold on;
        bar([1 2], [mean(g1) mean(g2)], 0.45, ...
            'FaceColor', bar_col, 'EdgeColor', bar_col*0.7, 'FaceAlpha', 0.75);
        sem1 = std(g1)/sqrt(n_subj);
        sem2 = std(g2)/sqrt(n_subj);
        errorbar([1 2], [mean(g1) mean(g2)], [sem1 sem2], ...
            'k', 'LineWidth', 1.5, 'LineStyle', 'none', 'CapSize', 8);
        for i = 1:n_subj
            plot([1 2], [g1(i) g2(i)], 'o-', 'Color', [0.65 0.65 0.65], ...
                'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor', [0.65 0.65 0.65]);
        end
        plot(1, mean(g1), 'o', 'Color', bar_col*0.7, 'MarkerSize', 10, ...
            'LineWidth', 2, 'MarkerFaceColor', bar_col*0.7);
        plot(2, mean(g2), 'o', 'Color', bar_col*0.7, 'MarkerSize', 10, ...
            'LineWidth', 2, 'MarkerFaceColor', bar_col*0.7);

        set(gca, 'XTick', [1 2], 'XTickLabel', labs, 'FontSize', 11);
        ylabel('Mean z-scored MEP amplitude', 'FontSize', 11);
        yline(0,'--k','LineWidth',0.8); xlim([0.5 2.5]); grid on;

        y_max = max([g1 g2]) + 0.15;
        plot([1 2], [y_max y_max], '-k', 'LineWidth', 1);
        text(1.5, y_max + 0.05, sig_p, 'HorizontalAlignment', 'center', 'FontSize', 13);

        title(sprintf('%s [%s] (n=%d)', ttl, mode, n_subj), 'FontSize', 11);
        subtitle(sprintf('%s, +/-%.0f deg | t(%d)=%.2f, p=%.4f', ...
            subtitle_band, rad2deg(narrow_hw), st_pan.df, st_pan.tstat, p_pan), ...
            'FontSize', 9, 'Color', [0.4 0.4 0.4]);

        exportgraphics(fig_grp, fullfile(base_path, fname), 'Resolution', 200);
        fprintf('Saved: %s\n', fname);
        close(fig_grp);
    end

    %% -- Rose plot for this band mode --
    all_phases   = [];
    all_mep_z    = [];

    for s = valid_subj
        if isempty(results(s).subject); continue; end
        all_phases = [all_phases, results(s).phase];
        all_mep_z  = [all_mep_z,  results(s).mep_z];
    end

    n_bins    = 24;
    bin_edges = linspace(-pi, pi, n_bins + 1);

    fig_rose = figure('Name', sprintf('T4TE Rose [%s]', mode), ...
        'Position', [100 100 600 550], 'NumberTitle','off');

    subplot(1,2,1);
    polarhistogram(all_phases, n_bins, 'FaceColor', bar_col, ...
        'EdgeColor', 'w', 'FaceAlpha', 0.85);
    hold on;
    theta_trough = linspace(trough_centre-narrow_hw, trough_centre+narrow_hw, 50);
    theta_peak   = linspace(peak_centre-narrow_hw,   peak_centre+narrow_hw,   50);
    r_max = max(histcounts(all_phases, bin_edges)) * 1.3;
    polarplot([theta_trough theta_trough(end)], [r_max*ones(1,50) 0], 'r-', 'LineWidth', 2.5);
    polarplot([theta_peak   theta_peak(end)],   [r_max*ones(1,50) 0], 'g-', 'LineWidth', 2.5);
    legend({'Trials','Trough window','Peak window'}, 'Location', 'southoutside', 'FontSize', 8);
    title(sprintf('Phase distribution\n(%s, n=%d trials)', subtitle_band, numel(all_phases)), ...
        'FontSize', 10);

    subplot(1,2,2);
    bin_centres  = bin_edges(1:end-1) + diff(bin_edges)/2;
    mean_mep_bin = nan(1, n_bins);
    for b = 1:n_bins
        in_bin = all_phases >= bin_edges(b) & all_phases < bin_edges(b+1);
        if sum(in_bin) >= 3
            mean_mep_bin(b) = mean(all_mep_z(in_bin), 'omitnan');
        end
    end
    offset = abs(min(mean_mep_bin,[],'omitnan')) + 0.1;
    r_vals = mean_mep_bin + offset;
    r_vals(isnan(r_vals)) = offset;
    polarplot([bin_centres bin_centres(1)], [r_vals r_vals(1)], ...
        'o-', 'Color', [0.8 0.4 0.4], 'LineWidth', 2, ...
        'MarkerFaceColor', [0.8 0.4 0.4], 'MarkerSize', 6);
    hold on;
    polarplot(linspace(-pi,pi,200), offset*ones(1,200), '--k', 'LineWidth', 0.8);
    title(sprintf('Mean z-MEP per phase bin\n(reference ring = z-MEP of 0)'), 'FontSize', 10);

    sgtitle(sprintf('T4TE — mu-alpha phase [%s] (n=%d subjects)', subtitle_band, n_subj), ...
        'FontSize', 11);

    rose_fname = sprintf('T4TE_group_rose_plot_v6_%s.png', mode);
    exportgraphics(fig_rose, fullfile(base_path, rose_fname), 'Resolution', 200);
    fprintf('Saved: %s\n', rose_fname);
    close(fig_rose);

    %% -- Save results for this band mode --
    mat_fname = sprintf('T4TE_phase_MEP_results_v7_%s.mat', mode);
    save(fullfile(base_path, mat_fname), ...
        'results', 'trough_grp', 'nontrough_grp', 'peak_grp', 'nonpeak_grp', ...
        'tp_trough_grp', 'tp_peak_grp', ...
        'p1','p2','p3','st1','st2','st3','ci1','ci2','ci3', ...
        'mode', 'iaf_halfwidth', 'fixed_alpha_band', ...
        'r_circ_grp', 'p_circ_grp', 'p_wilcox_1', 'p_wilcox_2');   % [v6]
    fprintf('Results saved: %s\n', mat_fname);

    %% -- Aggregate into all_results struct --
    all_results.(mode).results       = results;
    all_results.(mode).trough_grp    = trough_grp;
    all_results.(mode).nontrough_grp = nontrough_grp;
    all_results.(mode).peak_grp      = peak_grp;
    all_results.(mode).nonpeak_grp   = nonpeak_grp;
    all_results.(mode).tp_trough_grp = tp_trough_grp;
    all_results.(mode).tp_peak_grp   = tp_peak_grp;
    all_results.(mode).p1 = p1; all_results.(mode).p2 = p2; all_results.(mode).p3 = p3;
    all_results.(mode).st1 = st1; all_results.(mode).st2 = st2; all_results.(mode).st3 = st3;

end   % band mode loop

%% =========================================================================
%  COMPARISON FIGURE — IAF vs. fixed side by side (Analysis 3: Trough vs. Peak)
% =========================================================================
if isfield(all_results, 'IAF') && isfield(all_results, 'fixed')
    fprintf('\n=== IAF vs. fixed comparison figure ===\n');

    n_IAF   = numel(all_results.IAF.tp_trough_grp);
    n_fixed = numel(all_results.fixed.tp_trough_grp);

    fig_cmp = figure('Name', 'T4TE — IAF vs. Fixed band comparison', ...
        'Position', [100 100 900 450], 'NumberTitle','off');

    for pp = 1:2
        if pp == 1
            r      = all_results.IAF;
            g1     = r.tp_trough_grp;
            g2     = r.tp_peak_grp;
            col    = [0.3 0.5 0.8];
            ttl    = sprintf('IAF +/-%d Hz', iaf_halfwidth);
            n_s    = n_IAF;
        else
            r      = all_results.fixed;
            g1     = r.tp_trough_grp;
            g2     = r.tp_peak_grp;
            col    = [0.85 0.55 0.25];
            ttl    = 'Fixed 8–13 Hz';
            n_s    = n_fixed;
        end

        p_val = r.p3;
        st    = r.st3;
        if p_val < 0.001;    sig_p = '***';
        elseif p_val < 0.01; sig_p = '**';
        elseif p_val < 0.05; sig_p = '*';
        else;                sig_p = sprintf('p=%.3f', p_val);
        end

        subplot(1, 2, pp);
        hold on;
        bar([1 2], [mean(g1) mean(g2)], 0.45, ...
            'FaceColor', col, 'EdgeColor', col*0.7, 'FaceAlpha', 0.75);
        errorbar([1 2], [mean(g1) mean(g2)], ...
            [std(g1)/sqrt(n_s) std(g2)/sqrt(n_s)], ...
            'k', 'LineWidth', 1.5, 'LineStyle', 'none', 'CapSize', 8);
        for i = 1:n_s
            plot([1 2], [g1(i) g2(i)], 'o-', 'Color', [0.65 0.65 0.65], ...
                'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor', [0.65 0.65 0.65]);
        end
        plot(1, mean(g1), 'o', 'Color', col*0.7, 'MarkerSize', 10, 'LineWidth', 2, ...
            'MarkerFaceColor', col*0.7);
        plot(2, mean(g2), 'o', 'Color', col*0.7, 'MarkerSize', 10, 'LineWidth', 2, ...
            'MarkerFaceColor', col*0.7);

        set(gca, 'XTick', [1 2], 'XTickLabel', {'Trough','Peak'}, 'FontSize', 11);
        ylabel('Mean z-scored MEP amplitude', 'FontSize', 11);
        yline(0,'--k','LineWidth',0.8); xlim([0.5 2.5]); grid on;

        y_max = max([g1 g2]) + 0.15;
        plot([1 2], [y_max y_max], '-k', 'LineWidth', 1);
        text(1.5, y_max+0.05, sig_p, 'HorizontalAlignment','center','FontSize',13);

        title(sprintf('Trough vs. Peak — %s', ttl), 'FontSize', 11);
        subtitle(sprintf('n=%d subjects | t(%d)=%.2f, p=%.4f', ...
            n_s, st.df, st.tstat, p_val), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
    end

    sgtitle('T4TE — Alpha phase effect: IAF-based vs. Fixed band', 'FontSize', 12);
    exportgraphics(fig_cmp, fullfile(base_path, 'T4TE_phase_MEP_IAF_vs_fixed_v6.png'), ...
        'Resolution', 200);
    fprintf('Saved: T4TE_phase_MEP_IAF_vs_fixed_v6.png\n');
    close(fig_cmp);
end

%% -- Save combined results --
save(fullfile(base_path, 'T4TE_phase_MEP_comparison_v6.mat'), 'all_results', ...
    'iaf_halfwidth', 'fixed_alpha_band', 'narrow_hw', 'trough_centre', 'peak_centre');
fprintf('\nCombined results saved: T4TE_phase_MEP_comparison_v6.mat\n');
fprintf('\n=== T4TE_phase_MEP_analysis_v6 complete ===\n');
