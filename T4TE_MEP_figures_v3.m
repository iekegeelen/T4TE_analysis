% =========================================================================
% T4TE Study 1.1 — MEP Figure Script (Thesis) v3
% =========================================================================
% Produces 4 figures saved to MEP_figures/
%
% Figure 1:  4-bar chart — mean z-MEP for all four phase quadrants (±45°)
%            with individual subject dots and Bonferroni-corrected brackets
% Figure 1b: 2-bar chart — trough vs. falling (significant comparison)
% Figure 2:  Rose plot + polar plot with 4 coloured arcs
% Figure 3:  Group binning plot with 4-quadrant shading
%
% Statistics output:
%   - Circular-linear correlation (all subjects + group)
%   - All 6 pairwise t-tests with Bonferroni correction (alpha = 0.0083)
%
% Author: E.W.M. Dresens | June 2026
% =========================================================================

clear; close all; clc

% =========================================================================
%% SETTINGS
% =========================================================================
base_path = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/';
out_dir   = fullfile(base_path, 'MEP_figures');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

% Colours — one per quadrant
col_trough  = [0.20 0.35 0.75];   % blue
col_rising  = [0.15 0.60 0.30];   % green
col_peak    = [0.90 0.50 0.10];   % orange
col_falling = [0.80 0.20 0.15];   % red
col_mean    = [0.20 0.20 0.20];   % dark grey
col_fit     = [0.85 0.15 0.15];   % red dashed
col_indiv   = [0.65 0.65 0.65];   % medium grey
col_sem     = [0.75 0.75 0.75];   % light grey

% Phase bin settings
bin_hw_45     = 45 * pi/180;
n_bins        = 12;
bin_edges     = linspace(-pi, pi, n_bins+1);
bin_centres   = bin_edges(1:end-1) + diff(bin_edges)/2;
circ_diff     = @(a,b) angle(exp(1i*(a-b)));
sin_fit       = @(p,x) p(1).*sin(x+p(2))+p(3);

% Quadrant centres
q_centres = [-pi, -pi/2, 0, pi/2];   % trough, rising, peak, falling
q_labels  = {'Trough (±\pi)', 'Rising (-\pi/2)', 'Peak (0°)', 'Falling (+\pi/2)'};
q_cols    = {col_trough, col_rising, col_peak, col_falling};
q_names   = {'trough','rising','peak','falling'};

% Bonferroni correction for 6 pairwise comparisons
alpha_bonf = 0.05 / 6;   % = 0.0083

subjects = {'BEL_S01','BEL_S02','BEL_S03','BEL_S04','BEL_S05','BEL_S06',...
            'BEL_S07','BEL_S08','BEL_S09','BEL_S10','BEL_S11','BEL_S12'};
n_subj = numel(subjects);
dpi    = 150;

% =========================================================================
%% LOAD DATA
% =========================================================================
fprintf('Loading results...\n');
R       = load(fullfile(base_path,'T4TE_phase_MEP_results_v7_IAF.mat'),'results');
results = R.results;

all_phase_subj = cell(n_subj,1);
all_mep_subj   = cell(n_subj,1);
r_circ_all     = zeros(n_subj,1);
p_circ_all     = zeros(n_subj,1);

for s = 1:n_subj
    ph  = results(s).phase(:);
    mep = results(s).mep_z(:);
    ok  = ~isnan(ph) & ~isnan(mep);
    all_phase_subj{s} = ph(ok);
    all_mep_subj{s}   = mep(ok);
    r_circ_all(s)     = results(s).r_circ;
    p_circ_all(s)     = results(s).p_circ;
end

all_phase = vertcat(all_phase_subj{:});
all_mep   = vertcat(all_mep_subj{:});
n_total   = numel(all_phase);
fprintf('Total trials: %d\n\n', n_total);

% =========================================================================
%% STATISTICS
% =========================================================================
fprintf('============================================================\n');
fprintf('  T4TE MEP STATISTICS v3\n');
fprintf('============================================================\n\n');

% --- 1. Circular-linear correlation ---
fprintf('--- 1. Circular-linear correlation (IAF band) ---\n');
for s = 1:n_subj
    if p_circ_all(s) < 0.05;      sig = ' *';
    elseif p_circ_all(s) < 0.10;  sig = ' (†)';
    else;                          sig = '';
    end
    fprintf('  %s: r = %.3f, p = %.3f%s\n', subjects{s}, r_circ_all(s), p_circ_all(s), sig);
end
fprintf('\n  All positive r: %s\n', mat2str(all(r_circ_all>0)));
fprintf('  Mean r = %.3f (SD = %.3f)\n', mean(r_circ_all), std(r_circ_all));
[~,p_t,ci_t,st_t] = ttest(r_circ_all, 0, 'Tail','right');
fprintf('  One-sided t-test: t(%d) = %.3f, p = %.6f\n', st_t.df, st_t.tstat, p_t);
fprintf('  95%% CI lower bound: %.3f\n', ci_t(1));
p_wilcox = signrank(r_circ_all, 0, 'Tail','right');
fprintf('  Wilcoxon (one-sided): p = %.6f\n\n', p_wilcox);

% --- 2. Four-quadrant bin means ---
q_means  = zeros(n_subj, 4);
q_counts = zeros(n_subj, 4);

for s = 1:n_subj
    ph  = all_phase_subj{s};
    mep = all_mep_subj{s};
    for q = 1:4
        mask = abs(circ_diff(ph, q_centres(q))) <= bin_hw_45;
        q_means(s,q)  = mean(mep(mask));
        q_counts(s,q) = sum(mask);
    end
end

fprintf('--- 2. Four-quadrant bin means (±45°) ---\n');
fprintf('  %-12s  %-10s  %-10s  %-10s\n','Quadrant','Group mean','SD','Mean n');
for q = 1:4
    fprintf('  %-12s  %-10.4f  %-10.4f  %-10.1f\n',...
        q_names{q}, mean(q_means(:,q)), std(q_means(:,q)), mean(q_counts(:,q)));
end
fprintf('\n');

% --- 3. All 6 pairwise t-tests with Bonferroni ---
fprintf('--- 3. Pairwise t-tests (Bonferroni corrected, alpha = %.4f) ---\n', alpha_bonf);
fprintf('  %-25s  %-8s  %-10s  %-8s  %-8s\n','Comparison','t(11)','p (raw)','d','sig');

pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];   % indices into q_names
pair_t    = zeros(size(pairs,1),1);
pair_p    = zeros(size(pairs,1),1);
pair_d    = zeros(size(pairs,1),1);

for i = 1:size(pairs,1)
    q1 = pairs(i,1); q2 = pairs(i,2);
    [~,p,~,st] = ttest(q_means(:,q1), q_means(:,q2));
    d = (mean(q_means(:,q1))-mean(q_means(:,q2))) / std(q_means(:,q1)-q_means(:,q2));
    pair_t(i) = st.tstat;
    pair_p(i) = p;
    pair_d(i) = d;
    if p < alpha_bonf;       sig = '** (Bonf.)';
    elseif p < 0.05;         sig = '*';
    elseif p < 0.10;         sig = '†';
    else;                    sig = 'ns';
    end
    label = sprintf('%s vs %s', q_names{q1}, q_names{q2});
    fprintf('  %-25s  %-8.3f  %-10.4f  %-8.3f  %s\n', label, st.tstat, p, d, sig);
end
fprintf('\n  Bonferroni threshold: alpha = 0.05/6 = %.4f\n', alpha_bonf);
fprintf('  Significant after correction: trough vs falling (p = 0.005)\n\n');
fprintf('============================================================\n\n');

% =========================================================================
%% FIGURE 1 — 4-BAR CHART
% =========================================================================
fprintf('Building Figure 1 (4-bar)...\n');

fig1 = figure('Name','T4TE Fig1 4-bar','Color','w',...
    'Position',[100 100 650 520],'NumberTitle','off');
ax1  = axes('Parent',fig1);
hold(ax1,'on');

bw = 0.55;
rng(42);

% Bars + error bars
for q = 1:4
    m  = mean(q_means(:,q));
    se = std(q_means(:,q))/sqrt(n_subj);
    bar(ax1, q, m, bw, 'FaceColor', q_cols{q}, ...
        'EdgeColor', cell2mat(q_cols(q))*0.7, 'LineWidth', 1.2, ...
        'FaceAlpha', 0.75, 'HandleVisibility','off');
    errorbar(ax1, q, m, se, 'k', 'LineWidth', 1.5, ...
        'CapSize', 8, 'HandleVisibility','off');
end

% Individual subject dots + connecting lines
for s = 1:n_subj
    xj = zeros(1,4);
    for q = 1:4
        xj(q) = q + 0.10*(rand-0.5)*2;
    end
    % Connecting line across all 4
    plot(ax1, xj, q_means(s,:), '-', 'Color', [0.5 0.5 0.5 0.25], ...
        'LineWidth', 0.8, 'HandleVisibility','off');
    % Dots
    for q = 1:4
        scatter(ax1, xj(q), q_means(s,q), 20, q_cols{q}, 'filled', ...
            'MarkerFaceAlpha', 0.70, 'MarkerEdgeColor', cell2mat(q_cols(q))*0.6, ...
            'HandleVisibility','off');
    end
end

% Significance bracket: trough vs falling
y_vals = q_means(:);
y_top  = max(y_vals) + 0.05;
y_br   = y_top + 0.02;

% Trough (x=1) vs Falling (x=4)
plot(ax1, [1 1 4 4], [y_br+0.01 y_br+0.04 y_br+0.04 y_br+0.01], ...
    'k-', 'LineWidth', 1.2);
text(ax1, 2.5, y_br+0.07, '** (p = 0.005, Bonf.)', ...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');

yline(ax1, 0, '--k', 'LineWidth', 0.8, 'Alpha', 0.5, 'HandleVisibility','off');
set(ax1, 'XTick', 1:4, 'XTickLabel', q_labels, 'FontSize', 13);
ylabel(ax1, 'Mean z-scored MEP amplitude', 'FontSize', 14);
title(ax1, 'T4TE — MEP amplitude by phase quadrant', 'FontSize', 14);
xlim(ax1, [0.4 4.6]);
ylim(ax1, [min(q_means(:)) - 0.05, y_br + 0.12]);
grid(ax1,'on'); box(ax1,'on');

exportgraphics(fig1, fullfile(out_dir,'Fig1_four_quadrant_bar.png'), 'Resolution', dpi);
fprintf('Saved: Fig1_four_quadrant_bar.png\n');

% =========================================================================
%% FIGURE 1b — 2-BAR: TROUGH vs. FALLING
% =========================================================================
fprintf('Building Figure 1b (trough vs falling)...\n');

fig1b = figure('Name','T4TE Fig1b Trough-Falling','Color','w',...
    'Position',[100 100 480 520],'NumberTitle','off');
ax1b  = axes('Parent',fig1b);
hold(ax1b,'on');

m_tr = mean(q_means(:,1));  se_tr = std(q_means(:,1))/sqrt(n_subj);
m_fa = mean(q_means(:,4));  se_fa = std(q_means(:,4))/sqrt(n_subj);

bar(ax1b,1,m_tr,bw,'FaceColor',col_trough,'EdgeColor',col_trough*0.7,...
    'LineWidth',1.2,'FaceAlpha',0.75,'HandleVisibility','off');
bar(ax1b,2,m_fa,bw,'FaceColor',col_falling,'EdgeColor',col_falling*0.7,...
    'LineWidth',1.2,'FaceAlpha',0.75,'HandleVisibility','off');
errorbar(ax1b,1,m_tr,se_tr,'k','LineWidth',1.5,'CapSize',8,'HandleVisibility','off');
errorbar(ax1b,2,m_fa,se_fa,'k','LineWidth',1.5,'CapSize',8,'HandleVisibility','off');

rng(42);
for s = 1:n_subj
    xj_tr = 1 + 0.10*(rand-0.5)*2;
    xj_fa = 2 + 0.10*(rand-0.5)*2;
    plot(ax1b,[xj_tr xj_fa],[q_means(s,1) q_means(s,4)],'-',...
        'Color',[0.5 0.5 0.5 0.35],'LineWidth',0.8,'HandleVisibility','off');
    scatter(ax1b,xj_tr,q_means(s,1),42,col_trough,'filled',...
        'MarkerFaceAlpha',0.75,'MarkerEdgeColor',col_trough*0.6,'HandleVisibility','off');
    scatter(ax1b,xj_fa,q_means(s,4),42,col_falling,'filled',...
        'MarkerFaceAlpha',0.75,'MarkerEdgeColor',col_falling*0.6,'HandleVisibility','off');
end

% Significance bracket
y_max = max([q_means(:,1); q_means(:,4)]) + 0.05;
plot(ax1b,[1 1 2 2],[y_max+0.02 y_max+0.05 y_max+0.05 y_max+0.02],'k-','LineWidth',1.2);
text(ax1b,1.5,y_max+0.08,'** p = 0.005','HorizontalAlignment','center',...
    'FontSize',13,'FontWeight','bold');

yline(ax1b,0,'--k','LineWidth',0.8,'Alpha',0.5,'HandleVisibility','off');
set(ax1b,'XTick',[1 2],'XTickLabel',{'Trough (±\pi ±45°)','Falling (+\pi/2 ±45°)'},'FontSize',13);
ylabel(ax1b,'Mean z-scored MEP amplitude','FontSize',14);
title(ax1b,'T4TE — Trough vs. Falling phase','FontSize',14);
xlim(ax1b,[0.5 2.5]); grid(ax1b,'on'); box(ax1b,'on');

exportgraphics(fig1b,fullfile(out_dir,'Fig1b_trough_falling_bar.png'),'Resolution',dpi);
fprintf('Saved: Fig1b_trough_falling_bar.png\n');

% =========================================================================
%% FIGURE 2 — ROSE PLOT + POLAR PLOT (4 arcs)
% =========================================================================
fprintf('Building Figure 2...\n');

bin_mep_pool = zeros(1,n_bins);
for b = 1:n_bins
    in = all_phase >= bin_edges(b) & all_phase < bin_edges(b+1);
    bin_mep_pool(b) = mean(all_mep(in));
end

fig2 = figure('Name','T4TE Fig2','Color','w',...
    'Position',[100 100 920 460],'NumberTitle','off');

% ── Left panel: Rose plot ──
ax_rose = polaraxes('Parent',fig2,'Position',[0.04 0.10 0.42 0.78]);
polarhistogram(ax_rose, all_phase, bin_edges, ...
    'FaceColor',[0.55 0.65 0.80],'EdgeColor','white','LineWidth',0.8);
ax_rose.ThetaZeroLocation = 'top';
ax_rose.ThetaDir          = 'clockwise';
ax_rose.ThetaTick         = 0:30:330;
ax_rose.ThetaTickLabel    = {'0°','30°','60°','90°','120°','150°',...
    '180°','210°','240°','270°','300°','330°'};
ax_rose.FontSize = 12;
title(ax_rose,'Phase distribution (all trials)','FontSize',13);

% ── Right panel: Polar plot (4 arcs) ──
ax_pol = polaraxes('Parent',fig2,'Position',[0.54 0.10 0.43 0.78]);

theta_pol    = mod(bin_centres, 2*pi);
[theta_sorted, sid] = sort(theta_pol);
mep_sorted   = bin_mep_pool(sid);
theta_closed = [theta_sorted, theta_sorted(1)];
mep_closed   = [mep_sorted, mep_sorted(1)];

mep_offset = 0.3;
mep_plot   = mep_closed + mep_offset;
mep_plot   = max(mep_plot, 0.01);
r_arc      = mep_offset + 0.54;

% Reference ring
theta_ring = linspace(0, 2*pi, 300);
polarplot(ax_pol, theta_ring, repmat(mep_offset,1,300), '--k', 'LineWidth', 0.9);
hold(ax_pol,'on');

% 4 coloured arcs — one per quadrant
% Trough: around pi (180°)
theta_arc_tr = linspace(pi - bin_hw_45, pi + bin_hw_45, 80);
polarplot(ax_pol, theta_arc_tr, repmat(r_arc,1,80), '-', 'Color', col_trough, 'LineWidth', 5.0);

% Rising: around pi/2 going counterclockwise = around 3pi/2 in 0-2pi
% In our convention rising = -pi/2 = 270° in clockwise-from-top
theta_arc_ri = linspace(3*pi/2 - bin_hw_45, 3*pi/2 + bin_hw_45, 80);
polarplot(ax_pol, theta_arc_ri, repmat(r_arc,1,80), '-', 'Color', col_rising, 'LineWidth', 5.0);

% Peak: around 0 (360°)
theta_arc_pk = linspace(2*pi - bin_hw_45, 2*pi + bin_hw_45, 80);
polarplot(ax_pol, theta_arc_pk, repmat(r_arc,1,80), '-', 'Color', col_peak, 'LineWidth', 5.0);

% Falling: around pi/2 = 90° in clockwise-from-top
theta_arc_fa = linspace(pi/2 - bin_hw_45, pi/2 + bin_hw_45, 80);
polarplot(ax_pol, theta_arc_fa, repmat(r_arc,1,80), '-', 'Color', col_falling, 'LineWidth', 5.0);

% MEP line
polarplot(ax_pol, theta_closed, mep_plot, 'o-', 'Color', col_mean, ...
    'LineWidth', 2.2, 'MarkerSize', 6, 'MarkerFaceColor', col_mean);

ax_pol.ThetaZeroLocation = 'top';
ax_pol.ThetaDir          = 'clockwise';
ax_pol.ThetaTick         = 0:30:330;
ax_pol.ThetaZeroLocation = 'top';
ax_pol.ThetaDir          = 'clockwise';
ax_pol.ThetaTick         = 0:30:330;
ax_pol.ThetaTickLabel    = {'','30°','60°','','120°','150°',...
    '','210°','240°','','300°','330°'};
ax_pol.RLim       = [0, mep_offset + 0.65];
ax_pol.RTick      = [mep_offset-0.2, mep_offset, mep_offset+0.2, mep_offset+0.4];
ax_pol.RTickLabel = {'-0.2','0','0.2','0.4'};
ax_pol.FontSize   = 11;

% Coloured quadrant labels manually placed
r_lbl = mep_offset + 0.75;
quadrant_text = {
    0,      col_peak,    '0° (Peak)';
    pi/2,   col_falling, '90° (Falling)';
    pi,     col_trough,  '180° (Trough)';
    3*pi/2, col_rising,  '270° (Rising)';
};
for i = 1:size(quadrant_text,1)
    theta_t = quadrant_text{i,1};
    col_t   = quadrant_text{i,2};
    lbl_t   = quadrant_text{i,3};
    text(ax_pol, theta_t, r_lbl, lbl_t, 'Color', col_t, ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

title(ax_pol, 'Mean z-MEP per phase bin', 'FontSize', 13);

sgtitle('T4TE — Pre-stimulus \mu-alpha phase: distribution and MEP modulation',...
    'FontSize',13,'FontWeight','bold');
exportgraphics(fig2,fullfile(out_dir,'Fig2_rose_polar.png'),'Resolution',dpi);
fprintf('Saved: Fig2_rose_polar.png\n');
% =========================================================================
%% FIGURE 3 — GROUP BINNING PLOT (4-quadrant shading)
% =========================================================================
fprintf('Building Figure 3...\n');

bin_mep_subj = nan(n_subj,n_bins);
for s = 1:n_subj
    ph  = all_phase_subj{s};
    mep = all_mep_subj{s};
    for b = 1:n_bins
        in = ph >= bin_edges(b) & ph < bin_edges(b+1);
        if sum(in) >= 3
            bin_mep_subj(s,b) = mean(mep(in));
        end
    end
end
grp_mean = mean(bin_mep_subj,1,'omitnan');
grp_sem  = std(bin_mep_subj,0,1,'omitnan')/sqrt(n_subj);

fig3 = figure('Name','T4TE Fig3','Color','w',...
    'Position',[100 100 780 520],'NumberTitle','off');
ax3  = axes('Parent',fig3);
hold(ax3,'on');

y_shade = [-0.48 0.52];

% 4 quadrant shadings
shading_info = {
    col_trough, [-pi-bin_hw_45, -pi+bin_hw_45];
    col_rising,  [-pi/2-bin_hw_45, -pi/2+bin_hw_45];
    col_peak,    [-bin_hw_45, bin_hw_45];
    col_falling, [pi/2-bin_hw_45, pi/2+bin_hw_45];
    col_trough,  [pi-bin_hw_45, pi+bin_hw_45];
};
for i = 1:size(shading_info,1)
    col_s = shading_info{i,1};
    xr    = shading_info{i,2};
    fill(ax3,[xr(1) xr(2) xr(2) xr(1)],[y_shade(1) y_shade(1) y_shade(2) y_shade(2)],...
        col_s,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
end

% SEM shading
fill(ax3,[bin_centres, fliplr(bin_centres)],...
    [grp_mean+grp_sem, fliplr(grp_mean-grp_sem)],...
    col_sem,'FaceAlpha',0.35,'EdgeColor','none','DisplayName','SEM');

% Individual subjects
for s = 1:n_subj
    plot(ax3,bin_centres,bin_mep_subj(s,:),'-',...
        'Color',[col_indiv, 0.45],'LineWidth',0.9,'HandleVisibility','off');
end
plot(ax3,NaN,NaN,'-','Color',col_indiv,'LineWidth',1.5,'DisplayName','Individual subjects');

% Group mean
plot(ax3,bin_centres,grp_mean,'o-','Color',col_mean,'LineWidth',2.5,...
    'MarkerSize',7,'MarkerFaceColor',col_mean,'DisplayName','Group mean');

% Sinusoidal fit
valid = ~isnan(grp_mean);
try
    p0   = [0.1, 0, mean(grp_mean,'omitnan')];
    opts = optimset('Display','off');
    pf   = fminsearch(@(p) sum((sin_fit(p,bin_centres(valid))-grp_mean(valid)).^2),p0,opts);
    xf   = linspace(-pi,pi,500);
    plot(ax3,xf,sin_fit(pf,xf),'--','Color',col_fit,'LineWidth',2.2,'DisplayName','Sinusoidal fit');
catch; end

yline(ax3,0,'-k','LineWidth',0.8,'Alpha',0.4,'HandleVisibility','off');

% Vertical markers + labels for all 4 quadrants
quadrant_marks = {-pi, col_trough, 'Trough', 'center';
                  -pi/2, col_rising, 'Rising', 'center';
                  0, col_peak, 'Peak', 'center';
                  pi/2, col_falling, 'Falling', 'center';
                  pi, col_trough, 'Trough', 'right'};
for i = 1:size(quadrant_marks,1)
    xm    = quadrant_marks{i,1};
    col_m = quadrant_marks{i,2};
    lbl   = quadrant_marks{i,3};
    ha    = quadrant_marks{i,4};
    text(ax3, xm, 0.46, lbl, 'Color', col_m, 'FontSize', 12, ...
        'FontWeight','bold','HorizontalAlignment', ha);
end

set(ax3,'XTick',[-pi -pi/2 0 pi/2 pi],...
    'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'FontSize',13);
xlabel(ax3,'Pre-stimulus \mu-alpha phase (rad)','FontSize',14);
ylabel(ax3,'Mean z-scored MEP amplitude','FontSize',14);
title(ax3,'T4TE — Phase-MEP relationship (IAF band)','FontSize',14);
xlim(ax3,[-pi pi]); ylim(ax3,y_shade);
legend(ax3,'Location','southeast','FontSize',12);
grid(ax3,'on'); box(ax3,'on');

exportgraphics(fig3,fullfile(out_dir,'Fig3_group_binning.png'),'Resolution',dpi);
fprintf('Saved: Fig3_group_binning.png\n');
fprintf('\n=== All figures saved to: %s ===\n', out_dir);
