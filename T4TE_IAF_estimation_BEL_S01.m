% =========================================================================
% T4TE Study - IAF Estimation — BEL_S01
% =========================================================================
% Input:   BEL_S01_RS_clean.mat
% Output:  BEL_S01_IAF_results.mat + BEL_S01_IAF_PSD.png
%
% Method:  Centre-of-gravity (CoG) in 8-13 Hz
%          Hjorth-Laplacian on C3 (FC1, CP1, FC5, CP5)
% =========================================================================

clear; close all; clc

subj_code  = 'BEL_S01';
proc_path  = '/Users/e.w.m.dresens/Documents/master/Internship_Paolo/T4TE/data/BEL_S01/processed/';

hjorth_centre     = 'C3';
hjorth_neighbours = {'FC1', 'CP1', 'FC5', 'CP5'};
alpha_band        = [8 13];
fs                = 500;
win_len_s         = 2;
win_len           = round(win_len_s * fs);
n_overlap         = round(win_len / 2);
n_fft             = 2^nextpow2(win_len * 2);
max_cog_peak_diff = 1;

ft_path = '/Users/e.w.m.dresens/Documents/MATLAB/fieldtrip-20250106/';
addpath(ft_path); ft_defaults;

fprintf('\n=== STEP 1: Loading preprocessed data ===\n');
load_path = fullfile(proc_path, [subj_code '_RS_clean.mat']);
if ~exist(load_path, 'file')
    error('File not found: %s\nRun T4TE_RS_preprocessing_BEL_S01.m first.', load_path);
end
load(load_path, 'data_clean');
n_epochs = numel(data_clean.trial);
n_samp   = size(data_clean.trial{1}, 2);
fprintf('Loaded: %d epochs | %d channels\n', n_epochs, numel(data_clean.label));

fprintf('\n=== STEP 2: Hjorth-Laplacian on C3 ===\n');
missing = setdiff([{hjorth_centre}, hjorth_neighbours], data_clean.label);
if ~isempty(missing)
    error('Missing channels: %s', strjoin(missing, ', '));
end
idx_centre = find(strcmp(data_clean.label, hjorth_centre));
idx_neigh  = cellfun(@(c) find(strcmp(data_clean.label, c)), hjorth_neighbours);

hjorth_signal = zeros(1, n_epochs * n_samp);
ptr = 1;
for t = 1:n_epochs
    trial = data_clean.trial{t};
    hjorth_signal(ptr:ptr+n_samp-1) = trial(idx_centre,:) - mean(trial(idx_neigh,:), 1);
    ptr = ptr + n_samp;
end
fprintf('Total signal duration: %.1f s\n', length(hjorth_signal)/fs);

fprintf('\n=== STEP 3: Welch PSD ===\n');
[pxx, f] = pwelch(hjorth_signal, hanning(win_len), n_overlap, n_fft, fs);
freq_res = f(2) - f(1);
fprintf('Frequency resolution: %.4f Hz\n', freq_res);

fprintf('\n=== STEP 4: IAF estimation ===\n');
alpha_idx = f >= alpha_band(1) & f <= alpha_band(2);
f_alpha   = f(alpha_idx);
pxx_alpha = pxx(alpha_idx);
IAF_CoG   = sum(f_alpha .* pxx_alpha) / sum(pxx_alpha);
[peak_pwr, peak_idx] = max(pxx_alpha);
IAF_peak  = f_alpha(peak_idx);
delta_iaf = abs(IAF_CoG - IAF_peak);

fprintf('IAF (CoG):  %.3f Hz\n', IAF_CoG);
fprintf('IAF (peak): %.3f Hz\n', IAF_peak);
if delta_iaf > max_cog_peak_diff
    warning('CoG (%.2f Hz) and peak (%.2f Hz) differ %.2f Hz — inspect PSD.', IAF_CoG, IAF_peak, delta_iaf);
else
    fprintf('Quality check: CoG and peak within %.2f Hz — OK\n', delta_iaf);
end
if IAF_CoG < alpha_band(1) || IAF_CoG > alpha_band(2)
    warning('IAF CoG (%.2f Hz) outside [%d-%d Hz].', IAF_CoG, alpha_band(1), alpha_band(2));
end

fprintf('\n=== STEP 5: PSD figure ===\n');
fig = figure('Name', sprintf('PSD RS — %s', subj_code), ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [100 100 1000 460]);

subplot(1,2,1);
semilogy(f, pxx, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2); hold on;
yl = ylim;
fill([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.85 0.92 1.0], 'EdgeColor','none','FaceAlpha',0.45);
xline(IAF_CoG,  'b--', sprintf('CoG: %.2f Hz', IAF_CoG),  'LineWidth', 1.5);
xline(IAF_peak, 'r:',  sprintf('Peak: %.2f Hz', IAF_peak), 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Power density (uV2/Hz)');
title(sprintf('%s — full spectrum (1-45 Hz)', subj_code));
xlim([1 45]); grid on;
legend({'PSD','Alpha band','IAF CoG','IAF peak'}, 'Location','northeast');

subplot(1,2,2);
plot(f_alpha, pxx_alpha, 'Color', [0.2 0.4 0.8], 'LineWidth', 2); hold on;
fill([f_alpha; flipud(f_alpha)], [pxx_alpha; zeros(size(pxx_alpha))], ...
    [0.7 0.82 0.96], 'EdgeColor','none','FaceAlpha',0.4);
xline(IAF_CoG,  'b--', sprintf('CoG: %.2f Hz', IAF_CoG),  'LineWidth', 1.5);
xline(IAF_peak, 'r:',  sprintf('Peak: %.2f Hz', IAF_peak), 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Power density (uV2/Hz)');
title(sprintf('Alpha band (%d-%d Hz)', alpha_band(1), alpha_band(2)));
xlim(alpha_band); grid on;
sgtitle(sprintf('%s — Resting State IAF (Hjorth C3)', subj_code));

exportgraphics(fig, fullfile(proc_path, [subj_code '_IAF_PSD.png']), 'Resolution', 200);
fprintf('Figure saved.\n');

fprintf('\n=== STEP 6: Saving ===\n');
IAF_results.subject           = subj_code;
IAF_results.IAF_CoG           = IAF_CoG;
IAF_results.IAF_peak          = IAF_peak;
IAF_results.peak_power        = peak_pwr;
IAF_results.cog_peak_diff_hz  = delta_iaf;
IAF_results.f                 = f;
IAF_results.pxx               = pxx;
IAF_results.alpha_band        = alpha_band;
IAF_results.hjorth_centre     = hjorth_centre;
IAF_results.hjorth_neighbours = hjorth_neighbours;
IAF_results.win_len_s         = win_len_s;
IAF_results.n_fft             = n_fft;
IAF_results.freq_res_hz       = freq_res;
IAF_results.n_epochs          = n_epochs;
IAF_results.total_duration_s  = length(hjorth_signal)/fs;
IAF_results.date              = datestr(now, 'yyyy-mm-dd HH:MM');

save(fullfile(proc_path, [subj_code '_IAF_results.mat']), 'IAF_results', '-v7.3');
fprintf('IAF complete — %s: CoG = %.3f Hz | Peak = %.3f Hz\n', subj_code, IAF_CoG, IAF_peak);
