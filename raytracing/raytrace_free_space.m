% raytracing/raytrace_free_space.m
% フェーズ4①: 自由空間伝搬モデルによる識別性能評価
%
% Friis伝送公式で距離→受信SNRを計算し、3手法（尖度/プリアンブル/結合）の
% 識別精度を距離の関数として評価する。
% 反射・回折ゼロ、直接波のみ（フェーズ4の基準ライン）。

base_dir = fileparts(which('raytrace_free_space'));
addpath(fullfile(base_dir, '..', 'signals'));
addpath(fullfile(base_dir, '..', 'features'));
addpath(fullfile(base_dir, '..', 'classifier'));

%% 伝搬シナリオパラメータ
% ──────────────────────────────────────────
% ETC (PU)  : ARIB STD-T75 路側機
%   5.8GHz, Ptx=10dBm, 占有BW=4.4MHz, アンテナ利得 ≈ 0dBi（等方性近似）
% UAV (SU)  : 802.11a 5GHz帯（ch36 = 5.18GHz）
%   Ptx=20dBm, BW=20MHz, アンテナ利得 ≈ 0dBi
% 受信機雑音指数: 5dB（共通）
% 参照温度: 290K（IEEE標準）
% ──────────────────────────────────────────
f_etc   = 5.80e9;   % [Hz]  ETC中心周波数
ptx_etc = 10;       % [dBm] ETC送信電力
bw_etc  = 4.4e6;    % [Hz]  ETC受信バンド幅

f_uav   = 5.18e9;   % [Hz]  UAV中心周波数（802.11a ch36）
ptx_uav = 20;       % [dBm] UAV送信電力
bw_uav  = 20e6;     % [Hz]  UAV受信バンド幅

nf_db   = 5;        % [dB]  受信機雑音指数（ETC/UAV共通）

% 評価距離 [m]: 路側機〜UAV間の直線距離
dist_list = [10, 20, 50, 100, 200, 500, 1000, 2000];
N_dist    = length(dist_list);
N_trials  = 100;    % クラスあたり試行数
N_methods = 3;

%% テンプレートETC生成（UW・sps取得用）
[~, ep_ref] = gen_ETC('seed', 1);
uw  = ep_ref.uw;
sps = ep_ref.sps;

%% 距離 → 受信SNR変換（Friis公式）
snr_etc_vec = arrayfun(@(d) free_space_snr(d, f_etc, ptx_etc, nf_db, bw_etc), dist_list);
snr_uav_vec = arrayfun(@(d) free_space_snr(d, f_uav, ptx_uav, nf_db, bw_uav), dist_list);

fprintf('=== 自由空間 距離 → 受信SNR（Friis公式）===\n');
fprintf('%-10s | %-12s %-12s\n', '距離[m]', 'ETC SNR[dB]', 'UAV SNR[dB]');
fprintf('%s\n', repmat('-', 1, 38));
for di = 1:N_dist
    fprintf('%-10d | %-12.1f %-12.1f\n', dist_list(di), snr_etc_vec(di), snr_uav_vec(di));
end

%% 識別実験: 距離 × 3手法
% ETC試行 → ETC TX SNRで劣化, UAV試行 → UAV TX SNRで劣化（送信機依存）
thresh_fixed = [1.0, 0.5, 0.5];  % kurtosis / preamble / combined
acc_mat = zeros(N_dist, N_methods);
f1_mat  = zeros(N_dist, N_methods);

fprintf('\n識別実験中 (%d距離 × %d試行/クラス) ...\n', N_dist, N_trials);
for di = 1:N_dist
    snr_e = snr_etc_vec(di);
    snr_u = snr_uav_vec(di);

    scores_d    = zeros(N_trials*2, N_methods);
    labels_true = zeros(N_trials*2, 1);

    for ti = 1:N_trials
        seed = ti * 37;

        % ETC試行 (true label = 1): ETC送信パワー・周波数でのSNR
        [etc, ~] = gen_ETC('seed', seed);
        etc_n = awgn_add(etc, snr_e, true);
        k = calc_kurtosis(etc_n);
        p = calc_preamble(etc_n, uw, sps);
        scores_d(ti, :)     = make_scores(k, p);
        labels_true(ti)     = 1;

        % UAV試行 (true label = 0): UAV送信パワー・周波数でのSNR
        [uav, ~] = gen_UAV('seed', seed);
        uav_n = awgn_add(uav, snr_u, false);
        k = calc_kurtosis(uav_n);
        p = calc_preamble(uav_n, uw, sps);
        scores_d(N_trials+ti, :) = make_scores(k, p);
        labels_true(N_trials+ti) = 0;
    end

    for mi = 1:N_methods
        y_pred = scores_d(:, mi) >= thresh_fixed(mi);
        [acc_mat(di,mi), f1_mat(di,mi)] = calc_acc_f1(labels_true, y_pred);
    end
    fprintf('  dist=%5dm  (SNR_ETC=%5.1fdB, SNR_UAV=%5.1fdB) 完了\n', ...
        dist_list(di), snr_e, snr_u);
end

%% 結果表示
method_names = {'Kurtosis only', 'Preamble only', 'Combined    '};

fprintf('\n=== 識別精度 (Accuracy) vs 距離 ===\n');
fprintf('%-10s | %-14s %-14s %-14s\n', '距離[m]', method_names{:});
fprintf('%s\n', repmat('-', 1, 58));
for di = 1:N_dist
    fprintf('%-10d | %-14.4f %-14.4f %-14.4f\n', dist_list(di), acc_mat(di,:));
end

fprintf('\n=== F1スコア vs 距離 ===\n');
fprintf('%-10s | %-14s %-14s %-14s\n', '距離[m]', method_names{:});
fprintf('%s\n', repmat('-', 1, 58));
for di = 1:N_dist
    fprintf('%-10d | %-14.4f %-14.4f %-14.4f\n', dist_list(di), f1_mat(di,:));
end

% 90%精度を維持できる最大距離
fprintf('\n=== 識別精度90%%を維持できる最大距離 ===\n');
for mi = 1:N_methods
    idxs = find(acc_mat(:,mi) >= 0.90);
    if isempty(idxs)
        fprintf('  %-16s: 最短距離(10m)でも90%%未達\n', method_names{mi});
    else
        fprintf('  %-16s: %5d m まで (ETC SNR ≈ %.1f dB)\n', ...
            method_names{mi}, dist_list(idxs(end)), snr_etc_vec(idxs(end)));
    end
end

%% プロット①: 識別精度 vs 距離（対数軸）
lstyles = {'-', '--', '-.'};
colors  = lines(N_methods);

fig1 = figure('Visible', 'off');
hold on; grid on; box on;
for mi = 1:N_methods
    plot(dist_list, acc_mat(:,mi), lstyles{mi}, 'Color', colors(mi,:), ...
        'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5, ...
        'DisplayName', strtrim(method_names{mi}));
end
yline(0.9, 'k--', '90%', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
set(gca, 'XScale', 'log');
xlabel('ETC〜UAV 距離 [m]'); ylabel('Accuracy');
title('識別精度 vs 距離（自由空間）');
legend('Location', 'southwest'); ylim([0 1]);
saveas(fig1, fullfile(base_dir, 'accuracy_vs_distance_freespace.png'));
fprintf('\n精度グラフ保存: raytracing/accuracy_vs_distance_freespace.png\n');

%% プロット②: 識別精度 vs ETC受信SNR（距離から計算した物理SNR）
fig2 = figure('Visible', 'off');
hold on; grid on; box on;
for mi = 1:N_methods
    plot(snr_etc_vec, acc_mat(:,mi), lstyles{mi}, 'Color', colors(mi,:), ...
        'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5, ...
        'DisplayName', strtrim(method_names{mi}));
end
yline(0.9, 'k--', '90%', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
xlabel('ETC受信SNR [dB]（自由空間）'); ylabel('Accuracy');
title('識別精度 vs ETC受信SNR（自由空間）');
legend('Location', 'southeast'); ylim([0 1]);
saveas(fig2, fullfile(base_dir, 'accuracy_vs_snr_freespace.png'));
fprintf('SNRグラフ保存: raytracing/accuracy_vs_snr_freespace.png\n');

disp('Phase 4-① raytrace_free_space: OK');


%% ローカル関数
function snr_db = free_space_snr(d_m, f_hz, ptx_dbm, nf_db, bw_hz)
    % Friis伝送公式: SNR = Ptx - PL - N_floor
    % PL = 20*log10(4*pi*d*f/c) [dB], N_floor = 10*log10(kTB) + 30 + NF [dBm]
    c  = 3e8;
    pl_db = 20*log10(4*pi*d_m*f_hz/c);
    kB    = 1.38e-23;
    T0    = 290;
    noise_floor_dbm = 10*log10(kB*T0*bw_hz) + 30 + nf_db;
    snr_db = ptx_dbm - pl_db - noise_floor_dbm;
end

function sc = make_scores(kurt, p_peak)
    k_norm = max(0, min(1, (3 - kurt) / 2));
    sc = [3 - kurt, p_peak, 0.5*k_norm + 0.5*p_peak];
end

function y = awgn_add(x, snr_db, is_real)
    pwr  = mean(abs(x).^2);
    npwr = pwr / 10^(snr_db/10);
    if is_real
        y = x + sqrt(npwr) * randn(size(x));
    else
        y = x + sqrt(npwr/2) * (randn(size(x)) + 1j*randn(size(x)));
    end
end

function [acc, f1] = calc_acc_f1(y_true, y_pred)
    tp  = sum(y_true==1 & y_pred==1);
    tn  = sum(y_true==0 & y_pred==0);
    fp  = sum(y_true==0 & y_pred==1);
    fn  = sum(y_true==1 & y_pred==0);
    acc = (tp+tn) / (tp+tn+fp+fn);
    pre = tp / max(tp+fp, 1);
    rec = tp / max(tp+fn, 1);
    f1  = 2*pre*rec / max(pre+rec, eps);
end
