% raytracing/raytrace_two_ray.m
% フェーズ4②: 地面2波モデル（直接波＋地面反射波）による識別性能評価
%
% モデル: H = h_direct + h_reflect  (鏡像法, 反射係数 Γ=-1)
%   h_direct  = (λ/4πd1) * exp(-j*2π*fc*τ1)
%   h_reflect = Γ * (λ/4πd2) * exp(-j*2π*fc*τ2)
%
% 距離依存の位相差により SNR が自由空間から大きく乖離する（深いヌルは −∞ dB）。
% 自由空間との精度差を示すことで「どの伝搬現象が識別を壊すか」を切り分ける。

base_dir = fileparts(which('raytrace_two_ray'));
addpath(fullfile(base_dir, '..', 'signals'));
addpath(fullfile(base_dir, '..', 'features'));
addpath(fullfile(base_dir, '..', 'classifier'));

%% シナリオパラメータ
h_tx  = 5;    % [m] ETC路側機アンテナ高
h_rx  = 50;   % [m] UAV飛行高度
Gamma = -1;   % 反射係数（完全平面地, 近水平入射 → -1）

% ETC (PU): ARIB STD-T75
f_etc   = 5.80e9;  % Hz
ptx_etc = 10;      % dBm
bw_etc  = 4.4e6;   % Hz

% UAV (SU): 5.8GHz帯 空きチャネル（5.82GHz, DSRCサブバンド, プランB）
f_uav   = 5.82e9;  % Hz
ptx_uav = 20;      % dBm
bw_uav  = 20e6;    % Hz

nf_db   = 5;       % dB 受信機雑音指数（共通）

%% 電力パターン（細グリッド可視化: 干渉縞の確認）
dist_fine = 10:1:2000;
N_fine    = length(dist_fine);

snr_etc_fs_fine   = arrayfun(@(d) fs_snr_3d(d, h_tx, h_rx, f_etc, ptx_etc, nf_db, bw_etc),   dist_fine);
snr_etc_2ray_fine = arrayfun(@(d) tr_snr(d,   h_tx, h_rx, f_etc, ptx_etc, nf_db, bw_etc, Gamma), dist_fine);

%% 識別精度評価（自由空間と同じ8距離点）
dist_list = [10, 20, 50, 100, 200, 500, 1000, 2000];
N_dist    = length(dist_list);
N_trials  = 100;
N_methods = 3;

[~, ep_ref] = gen_ETC('seed', 1);
uw  = ep_ref.uw;
sps = ep_ref.sps;

thresh_fixed = [1.0, 0.5, 0.5];  % kurtosis / preamble / combined
acc_fs   = zeros(N_dist, N_methods);
acc_2ray = zeros(N_dist, N_methods);
f1_2ray  = zeros(N_dist, N_methods);

snr_etc_fs_c   = zeros(1, N_dist);
snr_etc_2ray_c = zeros(1, N_dist);

fprintf('2波モデル識別実験中 (%d距離 × %d試行/クラス) ...\n', N_dist, N_trials);
for di = 1:N_dist
    d = dist_list(di);
    snr_e_fs = fs_snr_3d(d, h_tx, h_rx, f_etc, ptx_etc, nf_db, bw_etc);
    snr_u_fs = fs_snr_3d(d, h_tx, h_rx, f_uav, ptx_uav, nf_db, bw_uav);
    snr_e_2r = tr_snr(d, h_tx, h_rx, f_etc, ptx_etc, nf_db, bw_etc, Gamma);
    snr_u_2r = tr_snr(d, h_tx, h_rx, f_uav, ptx_uav, nf_db, bw_uav, Gamma);

    snr_etc_fs_c(di)   = snr_e_fs;
    snr_etc_2ray_c(di) = snr_e_2r;

    scores_fs  = zeros(N_trials*2, N_methods);
    scores_2r  = zeros(N_trials*2, N_methods);
    labels_true = zeros(N_trials*2, 1);

    for ti = 1:N_trials
        seed = ti * 37;

        % ETC試行 (label=1)
        [etc, ~] = gen_ETC('seed', seed);
        k = calc_kurtosis(awgn_add(etc, snr_e_fs, true));
        p = calc_preamble(awgn_add(etc, snr_e_fs, true), uw, sps);
        scores_fs(ti,:) = make_scores(k, p);
        k = calc_kurtosis(awgn_add(etc, snr_e_2r, true));
        p = calc_preamble(awgn_add(etc, snr_e_2r, true), uw, sps);
        scores_2r(ti,:) = make_scores(k, p);
        labels_true(ti) = 1;

        % UAV試行 (label=0)
        [uav, ~] = gen_UAV('seed', seed);
        k = calc_kurtosis(awgn_add(uav, snr_u_fs, false));
        p = calc_preamble(awgn_add(uav, snr_u_fs, false), uw, sps);
        scores_fs(N_trials+ti,:) = make_scores(k, p);
        k = calc_kurtosis(awgn_add(uav, snr_u_2r, false));
        p = calc_preamble(awgn_add(uav, snr_u_2r, false), uw, sps);
        scores_2r(N_trials+ti,:) = make_scores(k, p);
        labels_true(N_trials+ti) = 0;
    end

    for mi = 1:N_methods
        [acc_fs(di,mi),   ~]         = calc_acc_f1(labels_true, scores_fs(:,mi)  >= thresh_fixed(mi));
        [acc_2ray(di,mi), f1_2ray(di,mi)] = calc_acc_f1(labels_true, scores_2r(:,mi) >= thresh_fixed(mi));
    end
    fprintf('  dist=%5dm  SNR_ETC:  自由空間=%5.1fdB  2波=%5.1fdB\n', d, snr_e_fs, snr_e_2r);
end

%% 結果表示
method_names = {'Kurtosis only', 'Preamble only', 'Combined    '};

fprintf('\n=== ETC受信SNR: 自由空間 vs 2波モデル ===\n');
fprintf('%-10s | %-14s %-14s %-14s\n', '距離[m]', 'SNR_fs[dB]', 'SNR_2ray[dB]', '差[dB]');
fprintf('%s\n', repmat('-',1,56));
for di = 1:N_dist
    diff_db = snr_etc_2ray_c(di) - snr_etc_fs_c(di);
    fprintf('%-10d | %-14.1f %-14.1f %-+14.1f\n', dist_list(di), snr_etc_fs_c(di), snr_etc_2ray_c(di), diff_db);
end

fprintf('\n=== 識別精度: 自由空間 vs 2波モデル (Preamble only) ===\n');
fprintf('%-10s | %-14s %-14s %-14s\n', '距離[m]', 'Acc_fs', 'Acc_2ray', '差');
fprintf('%s\n', repmat('-',1,56));
for di = 1:N_dist
    fprintf('%-10d | %-14.4f %-14.4f %-+14.4f\n', dist_list(di), ...
        acc_fs(di,2), acc_2ray(di,2), acc_2ray(di,2)-acc_fs(di,2));
end

fprintf('\n=== 識別精度 (2波モデル) vs 距離 ===\n');
fprintf('%-10s | %-14s %-14s %-14s\n', '距離[m]', method_names{:});
fprintf('%s\n', repmat('-',1,58));
for di = 1:N_dist
    fprintf('%-10d | %-14.4f %-14.4f %-14.4f\n', dist_list(di), acc_2ray(di,:));
end

fprintf('\n=== 精度差（2波 − 自由空間）===\n');
fprintf('%-10s | %-14s %-14s %-14s\n', '距離[m]', method_names{:});
fprintf('%s\n', repmat('-',1,58));
for di = 1:N_dist
    diff = acc_2ray(di,:) - acc_fs(di,:);
    fprintf('%-10d | %-+14.4f %-+14.4f %-+14.4f\n', dist_list(di), diff);
end

fprintf('\n=== 識別精度90%%を維持できる最大距離（2波モデル）===\n');
for mi = 1:N_methods
    idxs = find(acc_2ray(:,mi) >= 0.90);
    if isempty(idxs)
        fprintf('  %-16s: 最短10mでも90%%未達\n', method_names{mi});
    else
        fprintf('  %-16s: %5dm まで  (ETC 2波SNR ≈ %.1fdB)\n', ...
            method_names{mi}, dist_list(idxs(end)), snr_etc_2ray_c(idxs(end)));
    end
end

%% プロット①: SNR vs 距離（自由空間 vs 2波 — 干渉縞）
fig1 = figure('Visible','off');
hold on; grid on; box on;
plot(dist_fine, snr_etc_fs_fine,   'b-',  'LineWidth', 1.4, 'DisplayName', '自由空間 (Friis+3D)');
plot(dist_fine, snr_etc_2ray_fine, 'r-',  'LineWidth', 0.7, 'DisplayName', '2波モデル (干渉縞)');
yline(0,   'k--', '0 dB',   'LabelHorizontalAlignment','left','FontSize',8);
yline(-10, 'k:',  '-10 dB', 'LabelHorizontalAlignment','left','FontSize',8);
set(gca,'XScale','log');
xlabel('水平距離 [m]'); ylabel('ETC受信SNR [dB]');
title('ETC受信SNR vs 距離: 自由空間 vs 地面2波モデル (h_{tx}=5m, h_{rx}=50m)');
legend('Location','northeast'); xlim([10 2000]); ylim([-40 60]);
saveas(fig1, fullfile(base_dir, 'snr_vs_distance_two_ray.png'));
fprintf('\nSNR比較グラフ保存: raytracing/snr_vs_distance_two_ray.png\n');

%% プロット②: 識別精度 vs 距離（自由空間 vs 2波, 3手法）
lstyles = {'-','--','-.'};
colors_fs   = [0.2 0.6 0.9];  % 青系: 自由空間
colors_2ray = [0.9 0.3 0.1];  % 赤系: 2波

fig2 = figure('Visible','off');
hold on; grid on; box on;
for mi = 1:N_methods
    plot(dist_list, acc_fs(:,mi), lstyles{mi}, 'Color', colors_fs, ...
        'LineWidth',1.5,'Marker','o','MarkerSize',5, ...
        'DisplayName', sprintf('%s (自由空間)', strtrim(method_names{mi})));
    plot(dist_list, acc_2ray(:,mi), lstyles{mi}, 'Color', colors_2ray, ...
        'LineWidth',1.5,'Marker','s','MarkerSize',5, ...
        'DisplayName', sprintf('%s (2波)', strtrim(method_names{mi})));
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','left','FontSize',8);
set(gca,'XScale','log');
xlabel('水平距離 [m]'); ylabel('Accuracy');
title('識別精度 vs 距離: 自由空間 vs 地面2波モデル');
legend('Location','southwest','FontSize',7,'NumColumns',2);
ylim([0 1]);
saveas(fig2, fullfile(base_dir, 'accuracy_vs_distance_two_ray.png'));
fprintf('精度比較グラフ保存: raytracing/accuracy_vs_distance_two_ray.png\n');

disp('Phase 4-② raytrace_two_ray: OK');


%% ローカル関数
function snr_db = fs_snr_3d(d_m, h_tx, h_rx, f_hz, ptx_dbm, nf_db, bw_hz)
    % 直接波のみ（3D経路長使用）
    c  = 3e8;
    d1 = sqrt(d_m^2 + (h_tx - h_rx)^2);
    pl_db = 20*log10(4*pi*d1*f_hz/c);
    kB    = 1.38e-23; T0 = 290;
    nf_dbm = 10*log10(kB*T0*bw_hz) + 30 + nf_db;
    snr_db = ptx_dbm - pl_db - nf_dbm;
end

function snr_db = tr_snr(d_m, h_tx, h_rx, f_hz, ptx_dbm, nf_db, bw_hz, Gamma)
    % 2波モデル: H = h_direct + Gamma*h_reflect
    c      = 3e8;
    lambda = c / f_hz;
    d1 = sqrt(d_m^2 + (h_tx - h_rx)^2);   % 直接波経路長
    d2 = sqrt(d_m^2 + (h_tx + h_rx)^2);   % 反射波経路長（鏡像法）
    h1 = (lambda/(4*pi*d1)) * exp(-1j*2*pi*f_hz*d1/c);
    h2 = Gamma * (lambda/(4*pi*d2)) * exp(-1j*2*pi*f_hz*d2/c);
    H  = h1 + h2;
    % 有効経路損失（|H|→0の深いヌルは -50dB にクランプ）
    H_mag = max(abs(H), 10^(-50/20) * (lambda/(4*pi*d1)));
    pl_eff_db = -20*log10(H_mag);
    kB  = 1.38e-23; T0 = 290;
    nf_dbm = 10*log10(kB*T0*bw_hz) + 30 + nf_db;
    snr_db = ptx_dbm - pl_eff_db - nf_dbm;
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
