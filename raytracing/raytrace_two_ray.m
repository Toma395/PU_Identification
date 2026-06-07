% raytracing/raytrace_two_ray.m
% フェーズ4②: 地面2波モデル（直接波＋地面反射波）による識別性能評価
%
% モデル: H = h_direct + h_reflect  (鏡像法, 反射係数 Γ=-1)
%
% 改訂: モンテカルロ評価（N_mc=800 ランダム距離 × N_tmc=25 試行/点）
%       散布+移動平均プロット + ヒストグラム追加

base_dir = fileparts(which('raytrace_two_ray'));
addpath(fullfile(base_dir, '..', 'signals'));
addpath(fullfile(base_dir, '..', 'features'));
addpath(fullfile(base_dir, '..', 'classifier'));
addpath(fullfile(base_dir, '..'));   % apply_paper_style

%% シナリオパラメータ
h_tx  = 5;    % [m] ETC路側機アンテナ高
h_rx  = 50;   % [m] UAV飛行高度
Gamma = -1;   % 反射係数（完全平面地）

f_etc   = 5.80e9;  ptx_etc = 10;  bw_etc = 4.4e6;
f_uav   = 5.82e9;  ptx_uav = 20;  bw_uav = 20e6;
nf_db   = 5;

%% SNR細グリッド可視化（干渉縞確認）
dist_fine         = 10:1:2000;
snr_etc_fs_fine   = arrayfun(@(d) fs_snr_3d(d,h_tx,h_rx,f_etc,ptx_etc,nf_db,bw_etc), dist_fine);
snr_etc_2ray_fine = arrayfun(@(d) tr_snr(d,  h_tx,h_rx,f_etc,ptx_etc,nf_db,bw_etc,Gamma), dist_fine);

%% テンプレート
[~, ep_ref] = gen_ETC('seed', 1);
uw  = ep_ref.uw;
sps = ep_ref.sps;
thresh_fixed = [1.0, 0.5, 0.5];
method_names = {'Kurtosis','Preamble','Combined'};
N_methods    = 3;

%% ─── モンテカルロ評価 ────────────────────────────────────────────
N_mc   = 800;    % ランダム距離サンプル数
N_tmc  = 25;     % 信号試行数/距離点/クラス
rng(42);
d_mc = round(10 .^ (rand(1,N_mc) * log10(2000/10) + log10(10)));   % 対数一様 [10,2000]m

fprintf('モンテカルロ評価: %d距離 × %d試行/クラス ...\n', N_mc, N_tmc);

snr_etc_mc = arrayfun(@(d) tr_snr(d,h_tx,h_rx,f_etc,ptx_etc,nf_db,bw_etc,Gamma), d_mc);
snr_uav_mc = arrayfun(@(d) tr_snr(d,h_tx,h_rx,f_uav,ptx_uav,nf_db,bw_uav,Gamma), d_mc);

acc_mc = zeros(N_mc, N_methods);

for di = 1:N_mc
    snr_e = snr_etc_mc(di);
    snr_u = snr_uav_mc(di);

    scores_mc  = zeros(N_tmc*2, N_methods);
    labels_mc  = [ones(N_tmc,1); zeros(N_tmc,1)];

    for ti = 1:N_tmc
        seed_etc = ti * 37  + di * 10000;
        seed_uav = ti * 41  + di * 10000 + 9999983;

        [etc, ~] = gen_ETC('seed', seed_etc);
        k = calc_kurtosis(awgn_add(etc, snr_e, true));
        p = calc_preamble(awgn_add(etc, snr_e, true), uw, sps);
        scores_mc(ti,:) = make_scores(k, p);

        [uav, ~] = gen_UAV('seed', seed_uav);
        k = calc_kurtosis(awgn_add(uav, snr_u, false));
        p = calc_preamble(awgn_add(uav, snr_u, false), uw, sps);
        scores_mc(N_tmc+ti,:) = make_scores(k, p);
    end

    for mi = 1:N_methods
        [acc_mc(di,mi),~] = calc_acc_f1(labels_mc, scores_mc(:,mi) >= thresh_fixed(mi));
    end

    if mod(di, 200)==0
        fprintf('  %d / %d 完了\n', di, N_mc);
    end
end
fprintf('モンテカルロ完了\n');

%% サマリ
fprintf('\n=== 2波モデル MC評価サマリ ===\n');
fprintf('%-12s | %-10s %-10s %-10s\n','','Kurtosis','Preamble','Combined');
fprintf('%s\n', repmat('-',1,44));
fprintf('%-12s | %-10.3f %-10.3f %-10.3f\n','平均精度',    mean(acc_mc));
fprintf('%-12s | %-10.3f %-10.3f %-10.3f\n','中央値精度',  median(acc_mc));
fprintf('%-12s | %-10.3f %-10.3f %-10.3f\n','90%%ile点数', mean(acc_mc >= 0.9));

% 移動平均用ソート
[d_sort, sort_idx] = sort(d_mc);
ma_window = 80;   % 移動平均ウィンドウ幅
acc_ma = zeros(N_mc, N_methods);
for mi = 1:N_methods
    acc_ma(:,mi) = movmean(acc_mc(sort_idx,mi), ma_window);
end

%% ─── プロット ────────────────────────────────────────────────────
lstyles = {'-','--','-.'};
colors  = lines(N_methods);

%% ①: SNR vs 距離（干渉縞）
fig1 = figure('Visible','off','Color','w','Position',[0 0 900 500]);
hold on; grid on; box on;
plot(dist_fine, snr_etc_fs_fine,   '-',  'Color',[0.2 0.6 0.9],'LineWidth',1.8,...
    'DisplayName','自由空間 (Friis+3D)');
plot(dist_fine, snr_etc_2ray_fine, '-',  'Color',[0.9 0.3 0.1],'LineWidth',0.9,...
    'DisplayName','2波モデル (干渉縞)');
yline(0,   'k--','0 dB',   'LabelHorizontalAlignment','left','FontSize',11,'HandleVisibility','off');
yline(-10, 'k:', '-10 dB', 'LabelHorizontalAlignment','left','FontSize',11,'HandleVisibility','off');
set(gca,'XScale','log');
xlabel('水平距離 [m]','FontSize',14);
ylabel('ETC受信SNR [dB]','FontSize',14);
title(sprintf('ETC受信SNR vs 距離  h_{tx}=%.0fm, h_{rx}=%.0fm', h_tx, h_rx),'FontSize',16);
legend('Location','northeast','FontSize',12);
xlim([10 2000]); ylim([-50 60]);
apply_paper_style(fig1);
saveas(fig1, fullfile(base_dir,'snr_vs_distance_two_ray.png'));
fprintf('\nSNR比較保存: raytracing/snr_vs_distance_two_ray.png\n');

%% ②: MC散布 + 移動平均（3手法オーバーレイ）
fig2 = figure('Visible','off','Color','w','Position',[0 0 900 520]);
hold on; grid on; box on;
% 散布（薄く）
for mi = 1:N_methods
    scatter(d_mc, acc_mc(:,mi), 10, colors(mi,:), 'o', ...
        'MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.15,'HandleVisibility','off');
end
% 移動平均（線）
for mi = 1:N_methods
    plot(d_sort, acc_ma(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.5,'DisplayName',method_names{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
set(gca,'XScale','log');
xlabel('水平距離 [m]（対数軸）','FontSize',14);
ylabel('Accuracy','FontSize',14);
title(sprintf('識別精度 vs 距離  2波モデル (N_{mc}=%d, MA窓=%d点)', N_mc, ma_window),'FontSize',16);
legend('Location','southwest','FontSize',12);
ylim([0.3 1.05]); xlim([10 2000]);
apply_paper_style(fig2);
saveas(fig2, fullfile(base_dir,'accuracy_vs_distance_two_ray.png'));
fprintf('散布+移動平均保存: raytracing/accuracy_vs_distance_two_ray.png\n');

%% ③: 精度ヒストグラム（3手法）
fig3 = figure('Visible','off','Color','w','Position',[0 0 1100 420]);
for mi = 1:N_methods
    subplot(1,3,mi);
    histogram(acc_mc(:,mi), 21, 'FaceColor',colors(mi,:), ...
        'FaceAlpha',0.75,'EdgeColor','w');
    hold on; grid on; box on;
    xline(0.9,'k--','LineWidth',1.5,'HandleVisibility','off');
    xlabel('Accuracy','FontSize',13);
    ylabel('距離点数','FontSize',13);
    title(method_names{mi},'FontSize',14);
    xlim([0.3 1.05]);
    pct_above = 100 * mean(acc_mc(:,mi) >= 0.9);
    text(0.36, max(ylim)*0.85, sprintf('≥90%%: %.0f%%点', pct_above), ...
        'FontSize',11,'Color','k');
end
sgtitle(sprintf('識別精度分布  2波モデル (N_{mc}=%d距離点)', N_mc), ...
    'FontSize',14,'Color',[0.15 0.15 0.15]);
apply_paper_style(fig3);
saveas(fig3, fullfile(base_dir,'accuracy_mc_hist.png'));
fprintf('ヒストグラム保存: raytracing/accuracy_mc_hist.png\n');

disp('Phase 4-② raytrace_two_ray: OK');


%% ローカル関数
function snr_db = fs_snr_3d(d_m, h_tx, h_rx, f_hz, ptx_dbm, nf_db, bw_hz)
    c  = 3e8;
    d1 = sqrt(d_m^2 + (h_tx - h_rx)^2);
    pl_db = 20*log10(4*pi*d1*f_hz/c);
    kB    = 1.38e-23; T0 = 290;
    nf_dbm = 10*log10(kB*T0*bw_hz) + 30 + nf_db;
    snr_db = ptx_dbm - pl_db - nf_dbm;
end

function snr_db = tr_snr(d_m, h_tx, h_rx, f_hz, ptx_dbm, nf_db, bw_hz, Gamma)
    c      = 3e8;
    lambda = c / f_hz;
    d1 = sqrt(d_m^2 + (h_tx - h_rx)^2);
    d2 = sqrt(d_m^2 + (h_tx + h_rx)^2);
    h1 = (lambda/(4*pi*d1)) * exp(-1j*2*pi*f_hz*d1/c);
    h2 = Gamma * (lambda/(4*pi*d2)) * exp(-1j*2*pi*f_hz*d2/c);
    H  = h1 + h2;
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
