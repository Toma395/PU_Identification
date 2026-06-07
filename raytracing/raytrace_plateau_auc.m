% raytracing/raytrace_plateau_auc.m
% raytrace_plateau.m のコピー＋生スコア保存追加版
% 変更点: 格子点ごとの連続スコアを scores_{k,p,c}_map に保存し、
%         NLOS セル集約 AUC（perfcurve）と ROC 曲線を追加。
% 出力: plateau_results.mat（スコア配列追加）
%       experiments/plateau_nlos_auc.png

base_dir = fileparts(which('raytrace_plateau_auc'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..','citymodel'));
addpath(fullfile(base_dir,'..'));          % apply_paper_style

%% ─────────────────────────────────────────────────────────────────
%% パラメータ（raytrace_plateau.m と同一）
%% ─────────────────────────────────────────────────────────────────
etc_xy       = [-250, -300];
etc_z        =  5;
f_etc        = 5.80e9;
ptx_etc      = 10;
bw_etc       = 4.4e6;

f_uav        = 5.82e9;
ptx_uav      = 20;
bw_uav       = 20e6;
uav_z        = 100;

nf_db        = 5;
wall_loss_db = 12;

N_grid       = 50;
grid_range   = 500;

N_trials     = 30;

fprintf('=== raytrace_plateau_auc (Phase4-③ Step1 + AUC) ===\n');
fprintf('ETC: (%.0f, %.0f) m, z=%.0fm, f=%.2fGHz, Ptx=%.0fdBm\n', ...
    etc_xy(1),etc_xy(2), etc_z, f_etc/1e9, ptx_etc);
fprintf('UAV: z=%.0fm, 格子%d×%d, ±%.0fm, N_trials=%d/クラス\n', ...
    uav_z, N_grid, N_grid, grid_range, N_trials);

%% ─────────────────────────────────────────────────────────────────
%% 建物データ読込 & AABB
%% ─────────────────────────────────────────────────────────────────
city_dir  = fullfile(base_dir,'..','citymodel','data');
gml_files = dir(fullfile(city_dir,'*.gml'));
if isempty(gml_files); error('citymodel/data/ にGMLファイルが見つかりません'); end

fprintf('建物データ読込中...\n');
[buildings, info] = load_plateau(fullfile(city_dir, gml_files(1).name));
N_bldg = info.n_buildings;
fprintf('  %d棟 読込完了\n', N_bldg);

aabb = zeros(N_bldg, 5);
for i = 1:N_bldg
    fp = buildings(i).footprint;
    h  = buildings(i).height;
    if isnan(h); h = 0; end
    aabb(i,:) = [min(fp(:,1)), max(fp(:,1)), min(fp(:,2)), max(fp(:,2)), h];
end
fprintf('  AABB 完了\n\n');

%% ─────────────────────────────────────────────────────────────────
%% 格子 & テンプレート
%% ─────────────────────────────────────────────────────────────────
x_vec = linspace(etc_xy(1)-grid_range, etc_xy(1)+grid_range, N_grid);
y_vec = linspace(etc_xy(2)-grid_range, etc_xy(2)+grid_range, N_grid);
[X_grid, Y_grid] = meshgrid(x_vec, y_vec);

[~, ep] = gen_ETC('seed', 1);
uw  = ep.uw;
sps = ep.sps;
thresh = [1.0, 0.5, 0.5];

%% ─────────────────────────────────────────────────────────────────
%% メインループ
%% ─────────────────────────────────────────────────────────────────
n_walls_map  = zeros(N_grid, N_grid);
snr_etc_map  = zeros(N_grid, N_grid);
acc_map      = zeros(N_grid, N_grid, 3);

N_per_pt     = N_trials * 2;   % ETC×N_trials + UAV×N_trials
scores_k_map = zeros(N_grid, N_grid, N_per_pt);  % 尖度スコア (3-kurt)
scores_p_map = zeros(N_grid, N_grid, N_per_pt);  % プリアンブル相関ピーク
scores_c_map = zeros(N_grid, N_grid, N_per_pt);  % 結合スコア
labels_map   = zeros(N_grid, N_grid, N_per_pt);  % ラベル (ETC=1, UAV=0)

etc_pos3 = [etc_xy, etc_z];

fprintf('遮蔽判定 + 識別実験開始 (%d×%d格子, %d試行/クラス)...\n', ...
    N_grid, N_grid, N_trials);
t0 = tic;

for gi = 1:N_grid
    for gj = 1:N_grid
        uav_pos3 = [X_grid(gi,gj), Y_grid(gi,gj), uav_z];

        nw = count_walls(etc_pos3, uav_pos3, buildings, aabb);
        n_walls_map(gi,gj) = nw;

        d3d   = norm(uav_pos3 - etc_pos3);
        snr_e = fs_snr(d3d, f_etc, ptx_etc, nf_db, bw_etc) - wall_loss_db*nw;
        snr_u = fs_snr(d3d, f_uav, ptx_uav, nf_db, bw_uav) - wall_loss_db*nw;
        snr_etc_map(gi,gj) = snr_e;

        scores = zeros(N_per_pt, 3);
        labels = [ones(N_trials,1); zeros(N_trials,1)];
        for ti = 1:N_trials
            seed = ti*41 + gi*7 + gj*3;
            [es,~] = gen_ETC('seed', seed);
            k = calc_kurtosis(awgn_add(es, snr_e, true));
            p = calc_preamble( awgn_add(es, snr_e, true), uw, sps);
            scores(ti,:) = mk_score(k,p);

            [us,~] = gen_UAV('seed', seed);
            k = calc_kurtosis(awgn_add(us, snr_u, false));
            p = calc_preamble( awgn_add(us, snr_u, false), uw, sps);
            scores(N_trials+ti,:) = mk_score(k,p);
        end

        for mi = 1:3
            acc_map(gi,gj,mi) = mean((scores(:,mi) >= thresh(mi)) == labels);
        end

        % 生スコア保存（閾値判定前の連続値）
        scores_k_map(gi,gj,:) = scores(:,1);
        scores_p_map(gi,gj,:) = scores(:,2);
        scores_c_map(gi,gj,:) = scores(:,3);
        labels_map(gi,gj,:)   = labels;
    end

    if mod(gi,5)==0
        el = toc(t0);
        fprintf('  行 %2d/%d 完了  経過%.0fs  推定残%.0fs\n', ...
            gi, N_grid, el, el/gi*(N_grid-gi));
    end
end
fprintf('全完了  総計 %.1fs\n\n', toc(t0));

%% ─────────────────────────────────────────────────────────────────
%% 統計サマリ
%% ─────────────────────────────────────────────────────────────────
los_mask  = (n_walls_map == 0);
nlos_mask = (n_walls_map  > 0);

fprintf('=== LOS/NLOS 分布 ===\n');
fprintf('  LOS  (0棟): %4d点 (%.1f%%)\n', sum(los_mask(:)),  100*mean(los_mask(:)));
fprintf('  NLOS (1棟): %4d点 (%.1f%%)\n', sum(n_walls_map(:)==1), 100*mean(n_walls_map(:)==1));
fprintf('  NLOS (2棟): %4d点 (%.1f%%)\n', sum(n_walls_map(:)==2), 100*mean(n_walls_map(:)==2));
fprintf('  NLOS (3+棟):%4d点 (%.1f%%)\n', sum(n_walls_map(:)>=3), 100*mean(n_walls_map(:)>=3));

fprintf('\n=== 識別精度 (Accuracy) ===\n');
fprintf('%-12s | 全体    LOS     NLOS\n','手法');
fprintf('%s\n',repmat('-',1,40));
mnames = {'Kurtosis','Preamble','Combined'};
for mi = 1:3
    a = acc_map(:,:,mi);
    fprintf('%-12s | %.3f   %.3f   %.3f\n', mnames{mi}, ...
        mean(a(:)), mean(a(los_mask)), mean(a(nlos_mask)));
end

snr_los  = snr_etc_map(los_mask);
snr_nlos = snr_etc_map(nlos_mask);
fprintf('\n=== ETC受信SNR ===\n');
fprintf('  LOS  : 中央値=%.1fdB  10%%ile=%.1fdB\n', median(snr_los),  prctile(snr_los,10));
fprintf('  NLOS : 中央値=%.1fdB  10%%ile=%.1fdB\n', median(snr_nlos), prctile(snr_nlos,10));

%% ─────────────────────────────────────────────────────────────────
%% NLOS セル集約 AUC（perfcurve）
%% ─────────────────────────────────────────────────────────────────
fprintf('\nNLOS AUC 計算中 (%d格子点 × %d試行/点)...\n', sum(nlos_mask(:)), N_per_pt);

nlos_idx = find(nlos_mask(:));
N_nlos   = numel(nlos_idx);

all_labels   = zeros(N_nlos * N_per_pt, 1);
all_scores_k = zeros(N_nlos * N_per_pt, 1);
all_scores_p = zeros(N_nlos * N_per_pt, 1);
all_scores_c = zeros(N_nlos * N_per_pt, 1);

for ni = 1:N_nlos
    [gi, gj] = ind2sub([N_grid, N_grid], nlos_idx(ni));
    row_s = (ni-1)*N_per_pt + 1;
    row_e =  ni   *N_per_pt;
    all_labels(row_s:row_e)   = squeeze(labels_map(gi,gj,:));
    all_scores_k(row_s:row_e) = squeeze(scores_k_map(gi,gj,:));
    all_scores_p(row_s:row_e) = squeeze(scores_p_map(gi,gj,:));
    all_scores_c(row_s:row_e) = squeeze(scores_c_map(gi,gj,:));
end

[~,~,~, auc_k] = perfcurve(all_labels, all_scores_k, 1);
[~,~,~, auc_p] = perfcurve(all_labels, all_scores_p, 1);
[~,~,~, auc_c] = perfcurve(all_labels, all_scores_c, 1);

fprintf('\nPLATEAU NLOS AUC:\n');
fprintf('  Kurtosis  AUC = %.3f\n', auc_k);
fprintf('  Preamble  AUC = %.3f\n', auc_p);
fprintf('  Combined  AUC = %.3f\n', auc_c);
diff_pc = auc_p - auc_c;
fprintf('  Combined AUC < Preamble AUC の差 = %.3f\n', diff_pc);
if diff_pc > 0.005
    fprintf('  → 本物の情報破壊あり（固定重みが識別能力を損なっている）\n');
else
    fprintf('  → AUC差は誤差範囲内（0.005未満）\n');
end

%% ─────────────────────────────────────────────────────────────────
%% プロット: NLOS ROC 曲線
%% ─────────────────────────────────────────────────────────────────
score_arrays = {all_scores_k, all_scores_p, all_scores_c};
lstyles = {'-','--','-.'};
colors  = lines(3);

fig_roc = figure('Visible','off','Color','w','Position',[0 0 900 480]);
hold on; grid on; box on;
for mi = 1:3
    [fpr_r, tpr_r, ~, auc_r] = perfcurve(all_labels, score_arrays{mi}, 1);
    plot(fpr_r, tpr_r, lstyles{mi}, 'Color',colors(mi,:), 'LineWidth',2.0, ...
        'DisplayName', sprintf('%s (AUC=%.3f)', mnames{mi}, auc_r));
end
plot([0 1],[0 1],'k:','LineWidth',1.0,'HandleVisibility','off');
xlabel('FPR','FontSize',14);
ylabel('TPR','FontSize',14);
title(sprintf('ROC @ PLATEAU NLOS格子点  (%d点×%d試行/点)', N_nlos, N_per_pt), 'FontSize',15);
legend('Location','southeast','FontSize',12);
xlim([0 1]); ylim([0 1]);
apply_paper_style(fig_roc);

exp_dir = fullfile(base_dir, '..', 'experiments');
saveas(fig_roc, fullfile(exp_dir, 'plateau_nlos_auc.png'));
fprintf('\n保存: experiments/plateau_nlos_auc.png\n');

%% ─────────────────────────────────────────────────────────────────
%% .mat 保存（生スコア追加）
%% ─────────────────────────────────────────────────────────────────
mat_path = fullfile(base_dir, 'plateau_results.mat');
save(mat_path, ...
    'acc_map', 'n_walls_map', 'snr_etc_map', ...
    'los_mask', 'nlos_mask', ...
    'x_vec', 'y_vec', 'X_grid', 'Y_grid', ...
    'etc_xy', 'uav_z', 'wall_loss_db', 'N_grid', ...
    'scores_k_map', 'scores_p_map', 'scores_c_map', 'labels_map');
fprintf('計算結果保存（生スコア含む）: raytracing/plateau_results.mat\n');

fprintf('\nPhase 4-③ raytrace_plateau_auc: OK\n');


%% ─────────────────────────────────────────────────────────────────
%% ローカル関数（raytrace_plateau.m と同一）
%% ─────────────────────────────────────────────────────────────────
function nw = count_walls(P0, P1, buildings, aabb)
    dx = P1(1)-P0(1);  dy = P1(2)-P0(2);  dz = P1(3)-P0(3);
    x0 = P0(1);  y0 = P0(2);  z0 = P0(3);
    bldg_h = aabb(:,5);
    mask_h = bldg_h > z0;
    t_zhi  = (bldg_h - z0) / dz;
    t_zhi  = min(t_zhi, 1.0);
    mask_h = mask_h & (t_zhi > 1e-9);
    x_hi = x0 + t_zhi .* dx;
    y_hi = y0 + t_zhi .* dy;
    rx_lo = min(x0, x_hi);  rx_hi = max(x0, x_hi);
    ry_lo = min(y0, y_hi);  ry_hi = max(y0, y_hi);
    mask  = mask_h & ...
            (rx_hi >= aabb(:,1)) & (rx_lo <= aabb(:,2)) & ...
            (ry_hi >= aabb(:,3)) & (ry_lo <= aabb(:,4));
    cands = find(mask);
    nw = 0;
    for ci = 1:numel(cands)
        i  = cands(ci);
        th = t_zhi(i);
        A  = [x0,       y0      ];
        B  = [x0+th*dx, y0+th*dy];
        fp = buildings(i).footprint;
        if seg_x_poly(A, B, fp)
            nw = nw + 1;
        end
    end
end

function hit = seg_x_poly(A, B, poly)
    if inpolygon(A(1),A(2),poly(:,1),poly(:,2)) || ...
       inpolygon(B(1),B(2),poly(:,1),poly(:,2))
        hit = true; return;
    end
    M = size(poly,1);
    for ei = 1:M
        ej = mod(ei,M)+1;
        if seg_x_seg(A, B, poly(ei,:), poly(ej,:))
            hit = true; return;
        end
    end
    hit = false;
end

function r = seg_x_seg(A, B, C, D)
    AB = B-A;  CD = D-C;
    den = AB(1)*CD(2) - AB(2)*CD(1);
    if abs(den) < 1e-10; r=false; return; end
    AC = C-A;
    t  = (AC(1)*CD(2) - AC(2)*CD(1)) / den;
    u  = (AC(1)*AB(2) - AC(2)*AB(1)) / den;
    r  = (t>=0 && t<=1 && u>=0 && u<=1);
end

function snr = fs_snr(d_m, f_hz, ptx_dbm, nf_db, bw_hz)
    c  = 3e8;
    pl = 20*log10(4*pi*d_m*f_hz/c);
    N0 = 10*log10(1.38e-23 * 290 * bw_hz) + 30 + nf_db;
    snr = ptx_dbm - pl - N0;
end

function sc = mk_score(k, p)
    kn = max(0, min(1, (3-k)/2));
    sc = [3-k,  p,  0.5*kn + 0.5*p];
end

function y = awgn_add(x, snr_db, is_real)
    pwr  = mean(abs(x).^2);
    npwr = pwr / 10^(snr_db/10);
    if is_real
        y = x + sqrt(npwr)*randn(size(x));
    else
        y = x + sqrt(npwr/2)*(randn(size(x))+1j*randn(size(x)));
    end
end
