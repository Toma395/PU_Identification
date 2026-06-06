% raytracing/raytrace_plateau.m
% フェーズ4③ ステップ1: PLATEAU建物配置 × 幾何的LOS/NLOS伝搬モデル
%
% 伝搬モデル (ステップ1: 遮蔽・透過のみ):
%   LOS  (貫通 0棟): 自由空間パスロス (Friis, 3D経路長)
%   NLOS (貫通 n棟): 自由空間パスロス + wall_loss_db × n  (コンクリート透過)
%
% 3D遮蔽判定アルゴリズム:
%   ETC(z=5m)→UAV(z=100m) 線分に対し:
%   ① 高さスラブ判定 (h>z_ETC の建物のみ, ベクトル化)
%   ② XY-AABB判定 (ベクトル化, t∈[0,t_zhi] の包絡ボックス)
%   ③ 2Dポリゴン交差判定 (候補建物のみ)
%
% 出力:
%   raytracing/plateau_los_map.png
%   raytracing/plateau_walls_map.png
%   raytracing/plateau_snr_map.png
%   raytracing/plateau_acc_kurtosis.png
%   raytracing/plateau_acc_preamble.png
%   raytracing/plateau_acc_combined.png

base_dir = fileparts(which('raytrace_plateau'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..','citymodel'));

%% ─────────────────────────────────────────────────────────────────
%% パラメータ
%% ─────────────────────────────────────────────────────────────────
% ETC路側機 (PU)
etc_xy       = [-250, -300];  % ローカルXY [m]
etc_z        =  5;            % アンテナ高さ [m]
f_etc        = 5.80e9;        % Hz (ARIB STD-T75)
ptx_etc      = 10;            % dBm
bw_etc       = 4.4e6;         % Hz

% UAV (SU)
f_uav        = 5.82e9;        % Hz (DSRCサブバンド, プランB)
ptx_uav      = 20;            % dBm
bw_uav       = 20e6;          % Hz
uav_z        = 100;           % 飛行高度 [m]

% 共通受信機
nf_db        = 5;             % dB 雑音指数
wall_loss_db = 12;            % dB/棟 コンクリート壁透過損 (文献値)

% UAV格子 (ETC中心 ±grid_range m, N_grid × N_grid 点)
N_grid       = 50;
grid_range   = 500;           % m

% 識別試行数 (クラスあたり; ETC×N + UAV×N = 2N 試行/格子点)
N_trials     = 30;

fprintf('=== raytrace_plateau (Phase4-③ Step1) ===\n');
fprintf('ETC: (%.0f, %.0f) m, z=%.0fm, f=%.2fGHz, Ptx=%.0fdBm\n', ...
    etc_xy(1),etc_xy(2), etc_z, f_etc/1e9, ptx_etc);
fprintf('UAV: z=%.0fm, 格子%d×%d, ±%.0fm, N_trials=%d/クラス\n', ...
    uav_z, N_grid, N_grid, grid_range, N_trials);
fprintf('壁透過損: %.0f dB/棟\n\n', wall_loss_db);

%% ─────────────────────────────────────────────────────────────────
%% 建物データ読込 & AABB 事前計算
%% ─────────────────────────────────────────────────────────────────
city_dir  = fullfile(base_dir,'..','citymodel','data');
gml_files = dir(fullfile(city_dir,'*.gml'));
if isempty(gml_files); error('citymodel/data/ にGMLファイルが見つかりません'); end

fprintf('建物データ読込中...\n');
[buildings, info] = load_plateau(fullfile(city_dir, gml_files(1).name));
N_bldg = info.n_buildings;
fprintf('  %d棟 読込完了\n', N_bldg);

% AABB 行列 [N_bldg × 5]: xmin xmax ymin ymax height
aabb = zeros(N_bldg, 5);
for i = 1:N_bldg
    fp = buildings(i).footprint;
    h  = buildings(i).height;
    if isnan(h); h = 0; end
    aabb(i,:) = [min(fp(:,1)), max(fp(:,1)), min(fp(:,2)), max(fp(:,2)), h];
end
fprintf('  AABB 事前計算完了\n\n');

%% ─────────────────────────────────────────────────────────────────
%% UAV格子 & テンプレート準備
%% ─────────────────────────────────────────────────────────────────
x_vec = linspace(etc_xy(1)-grid_range, etc_xy(1)+grid_range, N_grid);
y_vec = linspace(etc_xy(2)-grid_range, etc_xy(2)+grid_range, N_grid);
[X_grid, Y_grid] = meshgrid(x_vec, y_vec);  % [N_grid×N_grid], row=Y, col=X

[~, ep] = gen_ETC('seed', 1);
uw  = ep.uw;
sps = ep.sps;
thresh = [1.0, 0.5, 0.5];  % kurtosis / preamble / combined 判定閾値

%% ─────────────────────────────────────────────────────────────────
%% メインループ
%% ─────────────────────────────────────────────────────────────────
n_walls_map = zeros(N_grid, N_grid);   % 貫通建物数
snr_etc_map = zeros(N_grid, N_grid);   % ETC受信SNR [dB]
acc_map     = zeros(N_grid, N_grid, 3);% 識別精度 (3手法)

etc_pos3 = [etc_xy, etc_z];

fprintf('遮蔽判定 + 識別実験開始 (%d×%d格子, %d試行/クラス)...\n', ...
    N_grid, N_grid, N_trials);
t0 = tic;

for gi = 1:N_grid          % row (Y方向)
    for gj = 1:N_grid      % col (X方向)
        uav_pos3 = [X_grid(gi,gj), Y_grid(gi,gj), uav_z];

        % ① 遮蔽建物数
        nw = count_walls(etc_pos3, uav_pos3, buildings, aabb);
        n_walls_map(gi,gj) = nw;

        % ② SNR (3D経路長 + 壁透過損)
        d3d    = norm(uav_pos3 - etc_pos3);
        snr_e  = fs_snr(d3d, f_etc, ptx_etc, nf_db, bw_etc) - wall_loss_db*nw;
        snr_u  = fs_snr(d3d, f_uav, ptx_uav, nf_db, bw_uav) - wall_loss_db*nw;
        snr_etc_map(gi,gj) = snr_e;

        % ③ 識別実験
        scores = zeros(N_trials*2, 3);
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
    end

    if mod(gi,10)==0
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
n_tot = N_grid^2;

fprintf('=== LOS/NLOS 分布 ===\n');
fprintf('  LOS  (0棟): %4d点 (%.1f%%)\n', sum(los_mask(:)),  100*mean(los_mask(:)));
fprintf('  NLOS (1棟): %4d点 (%.1f%%)\n', sum(n_walls_map(:)==1), 100*mean(n_walls_map(:)==1));
fprintf('  NLOS (2棟): %4d点 (%.1f%%)\n', sum(n_walls_map(:)==2), 100*mean(n_walls_map(:)==2));
fprintf('  NLOS (3+棟):%4d点 (%.1f%%)\n', sum(n_walls_map(:)>=3), 100*mean(n_walls_map(:)>=3));

fprintf('\n=== 識別精度 ===\n');
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
fprintf('  LOS  : 中央値=%.1fdB  10%%ile=%.1fdB  最大=%.1fdB\n', ...
    median(snr_los),  prctile(snr_los,10),  max(snr_los));
fprintf('  NLOS : 中央値=%.1fdB  10%%ile=%.1fdB  最大=%.1fdB\n', ...
    median(snr_nlos), prctile(snr_nlos,10), max(snr_nlos));

%% ─────────────────────────────────────────────────────────────────
%% プロット出力
%% ─────────────────────────────────────────────────────────────────
fa = {'Visible','off','Position',[0 0 1000 900]};
ap = [0.08 0.08 0.80 0.83];

% ① LOS/NLOSマップ (0=NLOS暗, 1=LOS明)
fig = figure(fa{:});  ax = axes(fig,'Position',ap);
imagesc(ax, x_vec, y_vec, double(los_mask));
colormap(ax, [0.12 0.15 0.40; 0.95 0.88 0.35]);
axis(ax,'xy'); axis(ax,'equal','tight'); grid(ax,'on');
hold(ax,'on');
plot(ax, etc_xy(1), etc_xy(2), 'r*','MarkerSize',14,'LineWidth',2.5);
text(ax, etc_xy(1)+12, etc_xy(2)-12,'ETC','Color','r','FontSize',10,'FontWeight','bold');
xlabel(ax,'East [m]'); ylabel(ax,'North [m]');
title(ax, sprintf('LOS/NLOSマップ (UAV z=%dm)  黄:LOS %.1f%% / 青:NLOS %.1f%%', ...
    uav_z, 100*mean(los_mask(:)), 100*mean(nlos_mask(:))), 'FontSize',11);
cb=colorbar(ax); cb.Ticks=[0.25 0.75]; cb.TickLabels={'NLOS','LOS'};
saveas(fig, fullfile(base_dir,'plateau_los_map.png')); close(fig);
fprintf('保存: raytracing/plateau_los_map.png\n');

% ② 貫通建物数マップ
fig = figure(fa{:});  ax = axes(fig,'Position',ap);
nw_clamp = min(n_walls_map, 5);
imagesc(ax, x_vec, y_vec, nw_clamp);
colormap(ax, [0 0 0; parula(5)]);   % 0棟=黒, 1〜5棟=parula
clim(ax,[0 5]); axis(ax,'xy'); axis(ax,'equal','tight'); grid(ax,'on');
hold(ax,'on');
plot(ax, etc_xy(1), etc_xy(2), 'w*','MarkerSize',14,'LineWidth',2.5);
text(ax, etc_xy(1)+12, etc_xy(2)-12,'ETC','Color','w','FontSize',10,'FontWeight','bold');
xlabel(ax,'East [m]'); ylabel(ax,'North [m]');
title(ax, sprintf('貫通建物数マップ (壁透過損 %ddB/棟)', wall_loss_db), 'FontSize',11);
cb=colorbar(ax); cb.Label.String='貫通建物数 [棟] (5+は5に丸め)';
saveas(fig, fullfile(base_dir,'plateau_walls_map.png')); close(fig);
fprintf('保存: raytracing/plateau_walls_map.png\n');

% ③ SNRマップ
fig = figure(fa{:});  ax = axes(fig,'Position',ap);
imagesc(ax, x_vec, y_vec, snr_etc_map);
colormap(ax, turbo(256)); axis(ax,'xy'); axis(ax,'equal','tight'); grid(ax,'on');
hold(ax,'on');
contour(ax, X_grid, Y_grid, snr_etc_map, [0 0], 'w-','LineWidth',1.5);   % 0dBライン
plot(ax, etc_xy(1), etc_xy(2), 'w*','MarkerSize',14,'LineWidth',2.5);
text(ax, etc_xy(1)+12, etc_xy(2)-12,'ETC','Color','w','FontSize',10,'FontWeight','bold');
xlabel(ax,'East [m]'); ylabel(ax,'North [m]');
title(ax,'ETC受信SNRマップ [dB]  (白線: 0 dB)','FontSize',11);
cb=colorbar(ax); cb.Label.String='受信SNR [dB]';
saveas(fig, fullfile(base_dir,'plateau_snr_map.png')); close(fig);
fprintf('保存: raytracing/plateau_snr_map.png\n');

% ④⑤⑥ 識別精度ヒートマップ (3手法)
pnames = {'plateau_acc_kurtosis.png','plateau_acc_preamble.png','plateau_acc_combined.png'};
for mi = 1:3
    a_mi = acc_map(:,:,mi);
    fig = figure(fa{:});  ax = axes(fig,'Position',ap);
    imagesc(ax, x_vec, y_vec, a_mi);
    colormap(ax, parula(256)); clim(ax,[0.4 1.0]);
    axis(ax,'xy'); axis(ax,'equal','tight'); grid(ax,'on');
    hold(ax,'on');
    contour(ax, X_grid, Y_grid, a_mi, [0.9 0.9], 'w-','LineWidth',2.0);  % 90%ライン
    plot(ax, etc_xy(1), etc_xy(2), 'r*','MarkerSize',14,'LineWidth',2.5);
    text(ax, etc_xy(1)+12, etc_xy(2)-12,'ETC','Color','r','FontSize',10,'FontWeight','bold');
    xlabel(ax,'East [m]'); ylabel(ax,'North [m]');
    title(ax, sprintf('%s\n全体=%.3f  LOS=%.3f  NLOS=%.3f  (白線:90%%ライン)', ...
        mnames{mi}, mean(a_mi(:)), mean(a_mi(los_mask)), mean(a_mi(nlos_mask))), ...
        'FontSize',10);
    cb=colorbar(ax); cb.Label.String='Accuracy';
    saveas(fig, fullfile(base_dir, pnames{mi})); close(fig);
    fprintf('保存: raytracing/%s\n', pnames{mi});
end

fprintf('\nPhase 4-③ raytrace_plateau (Step1): OK\n');


%% ─────────────────────────────────────────────────────────────────
%% ローカル関数
%% ─────────────────────────────────────────────────────────────────

function nw = count_walls(P0, P1, buildings, aabb)
% 3D線分 P0→P1 が貫通する建物数をカウント (AABB高速化付き)
    dx = P1(1)-P0(1);  dy = P1(2)-P0(2);  dz = P1(3)-P0(3);
    x0 = P0(1);  y0 = P0(2);  z0 = P0(3);

    bldg_h = aabb(:,5);

    % ─ Step1: 高さスラブフィルタ (ベクトル化) ─────────────────────
    % z(t)=z0+t*dz ≤ bldg_h  → t ≤ (bldg_h - z0)/dz
    % z0=5m, dz=95m → 単調増加; h>z0 の建物のみが候補
    mask_h    = bldg_h > z0;
    t_zhi     = (bldg_h - z0) / dz;         % [N×1] z(t)=bldg_h となるt
    t_zhi     = min(t_zhi, 1.0);
    mask_h    = mask_h & (t_zhi > 1e-9);

    % ─ Step2: XY-AABB フィルタ (ベクトル化) ───────────────────────
    % t∈[0, t_zhi] の XY 軌跡 包絡ボックス vs 建物AABB
    x_hi = x0 + t_zhi .* dx;
    y_hi = y0 + t_zhi .* dy;

    rx_lo = min(x0, x_hi);  rx_hi = max(x0, x_hi);
    ry_lo = min(y0, y_hi);  ry_hi = max(y0, y_hi);

    mask = mask_h & ...
           (rx_hi >= aabb(:,1)) & (rx_lo <= aabb(:,2)) & ...
           (ry_hi >= aabb(:,3)) & (ry_lo <= aabb(:,4));

    cands = find(mask);

    % ─ Step3: 2Dポリゴン交差判定 ──────────────────────────────────
    nw = 0;
    for ci = 1:numel(cands)
        i  = cands(ci);
        th = t_zhi(i);
        A  = [x0,        y0       ];
        B  = [x0+th*dx,  y0+th*dy ];
        fp = buildings(i).footprint;
        if seg_x_poly(A, B, fp)
            nw = nw + 1;
        end
    end
end

function hit = seg_x_poly(A, B, poly)
% 2D線分 A-B がポリゴン poly と交差するか
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
% 2D線分 A-B と C-D が交差するか (外積法)
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
