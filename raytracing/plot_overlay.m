% raytracing/plot_overlay.m
% 発表用：識別精度ヒートマップ × PLATEAU建物フットプリント 重ね合わせ図
%
% plateau_results.mat を読み込んで再計算なしで描画する。
% Combined手法の識別精度を背景に、建物輪郭・90%ライン・ETCマーカーを重ねる。
%
% 出力: raytracing/plateau_overlay.png

base_dir = fileparts(which('plot_overlay'));
addpath(fullfile(base_dir,'..','citymodel'));
addpath(fullfile(base_dir,'..'));          % apply_paper_style

%% ─────────────────────────────────────────────────────────────────
%% 計算結果の読み込み
%% ─────────────────────────────────────────────────────────────────
mat_path = fullfile(base_dir, 'plateau_results.mat');
if ~exist(mat_path, 'file')
    error(['plateau_results.mat が見つかりません。\n' ...
           'まず raytrace_plateau.m を実行してください。\n' ...
           '  python run_m.py raytracing/raytrace_plateau.m']);
end
load(mat_path, 'acc_map','n_walls_map','los_mask','nlos_mask', ...
               'x_vec','y_vec','X_grid','Y_grid', ...
               'etc_xy','uav_z','wall_loss_db','N_grid');
fprintf('計算結果読込完了: plateau_results.mat\n');

acc_combined = acc_map(:,:,3);   % Combined手法

%% ─────────────────────────────────────────────────────────────────
%% 建物データ読込（フットプリント輪郭用）
%% ─────────────────────────────────────────────────────────────────
city_dir  = fullfile(base_dir,'..','citymodel','data');
gml_files = dir(fullfile(city_dir,'*.gml'));
if isempty(gml_files)
    error('citymodel/data/ に .gml ファイルが見つかりません');
end
fprintf('建物データ読込中...\n');
[buildings, ~] = load_plateau(fullfile(city_dir, gml_files(1).name));
N_bldg = numel(buildings);
fprintf('  %d棟 読込完了\n', N_bldg);

%% ─────────────────────────────────────────────────────────────────
%% 描画範囲内の建物をフィルタしてNaN区切り輪郭ポリラインを構築
%% ─────────────────────────────────────────────────────────────────
% 格子範囲 + マージン（建物が少しはみ出していても描画）
x_lo = min(x_vec) - 20;  x_hi = max(x_vec) + 20;
y_lo = min(y_vec) - 20;  y_hi = max(y_vec) + 20;

% NaN区切りで全建物輪郭を1本のポリラインに連結（plot 1回で高速描画）
bx = [];  by = [];
n_drawn = 0;
for i = 1:N_bldg
    c  = buildings(i).center;
    if c(1) < x_lo || c(1) > x_hi || c(2) < y_lo || c(2) > y_hi; continue; end
    fp = buildings(i).footprint;   % [M×2]
    % 閉じたポリゴン: 始点を末尾に追加して輪郭を閉じる
    bx = [bx; fp(:,1); fp(1,1); NaN];  %#ok<AGROW>
    by = [by; fp(:,2); fp(1,2); NaN];  %#ok<AGROW>
    n_drawn = n_drawn + 1;
end
fprintf('輪郭描画対象: %d棟\n\n', n_drawn);

%% ─────────────────────────────────────────────────────────────────
%% 描画
%% ─────────────────────────────────────────────────────────────────
fig = figure('Visible','off', 'Color','w', 'Position',[0 0 1400 1200]);
ax  = axes(fig, 'Position',[0.07 0.08 0.80 0.85]);
hold(ax,'on');

% ── 背景: Combined 識別精度ヒートマップ ──────────────────────────
h_img = imagesc(ax, x_vec, y_vec, acc_combined);
colormap(ax, parula(256));
clim(ax, [0.4 1.0]);
axis(ax,'xy');

% ── 建物フットプリント輪郭（NaN区切り一括描画）─────────────────
plot(ax, bx, by, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7);

% ── 90% 精度ライン（等高線）──────────────────────────────────────
[~, hc] = contour(ax, X_grid, Y_grid, acc_combined, [0.9 0.9]);
hc.LineColor = [0.1 0.1 0.1];
hc.LineWidth = 2.5;
hc.LineStyle = '-';

% ── 精度50%ライン（参考: 識別崩壊境界）──────────────────────────
[~, hc2] = contour(ax, X_grid, Y_grid, acc_combined, [0.5 0.5]);
hc2.LineColor = [1.0 0.4 0.1];
hc2.LineWidth = 1.5;
hc2.LineStyle = '--';

% ── ETC位置マーカー（★）─────────────────────────────────────────
plot(ax, etc_xy(1), etc_xy(2), 'p', ...
    'MarkerSize', 20, 'LineWidth', 2.0, ...
    'MarkerFaceColor', [1.0 0.15 0.15], ...
    'MarkerEdgeColor', 'k');
text(ax, etc_xy(1)+18, etc_xy(2)+5, 'ETC', ...
    'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'Margin', 2);

% ── 凡例テキスト（手動配置）──────────────────────────────────────
x_leg = x_vec(end) - 180;
y_leg = y_vec(end) - 30;
dy    = 28;
% 凡例背景ボックス（白・黒枠）
lx1 = x_leg-32;  lx2 = x_leg+178;
ly1 = y_leg-3*dy-16;  ly2 = y_leg+16;
patch(ax, [lx1 lx2 lx2 lx1 lx1], [ly1 ly1 ly2 ly2 ly1], 'w', ...
    'EdgeColor',[0.15 0.15 0.15], 'LineWidth', 1.2, 'HandleVisibility','off');
text(ax, x_leg, y_leg,      '■ 識別精度 (Combined)', 'Color','k','FontSize',10,'FontWeight','bold');
text(ax, x_leg, y_leg-dy,   '── 90% ライン',          'Color','k','FontSize',10);
text(ax, x_leg, y_leg-2*dy, '-- 50% ライン',          'Color',[1.0 0.6 0.3],'FontSize',10);
text(ax, x_leg, y_leg-3*dy, '— 建物輪郭',             'Color',[0.3 0.3 0.3],'FontSize',10);
plot(ax, x_leg-18, y_leg,      's','MarkerSize',8,'MarkerFaceColor','none','MarkerEdgeColor','k');
plot(ax, x_leg-18, y_leg-dy,   '-','Color','k','LineWidth',2.5);
plot(ax, x_leg-18, y_leg-2*dy, '--','Color',[1.0 0.5 0.2],'LineWidth',1.5);
plot(ax, x_leg-18, y_leg-3*dy, '-','Color',[0.45 0.45 0.45],'LineWidth',0.8);

% ── 軸・ラベル設定 ─────────────────────────────────────────────
axis(ax,'equal','tight');
xlim(ax, [x_lo+20, x_hi-20]);
ylim(ax, [y_lo+20, y_hi-20]);
grid(ax,'on');
ax.GridColor     = [0.4 0.4 0.4];
ax.GridAlpha     = 0.4;
ax.GridLineStyle = ':';
ax.TickDir       = 'out';
ax.FontSize      = 11;
ax.Color         = 'w';

xlabel(ax, 'East [m]',  'FontSize', 13);
ylabel(ax, 'North [m]', 'FontSize', 13);
title(ax, sprintf(['識別精度マップ × PLATEAU建物配置\n', ...
    'Combined手法 | UAV高度%dm | ETC(−250,−300) | 壁透過損%ddB/棟 | 格子%d×%d'], ...
    uav_z, wall_loss_db, N_grid, N_grid), 'FontSize', 12);

% カラーバー
cb = colorbar(ax);
cb.Label.String   = '識別精度 Accuracy';
cb.Label.FontSize = 12;
cb.FontSize       = 10;

% 統計注釈（図内）
n_tot = N_grid^2;
ann_str = sprintf('LOS %.1f%%: Acc=1.000\nNLOS %.1f%%: Acc=%.3f\n全体平均: Acc=%.3f', ...
    100*mean(los_mask(:)),  100*mean(nlos_mask(:)), ...
    mean(acc_combined(nlos_mask)), mean(acc_combined(:)));
text(ax, x_vec(1)+15, y_vec(1)+60, ann_str, ...
    'Color','k','FontSize',10,'BackgroundColor','w', ...
    'Margin',5,'VerticalAlignment','bottom','FontName','FixedWidth');

%% ─────────────────────────────────────────────────────────────────
%% 保存
%% ─────────────────────────────────────────────────────────────────
out_path = fullfile(base_dir, 'plateau_overlay.png');
apply_paper_style(fig);
saveas(fig, out_path);
fprintf('保存: raytracing/plateau_overlay.png\n');
disp('plot_overlay: OK');
