% citymodel/find_road_space.m
% PLATEAU建物配置から「道路と思われる開けた空間」を目視確認するための座標付き2Dマップ出力
%
% 出力:
%   citymodel/road_map_overview.png  -- 全体俯瞰 (100mグリッド)
%   citymodel/road_map_center.png   -- 中心部拡大 ±350m (50mグリッド)

this_dir = fileparts(which('find_road_space'));
addpath(this_dir);

%% GML読込 ─────────────────────────────────────────────────────────
data_dir = fullfile(this_dir, 'data');
files    = dir(fullfile(data_dir, '*.gml'));
if isempty(files)
    error('citymodel/data/ に .gml ファイルが見つかりません');
end
fprintf('読込: %s\n', files(1).name);
[buildings, info] = load_plateau(fullfile(data_dir, files(1).name));
N = info.n_buildings;

fprintf('\n建物数      : %d棟\n', N);
fprintf('エリア(E)   : %.0f〜%.0f m  (幅 %.0f m)\n', info.x_range(1), info.x_range(2), diff(info.x_range));
fprintf('エリア(N)   : %.0f〜%.0f m  (幅 %.0f m)\n', info.y_range(1), info.y_range(2), diff(info.y_range));
fprintf('座標原点    : lat=%.6f°, lon=%.6f°\n\n', info.origin_latlon(1), info.origin_latlon(2));

%% ─────────────────────────────────────────────────────────────────
%% 図1: 全体俯瞰（100mグリッド）
%% ─────────────────────────────────────────────────────────────────
fig1 = figure('Visible','off', 'Position',[0 0 1600 1200]);
ax1  = axes(fig1, 'Position',[0.07 0.08 0.88 0.87]);
hold(ax1,'on');

for i = 1:N
    fp = buildings(i).footprint;
    fill(ax1, fp(:,1), fp(:,2), [0.52 0.54 0.64], 'EdgeColor','none', 'FaceAlpha',0.85);
end

% 原点十字
plot(ax1, 0, 0, 'r+', 'MarkerSize',14, 'LineWidth',2.5);
text(ax1, 20, 20, '(0,0)', 'Color','r', 'FontSize',9, 'FontWeight','bold');

grid(ax1,'on'); box(ax1,'on');
ax1.GridColor      = [0.15 0.15 0.15];
ax1.GridAlpha      = 0.60;
ax1.GridLineStyle  = '-';
x_lo = floor(info.x_range(1)/100)*100;
x_hi = ceil( info.x_range(2)/100)*100;
y_lo = floor(info.y_range(1)/100)*100;
y_hi = ceil( info.y_range(2)/100)*100;
ax1.XTick = x_lo:100:x_hi;
ax1.YTick = y_lo:100:y_hi;
ax1.TickDir = 'out';
ax1.FontSize = 8;
ax1.XTickLabelRotation = 0;

axis(ax1,'equal');
xlim(ax1, [info.x_range(1)-30, info.x_range(2)+30]);
ylim(ax1, [info.y_range(1)-30, info.y_range(2)+30]);
xlabel(ax1, sprintf('East [m]   ※原点: lat=%.5f°, lon=%.5f°', info.origin_latlon(1), info.origin_latlon(2)));
ylabel(ax1, 'North [m]');
title(ax1, sprintf('PLATEAU 建物配置 全体図  N=%d棟  |  グリッド 100m', N), 'FontSize',11);

saveas(fig1, fullfile(this_dir, 'road_map_overview.png'));
fprintf('保存: citymodel/road_map_overview.png  (全体俯瞰 100mグリッド)\n');
close(fig1);

%% ─────────────────────────────────────────────────────────────────
%% 図2: 中心部拡大（±350m, 50mグリッド + 10m補助線）
%% ─────────────────────────────────────────────────────────────────
Z = 350;   % 表示範囲 ±Z [m]

fig2 = figure('Visible','off', 'Position',[0 0 1600 1400]);
ax2  = axes(fig2, 'Position',[0.07 0.07 0.88 0.88]);
hold(ax2,'on');

n_drawn = 0;
for i = 1:N
    c = buildings(i).center;
    if abs(c(1)) > Z+80 || abs(c(2)) > Z+80; continue; end
    fp = buildings(i).footprint;
    h  = buildings(i).height;
    if isnan(h)
        fc = [0.75 0.65 0.55];  % 高さ不明 → 茶色で区別
    else
        fc = [0.45 0.48 0.60];  % 通常建物 → 青灰
    end
    fill(ax2, fp(:,1), fp(:,2), fc, ...
        'EdgeColor',[0.20 0.20 0.30], 'LineWidth',0.5, 'FaceAlpha',0.88);
    n_drawn = n_drawn + 1;
end

% 原点十字
plot(ax2, 0, 0, 'r+', 'MarkerSize',18, 'LineWidth',3);
text(ax2, 8, 8, '(0,0)', 'Color','r', 'FontSize',11, 'FontWeight','bold');

% 主グリッド 50m
grid(ax2,'on'); box(ax2,'on');
ax2.GridColor      = [0.10 0.10 0.10];
ax2.GridAlpha      = 0.65;
ax2.GridLineStyle  = '-';
ax2.XTick          = -500:50:500;
ax2.YTick          = -500:50:500;
ax2.TickDir        = 'out';
ax2.FontSize       = 9;

% 補助グリッド 10m（MinorTick）
ax2.XMinorTick       = 'on';
ax2.YMinorTick       = 'on';
ax2.MinorGridLineStyle = ':';
ax2.MinorGridAlpha   = 0.30;

axis(ax2,'equal');
xlim(ax2, [-Z, Z]);
ylim(ax2, [-Z, Z]);
xlabel(ax2, 'East [m]');
ylabel(ax2, 'North [m]');
title(ax2, sprintf(['中心部拡大 ±%dm  |  グリッド 50m (補助 10m)\n', ...
    '青灰: 建物  茶: 高さ不明棟  白/薄色: 道路・空白スペース  (描画%d棟)'], Z, n_drawn), ...
    'FontSize', 11);

saveas(fig2, fullfile(this_dir, 'road_map_center.png'));
fprintf('保存: citymodel/road_map_center.png  (中心部 ±%dm 50mグリッド)\n', Z);
close(fig2);

%% ─────────────────────────────────────────────────────────────────
%% 図3: 幹線道路拡大図（x=-400〜0, y=-350〜-150, 50m+10mグリッド）
%% ─────────────────────────────────────────────────────────────────
xlo3 = -420;  xhi3 =  20;   % 描画マージン付き
ylo3 = -370;  yhi3 = -130;

% 表示範囲（軸）
xview = [-400,  0];
yview = [-350, -150];

% 表示範囲内の建物を抽出
in3 = false(N,1);
for i = 1:N
    c = buildings(i).center;
    if c(1) >= xlo3 && c(1) <= xhi3 && c(2) >= ylo3 && c(2) <= yhi3
        in3(i) = true;
    end
end
idx3 = find(in3);
fprintf('幹線道路拡大: 描画対象 %d棟\n', numel(idx3));

% 幅=400m, 高さ=200m → 縦横比 2:1 → 1600×900
fig3 = figure('Visible','off', 'Position',[0 0 1600 900]);
ax3  = axes(fig3, 'Position',[0.07 0.10 0.88 0.83]);
hold(ax3,'on');

for i = idx3'
    fp = buildings(i).footprint;
    h  = buildings(i).height;
    fc = [0.40 0.43 0.58];
    if isnan(h); fc = [0.75 0.62 0.50]; end
    fill(ax3, fp(:,1), fp(:,2), fc, ...
        'EdgeColor',[0.15 0.15 0.25], 'LineWidth',0.6, 'FaceAlpha',0.90);
end

% 主グリッド 50m
grid(ax3,'on'); box(ax3,'on');
ax3.GridColor     = [0.85 0.85 0.85];
ax3.GridAlpha     = 0.70;
ax3.GridLineStyle = '-';
ax3.XTick         = -450:50:50;
ax3.YTick         = -400:50:-100;
ax3.TickDir       = 'out';
ax3.FontSize      = 10;

% 補助グリッド 10m
ax3.XMinorTick        = 'on';
ax3.YMinorTick        = 'on';
ax3.MinorGridLineStyle = ':';
ax3.MinorGridAlpha    = 0.35;

% グリッド交点に座標ラベル（100m刻みで表示、視認しやすいよう間引き）
for gx = -400:100:0
    for gy = -350:100:-150
        text(ax3, gx+4, gy+4, sprintf('(%d,%d)', gx, gy), ...
            'FontSize',7, 'Color',[0.9 0.9 0.3], 'FontWeight','bold');
    end
end

axis(ax3,'equal');
xlim(ax3, xview);
ylim(ax3, yview);
xlabel(ax3, 'East [m]');
ylabel(ax3, 'North [m]');
title(ax3, sprintf(['幹線道路拡大  x=[%d, %d] m, y=[%d, %d] m  |  グリッド 50m (補助 10m)\n', ...
    '黄ラベル: 100m交点座標 / 青灰: 建物 / 茶: 高さ不明 / 白黒: 道路・空白'], ...
    xview(1),xview(2),yview(1),yview(2)), 'FontSize', 11);

saveas(fig3, fullfile(this_dir, 'road_map_zoom.png'));
fprintf('保存: citymodel/road_map_zoom.png  (幹線道路拡大 x=[%d,%d] y=[%d,%d])\n', ...
    xview(1),xview(2),yview(1),yview(2));
close(fig3);

%% ─────────────────────────────────────────────────────────────────
%% コンソール: 道路座標の読み取り方ガイド
%% ─────────────────────────────────────────────────────────────────
fprintf('\n=== ETC設置座標の決め方 ===\n');
fprintf('road_map_center.png を画像ビューアで開き、拡大して確認してください。\n');
fprintf('\n【道路（ETC設置候補）の見分け方】\n');
fprintf('  白い隙間 (建物なし) が帯状に続く箇所 → 道路\n');
fprintf('  幅 5〜8m 程度の直線的な空白 → 日本の一般道\n');
fprintf('  幅 15〜25m 程度の広い空白  → 幹線道路 (ETCに適)\n');
fprintf('\n【座標→緯度経度 変換式】\n');
fprintf('  lat = %.6f + y_m / 6371000 * (180/pi)\n', info.origin_latlon(1));
fprintf('  lon = %.6f + x_m / (6371000 * cos(%.4f*pi/180)) * (180/pi)\n', ...
    info.origin_latlon(2), info.origin_latlon(1));
fprintf('\n【次ステップ】\n');
fprintf('  特定した道路座標を raytrace_plateau.m の etc_pos に設定して\n');
fprintf('  フェーズ4③ のレイトレーシングシミュレーションを実行する。\n');

disp('find_road_space: OK');
