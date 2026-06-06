% citymodel/check_etc_pos.m
% ETC設置候補座標が建物フットプリント内か判定し、±50m拡大図を出力する
%
% 判定: inpolygon で全建物フットプリントと照合
%   → 空白(道路)なら確定、建物内ならグリッドサーチで最寄り空白を提案
% 出力: citymodel/etc_pos_check.png

this_dir = fileparts(which('check_etc_pos'));
addpath(this_dir);

%% 候補座標（ここを変更して再実行）
etc_x = -250;
etc_y = -300;

%% GML読込 ─────────────────────────────────────────────────────────
data_dir = fullfile(this_dir, 'data');
files    = dir(fullfile(data_dir, '*.gml'));
if isempty(files); error('citymodel/data/ に .gml ファイルが見つかりません'); end
fprintf('読込: %s\n', files(1).name);
[buildings, info] = load_plateau(fullfile(data_dir, files(1).name));
N = info.n_buildings;

%% 建物内/外 判定 ─────────────────────────────────────────────────
fprintf('\n=== ETC設置座標チェック ===\n');
fprintf('候補点  : (x, y) = (%d, %d) m\n', etc_x, etc_y);
fprintf('建物数  : %d棟  (全棟照合)\n\n', N);

inside_idx = -1;
for i = 1:N
    fp = buildings(i).footprint;
    if inpolygon(etc_x, etc_y, fp(:,1), fp(:,2))
        inside_idx = i;
        break;
    end
end

if inside_idx < 0
    fprintf('[OK] (%d, %d) は建物外（道路/空白）です。ETC設置に問題なし。\n', etc_x, etc_y);
    pos_final  = [etc_x, etc_y];
    status_str = '道路上（空白）';
    is_alt     = false;
else
    b = buildings(inside_idx);
    fprintf('[NG] (%d, %d) は建物内です。\n', etc_x, etc_y);
    fprintf('     建物 #%d  高さ=%.1fm  重心=(%.1f, %.1f)\n', ...
        inside_idx, b.height, b.center(1), b.center(2));
    fprintf('\n最寄り空白座標を探索中（5mグリッド, 半径50m）...\n');
    pos_final  = nearest_open(etc_x, etc_y, buildings, N, 5, 50);
    fprintf('[代替] → (%.0f, %.0f) m\n', pos_final(1), pos_final(2));
    status_str = sprintf('建物内→代替(%.0f,%.0f)', pos_final(1), pos_final(2));
    is_alt     = true;
end

%% 緯度経度変換 ──────────────────────────────────────────────────
R    = 6371000;
lat0 = info.origin_latlon(1);
lon0 = info.origin_latlon(2);
lat_etc = lat0 + pos_final(2) / R * (180/pi);
lon_etc = lon0 + pos_final(1) / (R * cos(lat0*pi/180)) * (180/pi);
fprintf('\n確定座標: (x,y) = (%.0f, %.0f) m\n', pos_final(1), pos_final(2));
fprintf('緯度経度: lat=%.6f°, lon=%.6f°\n', lat_etc, lon_etc);

%% ±50m 拡大図 ────────────────────────────────────────────────────
margin = 50;
cx0 = pos_final(1);  cy0 = pos_final(2);
xlo = cx0 - margin;  xhi = cx0 + margin;
ylo = cy0 - margin;  yhi = cy0 + margin;

% 表示範囲内の建物を抽出（マージン+30m）
idx_plot = [];
for i = 1:N
    c = buildings(i).center;
    if c(1) >= xlo-30 && c(1) <= xhi+30 && c(2) >= ylo-30 && c(2) <= yhi+30
        idx_plot(end+1) = i;  %#ok<AGROW>
    end
end
fprintf('\n拡大図: 描画対象 %d棟\n', numel(idx_plot));

fig = figure('Visible','off', 'Position',[0 0 1200 1200]);
ax  = axes(fig, 'Position',[0.09 0.09 0.85 0.84]);
hold(ax,'on');

for i = idx_plot
    fp = buildings(i).footprint;
    h  = buildings(i).height;
    fc = [0.38 0.42 0.58];
    if isnan(h); fc = [0.75 0.62 0.50]; end
    fill(ax, fp(:,1), fp(:,2), fc, ...
        'EdgeColor',[0.15 0.15 0.25], 'LineWidth', 0.9, 'FaceAlpha', 0.90);
end

% 主グリッド 10m
grid(ax,'on'); box(ax,'on');
ax.GridColor      = [0.80 0.80 0.80];
ax.GridAlpha      = 0.65;
ax.GridLineStyle  = '-';
ax.XTick          = floor(xlo/10)*10 : 10 : ceil(xhi/10)*10;
ax.YTick          = floor(ylo/10)*10 : 10 : ceil(yhi/10)*10;
ax.TickDir        = 'out';
ax.FontSize       = 10;

% 補助グリッド 5m
ax.XMinorTick        = 'on';
ax.YMinorTick        = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridAlpha    = 0.30;

% 候補点マーカー（黄 +）
plot(ax, etc_x, etc_y, '+', 'Color',[1.0 0.9 0.0], 'MarkerSize',22, 'LineWidth',2.5);
text(ax, etc_x+1.5, etc_y+2.0, sprintf('候補 (%d,%d)', etc_x, etc_y), ...
    'Color',[1.0 0.95 0.2], 'FontSize',10, 'FontWeight','bold');

% 確定点マーカー
if is_alt
    % 代替座標（緑 ★）
    plot(ax, pos_final(1), pos_final(2), '*', 'Color',[0.2 1.0 0.4], ...
        'MarkerSize',20, 'LineWidth',2.5);
    text(ax, pos_final(1)+1.5, pos_final(2)-3.5, ...
        sprintf('代替確定 (%.0f,%.0f)', pos_final(1), pos_final(2)), ...
        'Color',[0.2 1.0 0.4], 'FontSize',10, 'FontWeight','bold');
else
    % そのまま確定（赤 ★）
    plot(ax, pos_final(1), pos_final(2), '*', 'Color',[1.0 0.25 0.25], ...
        'MarkerSize',24, 'LineWidth',3.0);
    text(ax, pos_final(1)+1.5, pos_final(2)-3.5, ...
        sprintf('ETC確定 (%d,%d)', etc_x, etc_y), ...
        'Color',[1.0 0.3 0.3], 'FontSize',11, 'FontWeight','bold');
end

% ±20m / ±50m の参照円
th = linspace(0, 2*pi, 200);
plot(ax, cx0+20*cos(th), cy0+20*sin(th), '--', 'Color',[0.7 0.7 0.7], 'LineWidth',0.8);
plot(ax, cx0+50*cos(th), cy0+50*sin(th), '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
text(ax, cx0+21, cy0+0, '±20m', 'Color',[0.7 0.7 0.7], 'FontSize',8);
text(ax, cx0+51, cy0+0, '±50m', 'Color',[0.5 0.5 0.5], 'FontSize',8);

axis(ax,'equal');
xlim(ax, [xlo, xhi]); ylim(ax, [ylo, yhi]);
xlabel(ax, 'East [m]'); ylabel(ax, 'North [m]');
title(ax, sprintf('ETC座標確認  中心(%d,%d) m  ±%dm  |  グリッド10m (補助5m)  |  判定: %s', ...
    etc_x, etc_y, margin, status_str), 'FontSize', 11);

saveas(fig, fullfile(this_dir, 'etc_pos_check.png'));
fprintf('保存: citymodel/etc_pos_check.png\n');

fprintf('\n=== raytrace_plateau.m への設定値 ===\n');
fprintf('  etc_pos = [%.0f, %.0f];   %% ETC路側機 ローカルXY [m]\n', pos_final(1), pos_final(2));
fprintf('  etc_lat = %.6f;  etc_lon = %.6f;  %% 参考: 緯度経度\n', lat_etc, lon_etc);

disp('check_etc_pos: OK');


%% ─── ローカル関数: 最寄り空白座標サーチ ─────────────────────────

function pos = nearest_open(x0, y0, buildings, N, step, max_r)
% step[m] グリッドで max_r[m] 以内の最寄り建物外座標を返す

    % 探索範囲を絞った近傍建物リスト（速度向上）
    margin = max_r + 60;
    near = false(N,1);
    for i = 1:N
        c = buildings(i).center;
        if abs(c(1)-x0) <= margin && abs(c(2)-y0) <= margin
            near(i) = true;
        end
    end
    near_idx = find(near);

    best_d2  = Inf;
    best_pos = [x0, y0];

    for dx = -max_r : step : max_r
        for dy = -max_r : step : max_r
            d2 = dx^2 + dy^2;
            if d2 > max_r^2 || d2 >= best_d2; continue; end
            cx = x0 + dx;
            cy = y0 + dy;
            in_any = false;
            for i = near_idx'
                fp = buildings(i).footprint;
                if inpolygon(cx, cy, fp(:,1), fp(:,2))
                    in_any = true; break;
                end
            end
            if ~in_any
                best_d2  = d2;
                best_pos = [cx, cy];
            end
        end
    end
    pos = best_pos;
end
