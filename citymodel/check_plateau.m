% citymodel/check_plateau.m
% PLATEAU建物 3D目視確認スクリプト
% 実行: check_plateau
%
% load_plateau で読んだ建物フットプリントを高さ方向に押し出して 3D patch 表示。
% MATLAB GUI でマウス回転・ズーム可能。
%
% 描画範囲は plot_overlay.m と同一座標系 (plateau_results.mat の x_vec/y_vec)。

this_dir  = fileparts(which('check_plateau'));
addpath(this_dir);                   % load_plateau.m と同じフォルダ

%% 1. GML読込 ──────────────────────────────────────────────────────
data_dir = fullfile(this_dir, 'data');
files    = dir(fullfile(data_dir, '*.gml'));
if isempty(files)
    error('citymodel/data/ に .gml ファイルが見つかりません');
end
fprintf('読込: %s\n', files(1).name);
[buildings, info] = load_plateau(fullfile(data_dir, files(1).name));
fprintf('読込完了: %d棟\n\n', info.n_buildings);

%% 2. 描画範囲: plateau_results.mat を単一出どころとする ─────────────
% → plot_overlay.m も同ファイルから x_vec/y_vec を読むため座標系が一致する
mat_path = fullfile(this_dir, '..', 'raytracing', 'plateau_results.mat');
if exist(mat_path, 'file')
    tmp    = load(mat_path, 'x_vec','y_vec','etc_xy');
    x_vec  = tmp.x_vec;
    y_vec  = tmp.y_vec;
    etc_xy = tmp.etc_xy;
    fprintf('描画範囲: plateau_results.mat から取得\n');
else
    % raytrace_plateau.m と同一パラメータ（フォールバック）
    etc_xy     = [-250, -300];
    grid_range = 500;
    N_grid     = 50;
    x_vec = linspace(etc_xy(1)-grid_range, etc_xy(1)+grid_range, N_grid);
    y_vec = linspace(etc_xy(2)-grid_range, etc_xy(2)+grid_range, N_grid);
    fprintf('描画範囲: フォールバック（plateau_results.mat 未検出）\n');
end

% plot_overlay.m と完全一致する xlim/ylim の導出
x_lo = min(x_vec) - 20;  x_hi = max(x_vec) + 20;
y_lo = min(y_vec) - 20;  y_hi = max(y_vec) + 20;
xlim_range = [min(x_vec), max(x_vec)];   % = [x_lo+20, x_hi-20]
ylim_range = [min(y_vec), max(y_vec)];

fprintf('xlim = [%.0f, %.0f]\n', xlim_range);
fprintf('ylim = [%.0f, %.0f]\n', ylim_range);

%% 3. 表示対象: 描画範囲内の建物（plot_overlay.m と同一フィルタ）──────
N_bldg = info.n_buildings;
idx = [];
for i = 1:N_bldg
    c = buildings(i).center;
    if c(1) >= x_lo && c(1) <= x_hi && c(2) >= y_lo && c(2) <= y_hi
        idx(end+1) = i;  %#ok<AGROW>
    end
end
N_plot = numel(idx);
fprintf('\n描画対象: %d棟\n', N_plot);

%% 4. 高さ確定（NaN / 非正 → 3m 仮置き）─────────────────────────
hh = arrayfun(@(i) buildings(i).height, idx);
nan_m = isnan(hh) | hh <= 0;
hh(nan_m) = 3;
h_max = max(hh);  if h_max < 1; h_max = 1; end

fprintf('高さ %.1f〜%.1f m  (高さ不明 %d棟 → 3m 仮置き)\n', ...
    min(hh), max(hh), sum(nan_m));

cmap = parula(256);   % 低→高 = 青→黄

%% 5. 頂点・面データを事前収集（patch は上面/側面の2回だけ）────────
top_verts = [];  top_faces = cell(N_plot,1);  top_c = zeros(N_plot,3);
sid_verts = [];  sid_faces = [];  sid_c = [];
tv_base = 0;  sv_base = 0;

for k = 1:N_plot
    bi  = idx(k);
    fp  = buildings(bi).footprint;   % [M×2] ローカルXY [m]
    h   = hh(k);
    N_v = size(fp, 1);

    ci       = max(1, min(256, round(h/h_max*255)+1));
    col      = cmap(ci, :);
    col_side = col * 0.62;           % 側面を暗くして立体感

    % ── 上面 ────────────────────────────────────────
    top_verts          = [top_verts; fp, h*ones(N_v,1)];  %#ok<AGROW>
    top_faces{k}       = tv_base + (1:N_v);
    top_c(k,:)         = col;
    tv_base            = tv_base + N_v;

    % ── 側面クワッド × N_v 枚 ─────────────────────
    for vi = 1:N_v
        vj = mod(vi, N_v) + 1;
        sv = [fp(vi,1) fp(vi,2) 0;
              fp(vj,1) fp(vj,2) 0;
              fp(vj,1) fp(vj,2) h;
              fp(vi,1) fp(vi,2) h];
        sid_verts = [sid_verts; sv];                       %#ok<AGROW>
        sid_faces = [sid_faces; sv_base + [1 2 3 4]];     %#ok<AGROW>
        sid_c     = [sid_c; col_side];                     %#ok<AGROW>
        sv_base   = sv_base + 4;
    end
end

% NaN-padded face matrix（MATLAB patch の可変頂点数対応）
max_nv = max(cellfun(@numel, top_faces));
top_F  = NaN(N_plot, max_nv);
for k = 1:N_plot
    fi = top_faces{k};
    top_F(k, 1:numel(fi)) = fi;
end

fprintf('描画データ準備完了 (上面 %d面, 側面 %d面)\n', N_plot, size(sid_faces,1));

%% 6. 3D Figure ────────────────────────────────────────────────────
fig = figure('Name','PLATEAU 3D建物確認 — マウスで回転・ズーム可', ...
             'NumberTitle','off', ...
             'Position', [100 100 1200 800]);
ax = axes(fig, 'Position', [0.05 0.08 0.82 0.86]);
hold(ax,'on');  grid(ax,'on');  box(ax,'on');
xlabel(ax,'East [m]');  ylabel(ax,'North [m]');  zlabel(ax,'高さ [m]');
title(ax, sprintf('PLATEAU 3D建物  N=%d棟  xlim=[%.0f,%.0f]  ylim=[%.0f,%.0f]', ...
    N_plot, xlim_range(1), xlim_range(2), ylim_range(1), ylim_range(2)));
view(ax, 45, 30);

% 側面 (一括 patch)
patch(ax, ...
    'Vertices',       sid_verts, ...
    'Faces',          sid_faces, ...
    'FaceColor',      'flat', ...
    'FaceVertexCData', sid_c, ...
    'EdgeColor',      [0.25 0.25 0.25], ...
    'LineWidth',      0.3, ...
    'FaceAlpha',      0.88);

% 上面 (一括 patch — 最後に描いてエッジが前面に出る)
patch(ax, ...
    'Vertices',       top_verts, ...
    'Faces',          top_F, ...
    'FaceColor',      'flat', ...
    'FaceVertexCData', top_c, ...
    'EdgeColor',      [0.25 0.25 0.25], ...
    'LineWidth',      0.5, ...
    'FaceAlpha',      0.95);

%% 7. 仕上げ・軸設定 ────────────────────────────────────────────────
colormap(ax, parula);
clim(ax, [0 h_max]);
cb = colorbar(ax);  cb.Label.String = '建物高さ [m]';

lighting(ax, 'gouraud');
light(ax, 'Position', [400 400 600], 'Style', 'infinite');
material(ax, 'dull');

% plot_overlay.m と一致する XY 範囲・アスペクト比を設定
axis(ax, 'equal');
xlim(ax, xlim_range);   % [min(x_vec), max(x_vec)]
ylim(ax, ylim_range);   % [min(y_vec), max(y_vec)]
set(ax, 'ZLim', [0, h_max * 1.1]);
% YDir: 3D axes では 'normal' が既定（Y 上向き）。axis equal 後も維持される。

fprintf('\n3D表示完了。\n');
fprintf('  回転 : ツールバー [Rotate 3D] ボタン → マウスドラッグ\n');
fprintf('  ズーム: スクロールホイール\n');
fprintf('  パン  : Shift + ドラッグ\n');
