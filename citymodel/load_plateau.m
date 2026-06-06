function [buildings, info] = load_plateau(gml_file)
% PLATEAU CityGML LOD1 から建物フットプリントと高さを抽出する
%
% 使い方:
%   load_plateau                   % 引数なし: citymodel/data/*.gml を自動検索して自己テスト
%   [b, info] = load_plateau(path) % 指定ファイルを読み込んでデータ返却
%
% 出力:
%   buildings(i).footprint  [M×2] ローカルXY座標 [m] (East/North, 閉じていない)
%   buildings(i).height     建物高さ [m] (bldg:measuredHeight)
%   buildings(i).base_z     底面標高 [m] (海面上)
%   buildings(i).center     [1×2] 重心 [m]
%   info.n_buildings, .origin_latlon, .x_range, .y_range, .height_range
%
% 座標系: EPSG:6697 (JGD2011 地理3D) → Equirectangular でローカルXY [m] に変換

if nargin == 0
    % ─── 自己テストモード ─────────────────────────────────────────
    this_dir = fileparts(mfilename('fullpath'));
    data_dir = fullfile(this_dir, 'data');
    gml_files = dir(fullfile(data_dir, '*.gml'));
    if isempty(gml_files)
        error('citymodel/data/ にGMLファイルが見つかりません。');
    end
    fname = gml_files(1).name;
    gml_path = fullfile(data_dir, fname);
    fprintf('対象ファイル: %s (%.1f MB)\n\n', fname, dir(gml_path).bytes/1e6);

    [buildings, info] = load_plateau(gml_path);

    % ── 結果表示 ────────────────────────────────────────────────
    fprintf('\n=== load_plateau 読込結果 ===\n');
    fprintf('建物数         : %d 棟\n',       info.n_buildings);
    fprintf('座標原点       : lat=%.6f°, lon=%.6f°\n', info.origin_latlon(1), info.origin_latlon(2));
    fprintf('X範囲 (East)   : %7.1f〜%7.1f m  (%.0f m)\n', info.x_range(1),  info.x_range(2),  diff(info.x_range));
    fprintf('Y範囲 (North)  : %7.1f〜%7.1f m  (%.0f m)\n', info.y_range(1),  info.y_range(2),  diff(info.y_range));
    fprintf('建物高さ範囲   : %.1f〜%.1f m\n', info.height_range(1), info.height_range(2));

    heights = [buildings.height];
    heights = heights(~isnan(heights));
    fprintf('高さ 中央値    : %.1f m\n', median(heights));
    fprintf('高さ 平均      : %.1f m\n', mean(heights));
    fprintf('高さ分布       : 〜10m: %d棟 / 10〜30m: %d棟 / 30m+: %d棟\n', ...
        sum(heights < 10), sum(heights >= 10 & heights < 30), sum(heights >= 30));

    n_verts = cellfun(@(fp) size(fp,1), {buildings.footprint});
    fprintf('頂点数/棟      : min=%d / median=%d / max=%d\n', ...
        min(n_verts), round(median(n_verts)), max(n_verts));

    % ── フットプリント俯瞰プロット（最大5000棟）──────────────────
    N_plot = min(5000, info.n_buildings);
    fig = figure('Visible', 'off');
    hold on; axis equal; grid on; box on;
    for i = 1:N_plot
        fp = buildings(i).footprint;
        h  = buildings(i).height;
        if ~isnan(h)
            c = min(h / 40, 1);  % 高さで色付け (40m=最大)
        else
            c = 0.5;
        end
        fill(fp(:,1), fp(:,2), c, 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.3, 'FaceAlpha', 0.7);
    end
    colormap(parula); cb = colorbar; cb.Label.String = '正規化高さ (max 40m)';
    xlabel('East [m]'); ylabel('North [m]');
    title(sprintf('PLATEAU 建物フットプリント  N=%d棟 (最大%d棟表示)', info.n_buildings, N_plot));
    out_png = fullfile(this_dir, 'plateau_footprints.png');
    saveas(fig, out_png);
    fprintf('\nフットプリント図保存: citymodel/plateau_footprints.png\n');

    disp('load_plateau: OK');
    return;
end

%% ─── メイン読込処理 ─────────────────────────────────────────────

fprintf('CityGML読込中: %s\n', gml_file);
t0 = tic;

% ファイル全体を char 配列として読み込む（DOM より高速）
fid = fopen(gml_file, 'r', 'n', 'UTF-8');
if fid == -1; error('ファイルを開けません: %s', gml_file); end
raw = fread(fid, '*char')';
fclose(fid);
fprintf('  読込完了: %.1f MB, %.1fs\n', length(raw)/1e6, toc(t0));

% <core:cityObjectMember>...</core:cityObjectMember> でブロック分割
pat_s = '<core:cityObjectMember>';
pat_e = '</core:cityObjectMember>';
starts = strfind(raw, pat_s);
ends   = strfind(raw, pat_e);
n_blk  = min(length(starts), length(ends));
fprintf('  cityObjectMember 総数: %d\n', n_blk);

% 座標変換基準点（最初の有効建物の重心を原点に設定）
lat0 = NaN; lon0 = NaN;
R    = 6371000;  % 地球半径 [m]

buildings = struct('footprint',{},'height',{},'base_z',{},'center',{});
n_ok   = 0;
n_skip = 0;
t1 = tic;

for bi = 1:n_blk
    blk = raw(starts(bi) : ends(bi) + length(pat_e) - 1);

    % bldg:Building 以外（道路・植生・地形等）をスキップ
    if isempty(strfind(blk, '<bldg:Building'))
        n_skip = n_skip + 1;
        continue;
    end

    % ── 建物高さ ────────────────────────────────────────────────
    hm = regexp(blk, '<bldg:measuredHeight[^>]*>([0-9.eE+-]+)', 'tokens', 'once');
    height = tok2num(hm);
    if ~isnan(height) && height < 0; height = NaN; end  % -9999 は PLATEAU の「不明」マーカー

    % ── フットプリント座標 ──────────────────────────────────────
    % 優先1: lod0RoofEdge  (z=0 の底面2Dポリゴン、最も確実)
    % 優先2: lod1Solid の最初の posList (底面を含むいずれかの面)
    pos_str = first_poslist(blk, '<bldg:lod0RoofEdge>');
    if isempty(pos_str)
        pos_str = first_poslist(blk, '<bldg:lod1Solid>');
    end
    if isempty(pos_str)
        n_skip = n_skip + 1;
        continue;
    end

    % "lat lon z ..." → 数値配列
    vals = sscanf(strtrim(pos_str), '%f');
    if length(vals) < 9 || mod(length(vals), 3) ~= 0
        n_skip = n_skip + 1;
        continue;
    end
    n_pts  = length(vals) / 3;
    coords = reshape(vals, 3, n_pts)';   % [N×3]: lat lon elev

    % CityGML は始点＝終点の閉ループ → 最後の重複点を除去
    if size(coords,1) > 1 && norm(coords(1,:) - coords(end,:)) < 1e-10
        coords = coords(1:end-1, :);
    end

    lat    = coords(:,1);
    lon    = coords(:,2);
    base_z = coords(1,3);

    % 基準点設定（最初の有効建物）
    if isnan(lat0)
        lat0 = mean(lat);
        lon0 = mean(lon);
        fprintf('  座標原点設定: lat=%.6f°, lon=%.6f°\n', lat0, lon0);
    end

    % Equirectangular 投影: lat/lon → ローカルXY [m]
    x = (lon - lon0) .* cos(lat0 * pi/180) .* (R * pi/180);  % East  [m]
    y = (lat - lat0)                        .* (R * pi/180);  % North [m]

    n_ok = n_ok + 1;
    buildings(n_ok).footprint = [x, y];
    buildings(n_ok).height    = height;
    buildings(n_ok).base_z    = base_z;
    buildings(n_ok).center    = [mean(x), mean(y)];

    if mod(n_ok, 3000) == 0
        fprintf('  ... %d棟完了 (%.1fs)\n', n_ok, toc(t1));
    end
end

fprintf('抽出完了: %d棟 (スキップ %d, 計 %.1fs)\n', n_ok, n_skip, toc(t0));

% info 構造体
cx      = reshape([buildings.center], 2, [])';
h_all   = [buildings.height];
h_valid = h_all(~isnan(h_all));
info.n_buildings   = n_ok;
info.origin_latlon = [lat0, lon0];
info.x_range       = [min(cx(:,1)), max(cx(:,1))];
info.y_range       = [min(cx(:,2)), max(cx(:,2))];
info.height_range  = [min(h_valid),  max(h_valid)];

end  % function load_plateau


%% ─── ローカル関数 ────────────────────────────────────────────────

function pos_str = first_poslist(blk, section_tag)
% section_tag セクション内の最初の gml:posList テキストを返す
    pat = [section_tag, '.*?<gml:posList[^>]*>(.*?)</gml:posList>'];
    m   = regexp(blk, pat, 'tokens', 'once', 'dotall');
    if isempty(m); pos_str = ''; else; pos_str = m{1}; end
end

function v = tok2num(tok)
% regexp tokens から数値変換
    if isempty(tok); v = NaN; else; v = str2double(tok{1}); end
end
