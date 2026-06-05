function [peak_corr, peak_pos, detected] = calc_preamble(signal, uw, sps, threshold)
% ETC Unique Word との正規化相互相関を計算しプリアンブル検出を行う
%
% 入力:
%   signal    - 受信信号（実数 or 複素数）
%   uw        - Unique Wordビット列（論理 or 0/1ベクトル、gen_ETC と共通）
%   sps       - 1ビットあたりのサンプル数（gen_ETC の params.sps）
%   threshold - 正規化ピーク検出閾値（デフォルト: 0.6）
% 出力:
%   peak_corr - 正規化相関ピーク値 [0, 1]
%   peak_pos  - ピークの遅延サンプル数（0始まり）
%   detected  - true: UW検出（ETC候補）/ false: 非ETC
%
% 備考:
%   マルチパス環境でも相関ピークは残りやすい（CLAUDE.md 識別手法参照）

if nargin < 4
    threshold = 0.5;  % 参照エネルギー正規化: ETC≈1.0, UAV≪0.5
end

% OOK参照波形: gen_ETC と同じ 0/1 矩形パルス
ref = repelem(double(uw(:)), sps);   % 長さ length(uw)*sps

sig = real(signal(:));               % 実部を使用

% 参照エネルギー正規化: perfect match → 1.0, 非マッチ → 0 付近
% Pearson正規化（sqrt(||sig||*||ref||)）は sig が長い場合に過剰補正になるため不使用
corr_full  = xcorr(sig, ref);
ref_energy = sum(ref.^2);
corr_full  = corr_full / (ref_energy + eps);
% xcorr(x,y) のゼロラグインデックス = length(x)
zero_lag  = length(sig);
corr_pos  = corr_full(zero_lag:end); % lag=0 以降（ref が sig の先頭に整合する方向）

[peak_corr, idx] = max(abs(corr_pos));
peak_pos = idx - 1;                  % 0-based 遅延サンプル数
detected = peak_corr >= threshold;
end
