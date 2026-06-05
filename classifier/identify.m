function [label, score, features] = identify(signal, uw, sps, method)
% ETC/UAV識別器: 尖度・プリアンブル相関・結合の3手法に対応
%
% 入力:
%   signal  - 受信信号（実数 or 複素数）
%   uw      - Unique Wordビット列（gen_ETC params.uw と共通）
%   sps     - 1ビットあたりサンプル数（gen_ETC params.sps）
%   method  - 'kurtosis' / 'preamble' / 'combined'（デフォルト: 'combined'）
% 出力:
%   label    - 'ETC' or 'UAV'
%   score    - 識別スコア（高いほどETC候補）
%   features - 構造体: .kurt, .preamble_peak

if nargin < 4
    method = 'combined';
end

kurt          = calc_kurtosis(signal);
preamble_peak = calc_preamble(signal, uw, sps);

% [0,1]正規化: ETC→1, UAV→0
k_norm = max(0, min(1, (3 - kurt) / 2));

features.kurt          = kurt;
features.preamble_peak = preamble_peak;

switch method
    case 'kurtosis'
        score = 3 - kurt;           % ETC(clean)≈2, UAV≈0
        label = pick(score >= 1.0, 'ETC', 'UAV');
    case 'preamble'
        score = preamble_peak;      % ETC≈1.0, UAV≈0.1
        label = pick(score >= 0.5, 'ETC', 'UAV');
    case 'combined'
        % 低SNRではプリアンブル相関が支配、高SNRでは尖度が補完
        score = 0.5 * k_norm + 0.5 * preamble_peak;
        label = pick(score >= 0.5, 'ETC', 'UAV');
    otherwise
        error('method は kurtosis / preamble / combined のいずれかを指定');
end
end

function s = pick(cond, a, b)
if cond; s = a; else; s = b; end
end
