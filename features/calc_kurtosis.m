function [kurt, label] = calc_kurtosis(signal, threshold)
% 受信信号の尖度を計算し ETC / UAV ラベルを返す
% 入力:
%   signal    - 受信信号（実数 or 複素数）
%   threshold - 判別閾値（デフォルト 2.0）
%               ETC(OOK) ≈ 1 / Gaussian ≈ 3 の中点
% 出力:
%   kurt  - 尖度値（MATLAB kurtosis: Gaussian=3, OOK=1）
%   label - 'ETC' or 'UAV'

if nargin < 2
    threshold = 2.0;
end

% 複素信号は実部を使う（OFDM実部はCLTによりGaussian）
if isreal(signal)
    x = signal(:);
else
    x = real(signal(:));
end

kurt = kurtosis(x);
label = choose_label(kurt, threshold);
end

function label = choose_label(kurt, threshold)
if kurt < threshold
    label = 'ETC';
else
    label = 'UAV';
end
end
