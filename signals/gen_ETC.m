function [signal, params] = gen_ETC(varargin)
% ETC路側機信号（DSRC / ARIB STD-T75）のASKベースバンド信号を生成する
% 周波数: 5.8GHz帯, 占有帯域幅 4.4MHz, 変調: ASK (OOK), Unique Word付き

p = inputParser;
addParameter(p, 'fs',       44e6);   % サンプリング周波数 [Hz]  (帯域幅の10倍)
addParameter(p, 'bitrate',  1e6);    % ビットレート [bps]
addParameter(p, 'nbits',    256);    % ペイロードビット数
addParameter(p, 'modDepth', 1.0);    % 変調度 (1.0 = OOK: 0→0, 1→1)
addParameter(p, 'seed',     42);
parse(p, varargin{:});
opt = p.Results;

rng(opt.seed);

% Unique Word: ARIB STD-T75 フレーム同期用16ビットパターン
% (calc_preamble.m での相関検出と対応させること)
uw = logical([1 1 1 0 1 0 1 1 0 1 0 0 0 0 0 1]);

% ランダムペイロード
data = logical(randi([0 1], 1, opt.nbits));

% フレーム: [UW | payload]
bits = [uw, data];

% 1ビットあたりのサンプル数（矩形パルス成形）
sps = round(opt.fs / opt.bitrate);

% ASK振幅マッピング: 0 → (1 - modDepth), 1 → 1.0
amp = (1 - opt.modDepth) + opt.modDepth * double(bits);

% パルス成形（矩形）→ 列ベクトル
signal = repelem(amp, sps).';

params.fs         = opt.fs;
params.bitrate    = opt.bitrate;
params.uw         = uw;
params.sps        = sps;
params.nbits      = opt.nbits;
params.n_total    = length(bits);
params.signal_len = length(signal);
params.duration_s = length(signal) / opt.fs;
end
