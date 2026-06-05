function [signal, params] = gen_UAV(varargin)
% UAV通信信号（IEEE 802.11a系 OFDM）のベースバンド複素信号を生成する
% FFT64 / CP16 / 52サブキャリア / サブキャリア変調: QPSK または 16QAM
% 注: OFDMは多重化方式, QPSK/16QAMはサブキャリア変調（階層関係）

p = inputParser;
addParameter(p, 'fs',       20e6);   % サンプリング周波数 [Hz]
addParameter(p, 'Nfft',     64);     % FFTサイズ
addParameter(p, 'Ncp',      16);     % サイクリックプレフィックス長
addParameter(p, 'Nsym',     40);     % OFDMシンボル数
addParameter(p, 'modOrder', 4);      % サブキャリア変調多値数 (4=QPSK, 16=16QAM)
addParameter(p, 'seed',     42);
parse(p, varargin{:});
opt = p.Results;

rng(opt.seed);

Nfft = opt.Nfft;
Ncp  = opt.Ncp;
Nsym = opt.Nsym;

% 802.11a サブキャリア割り当て: ±1〜±26 (DCゼロ, 計52本)
% MATLAB 1-based・centered表現: DC = Nfft/2+1 = 33
center = Nfft/2 + 1;
sc_idx = [(center-26):(center-1), (center+1):(center+26)];  % 52本
Nsc    = length(sc_idx);

% QAMシンボル生成
switch opt.modOrder
    case 4   % QPSK: 単位円上に正規化
        re  = 2*randi([0 1], Nsc*Nsym, 1) - 1;
        im  = 2*randi([0 1], Nsc*Nsym, 1) - 1;
        qam = (re + 1j*im) / sqrt(2);
    case 16  % 16QAM: 平均電力1に正規化
        alph = [-3 -1 1 3] / sqrt(10);
        re   = alph(randi([1 4], Nsc*Nsym, 1));
        im   = alph(randi([1 4], Nsc*Nsym, 1));
        qam  = re(:) + 1j*im(:);
    otherwise
        error('modOrder は 4 または 16 を指定');
end

% OFDMシンボル生成
sym_len = Nfft + Ncp;
signal  = zeros(Nsym * sym_len, 1);

for k = 1:Nsym
    % 周波数領域: centered形式（DC at center）でサブキャリア割り当て
    fd = zeros(Nfft, 1);
    fd(sc_idx) = qam((k-1)*Nsc+1 : k*Nsc);

    % ifftshift: centered→標準順序 (DC→index 1) → IFFT → 時間領域
    td = ifft(ifftshift(fd)) * sqrt(Nfft);

    % サイクリックプレフィックス付加
    cp = td(end-Ncp+1:end);
    signal((k-1)*sym_len+1 : k*sym_len) = [cp; td];
end

params.fs         = opt.fs;
params.Nfft       = Nfft;
params.Ncp        = Ncp;
params.Nsym       = Nsym;
params.Nsc        = Nsc;
params.modOrder   = opt.modOrder;
params.signal_len = length(signal);
params.duration_s = length(signal) / opt.fs;
end
