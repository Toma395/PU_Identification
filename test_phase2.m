% test_phase2.m — フェーズ2動作確認: calc_kurtosis / calc_preamble
% SNRを振って特徴量の劣化傾向を確認（比較実験Bの予備検証）
addpath('signals');
addpath('features');

%% 信号生成
[etc, ep] = gen_ETC();
[uav, ~]  = gen_UAV();

snr_list = [20, 10, 5, 0, -5, -10];

fprintf('%-8s | %-16s %-8s | %-16s %-8s\n', ...
    'SNR[dB]', 'ETC kurt', 'ETC lbl', 'UAV kurt', 'UAV lbl');
fprintf('%s\n', repmat('-', 1, 62));

for snr = snr_list
    etc_n = awgn_add(etc, snr, true);   % 実数AWGN
    uav_n = awgn_add(uav, snr, false);  % 複素AWGN

    [k_etc, l_etc] = calc_kurtosis(etc_n);
    [k_uav, l_uav] = calc_kurtosis(uav_n);

    fprintf('%-8d | %-10.4f %-14s | %-10.4f %-8s\n', ...
        snr, k_etc, l_etc, k_uav, l_uav);
end

%% プリアンブル相関テスト
fprintf('\n--- プリアンブル相関 (ETC Unique Word) ---\n');
fprintf('%-8s | %-12s %-10s | %-12s %-10s\n', ...
    'SNR[dB]', 'ETC peak', 'detected', 'UAV peak', 'detected');
fprintf('%s\n', repmat('-', 1, 56));

for snr = snr_list
    etc_n = awgn_add(etc, snr, true);
    uav_n = awgn_add(uav, snr, false);

    [pk_etc, ~, det_etc] = calc_preamble(etc_n, ep.uw, ep.sps);
    [pk_uav, ~, det_uav] = calc_preamble(uav_n, ep.uw, ep.sps);

    fprintf('%-8d | %-12.4f %-10s | %-12.4f %-10s\n', ...
        snr, pk_etc, tf2str(det_etc), pk_uav, tf2str(det_uav));
end

%% 結合特徴量テスト（クリーン信号のみ）
fprintf('\n--- 結合特徴量（SNR=20dB）---\n');
etc_n = awgn_add(etc, 20, true);
uav_n = awgn_add(uav, 20, false);

[k_e, l_e] = calc_kurtosis(etc_n);
[p_e, ~, d_e] = calc_preamble(etc_n, ep.uw, ep.sps);
fprintf('ETC: kurt=%.4f(%s)  preamble_peak=%.4f(%s)  → %s\n', ...
    k_e, l_e, p_e, tf2str(d_e), combined_label(l_e, d_e));

[k_u, l_u] = calc_kurtosis(uav_n);
[p_u, ~, d_u] = calc_preamble(uav_n, ep.uw, ep.sps);
fprintf('UAV: kurt=%.4f(%s)  preamble_peak=%.4f(%s)  → %s\n', ...
    k_u, l_u, p_u, tf2str(d_u), combined_label(l_u, d_u));

disp('Phase 2 calc_kurtosis / calc_preamble: OK');


%% ローカル関数
function y = awgn_add(x, snr_db, is_real)
    pwr   = mean(abs(x).^2);
    npwr  = pwr / 10^(snr_db/10);
    if is_real
        y = x + sqrt(npwr) * randn(size(x));
    else
        y = x + sqrt(npwr/2) * (randn(size(x)) + 1j*randn(size(x)));
    end
end

function s = tf2str(b)
    if b; s = 'DETECT'; else; s = '---'; end
end

function lbl = combined_label(kurt_label, preamble_detected)
    % 尖度がETCかつUW検出 → ETC確定; どちらかでもUAV → UAV
    if strcmp(kurt_label, 'ETC') && preamble_detected
        lbl = 'ETC';
    else
        lbl = 'UAV';
    end
end
