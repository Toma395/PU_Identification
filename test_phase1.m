% test_phase1.m — フェーズ1動作確認: gen_ETC / gen_UAV
addpath('signals');

%% ETC (ASK) 信号生成
[etc, ep] = gen_ETC();
k_etc = kurtosis(etc);
fprintf('--- ETC (ASK / ARIB STD-T75) ---\n');
fprintf('  サンプル数  : %d\n',      ep.signal_len);
fprintf('  継続時間    : %.2f us\n', ep.duration_s * 1e6);
fprintf('  平均電力    : %.4f\n',    mean(etc.^2));
fprintf('  尖度        : %.4f  (期待値 ≈ 1, 二値離散)\n', k_etc);

%% UAV (OFDM) 信号生成
[uav, up] = gen_UAV();
k_uav = kurtosis(real(uav));
fprintf('--- UAV (OFDM / 802.11a系) ---\n');
fprintf('  サンプル数  : %d\n',      up.signal_len);
fprintf('  継続時間    : %.2f us\n', up.duration_s * 1e6);
fprintf('  平均電力    : %.4f\n',    mean(abs(uav).^2));
fprintf('  尖度(実部)  : %.4f  (期待値 ≈ 3, ガウス的)\n', k_uav);

%% 尖度マージン確認
margin = abs(k_etc - k_uav);
fprintf('\n--- 尖度による判別マージン ---\n');
fprintf('  ETC 尖度 = %.4f,  UAV 尖度 = %.4f,  差 = %.4f\n', k_etc, k_uav, margin);
assert(margin > 1.0, 'ERROR: 判別マージンが 1.0 未満 — 識別が困難');
fprintf('  判別マージン %.2f > 1.0 → PASS\n', margin);

disp('Phase 1 gen_ETC / gen_UAV: OK');
