% experiments/compare_multipath.m
% 比較実験: マルチパス環境でのETC検出性能
%
% 目的: タップ遅延線チャネルでプリアンブル相関のピークがにじむとき
%       尖度・プリアンブル・結合の3手法がどこまで耐えるかを検証する。
%
% マルチパスモデル:
%   受信 = sum_{l=0}^{L-1} h_l * s(t - tau_l) + AWGN
%   h_l  = 複素ガウスタップ（振幅がレイリー分布）
%   power_l = exp(-l/decay)、全パス正規化（受信電力一定）
%   tau_0 = 0（直接波）、tau_l>0 はランダム [1, tau_max]
%   tau_max = 5*sps = 220 サンプル（≈5 μs @ 44MHz = 5ビット幅）
%
% 評価軸:
%   主実験  : パス数 L を振る（L=1,2,3,5,10）
%   副実験  : マルチパス + 干渉UAV（SIR=0dB）の複合条件
%
% 出力:
%   experiments/mp_accuracy.png   -- パス数 vs 精度（3手法）
%   experiments/mp_corr_shape.png -- 相関波形例（L=1 vs L=10）
%   experiments/mp_sir_combined.png -- 複合条件（マルチパス×SIR）

base_dir = fileparts(which('compare_multipath'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));

%% ─────────────────────────────────────────────────────────────────
%% パラメータ
%% ─────────────────────────────────────────────────────────────────
L_list   = [1, 2, 3, 5, 10];   % パス数
N_L      = length(L_list);
N_trials = 300;   % 試行数/クラス
snr_db   = 20;    % AWGN SNR [dB]（干渉なし主実験）
thresh   = [1.0, 0.5, 0.5];    % kurtosis / preamble / combined

% タップ遅延線パラメータ
[~, ep]  = gen_ETC('seed', 1);
uw  = ep.uw;
sps = ep.sps;       % 44
tau_max = 5 * sps;  % 220 サンプル ≈ 5 μs（5ビット幅）

% UAV混在条件用
RS_P = 11; RS_Q = 5;  % 20MHz → 44MHz リサンプル比
sir_list_cmb = [10, 5, 0, -5];   % dB

mnames = {'Kurtosis','Preamble','Combined'};

%% ─────────────────────────────────────────────────────────────────
%% 主実験: パス数 L vs 検出精度
%% ─────────────────────────────────────────────────────────────────
fprintf('=== 主実験: パス数 vs 検出精度 (%d L × %d試行/クラス) ===\n', N_L, N_trials);

scores_mp = zeros(N_trials*2, 3, N_L);
labels_mp = [ones(N_trials,1); zeros(N_trials,1)];
acc_mp    = zeros(N_L, 3);
f1_mp     = zeros(N_L, 3);

% 相関波形サンプル（L=1とL=10の比較用）
corr_example = cell(2,1);   % {L=1, L=10}

for li = 1:N_L
    L_paths = L_list(li);

    for ti = 1:N_trials
        seed = ti*37 + li*200;

        % ── 正例: ETC over multipath ──────────────────────────────
        [s_etc, ~] = gen_ETC('seed', seed);
        y_etc_mp = apply_multipath(s_etc, L_paths, tau_max, seed);
        y_pos = add_awgn(y_etc_mp, snr_db);
        k = calc_kurtosis(y_pos);
        p = calc_preamble(y_pos, uw, sps);
        scores_mp(ti,:,li) = make_scores(k,p);

        % ── 負例: UAV（マルチパスなし）──────────────────────────
        [s_uav, ~] = gen_UAV('seed', seed + 1000);
        y_neg = add_awgn(s_uav, snr_db);
        k = calc_kurtosis(y_neg);
        p = calc_preamble(y_neg, uw, sps);
        scores_mp(N_trials+ti,:,li) = make_scores(k,p);

        % 相関波形サンプル保存（各Lの最初の1試行）
        if ti==1 && (L_paths==1 || L_paths==10)
            store_idx = 1 + (L_paths==10);
            [~, ~, ~, corr_full] = calc_preamble_debug(y_pos, uw, sps);
            corr_example{store_idx} = corr_full;
        end
    end

    for mi = 1:3
        [acc_mp(li,mi), f1_mp(li,mi)] = calc_acc_f1(labels_mp, scores_mp(:,mi,li) >= thresh(mi));
    end
    fprintf('  L=%2d  完了\n', L_paths);
end

%% ─────────────────────────────────────────────────────────────────
%% 副実験: マルチパス（L=5） + 干渉UAV  複合条件
%% ─────────────────────────────────────────────────────────────────
L_cmb    = 5;
N_sir    = length(sir_list_cmb);
fprintf('\n=== 副実験: マルチパス(L=%d) + 干渉UAV  SIR sweeps ===\n', L_cmb);

acc_cmb  = zeros(N_sir, 3);

for si = 1:N_sir
    sir_db_val = sir_list_cmb(si);
    sir_lin    = 10^(sir_db_val/10);

    scores_cmb = zeros(N_trials*2, 3);
    labels_cmb = [ones(N_trials,1); zeros(N_trials,1)];

    for ti = 1:N_trials
        seed = ti*41 + si*300 + 9000;

        % ── 正例: ETC(マルチパス) + UAV干渉 ─────────────────────
        [s_etc, ~] = gen_ETC('seed', seed);
        y_etc_mp   = apply_multipath(s_etc, L_cmb, tau_max, seed);
        [s_uav, ~] = gen_UAV('seed', seed + 500);
        y_pos      = mix_etc_uav(y_etc_mp, s_uav, sir_lin, snr_db, RS_P, RS_Q);
        k = calc_kurtosis(y_pos);
        p = calc_preamble(y_pos, uw, sps);
        scores_cmb(ti,:) = make_scores(k,p);

        % ── 負例: UAV干渉のみ ────────────────────────────────────
        [s_uav2, ~] = gen_UAV('seed', seed + 2000);
        y_neg = uav_only(s_uav2, length(s_etc), snr_db, RS_P, RS_Q);
        k = calc_kurtosis(y_neg);
        p = calc_preamble(y_neg, uw, sps);
        scores_cmb(N_trials+ti,:) = make_scores(k,p);
    end

    for mi = 1:3
        [acc_cmb(si,mi), ~] = calc_acc_f1(labels_cmb, scores_cmb(:,mi) >= thresh(mi));
    end
    fprintf('  SIR=%4ddB  完了\n', sir_db_val);
end

%% ─────────────────────────────────────────────────────────────────
%% 統計サマリ
%% ─────────────────────────────────────────────────────────────────
fprintf('\n=== 主実験: パス数 vs 検出精度 ===\n');
fprintf('%-6s | %-12s %-12s %-12s\n', 'L_paths', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for li = 1:N_L
    fprintf('%-6d | %-12.3f %-12.3f %-12.3f\n', L_list(li), acc_mp(li,:));
end

fprintf('\n=== 副実験: マルチパス(L=%d) + UAV干渉 ===\n', L_cmb);
fprintf('%-8s | %-12s %-12s %-12s\n', 'SIR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_sir
    fprintf('%-8d | %-12.3f %-12.3f %-12.3f\n', sir_list_cmb(si), acc_cmb(si,:));
end

fprintf('\n=== 90%%精度維持（主実験: SNR=20dB） ===\n');
for mi = 1:3
    idx = find(acc_mp(:,mi) >= 0.90, 1, 'last');
    if isempty(idx)
        fprintf('  %-10s: 全条件で90%%未達\n', mnames{mi});
    else
        fprintf('  %-10s: L = %d まで 90%%維持\n', mnames{mi}, L_list(idx));
    end
end

%% ─────────────────────────────────────────────────────────────────
%% プロット
%% ─────────────────────────────────────────────────────────────────
lstyles = {'-','--','-.'};
colors  = lines(3);

%% ① パス数 vs 精度・F1（2段）
fig1 = figure('Visible','off','Position',[0 0 900 700]);

subplot(2,1,1);
hold on; grid on; box on;
for mi = 1:3
    plot(L_list, acc_mp(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',1.8,'Marker','o','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','left','FontSize',8);
yline(0.5,'k:','50%','LabelHorizontalAlignment','left','FontSize',8);
xlabel('マルチパス数 L（1=直接波のみ）');
ylabel('Accuracy');
title(sprintf('マルチパス vs ETC検出精度  (SNR=%ddB, τ_{max}=%dサンプル=%.1fμs)', ...
    snr_db, tau_max, tau_max/44e6*1e6));
legend('Location','southwest'); ylim([0.4 1.05]);
set(gca,'XTick',L_list);

subplot(2,1,2);
hold on; grid on; box on;
for mi = 1:3
    plot(L_list, f1_mp(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',1.8,'Marker','s','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','F1=0.9','LabelHorizontalAlignment','left','FontSize',8);
xlabel('マルチパス数 L');
ylabel('F1 Score');
title('マルチパス vs F1スコア');
legend('Location','southwest'); ylim([0.4 1.05]);
set(gca,'XTick',L_list);

saveas(fig1, fullfile(base_dir,'mp_accuracy.png'));
fprintf('\n保存: experiments/mp_accuracy.png\n');

%% ② 相関波形の可視化（L=1 vs L=10）
fig2 = figure('Visible','off','Position',[0 0 1100 500]);
L_labels  = {sprintf('L=1（直接波のみ、τ_{max}=%dサンプル）', tau_max), ...
             sprintf('L=10（10パス、τ_{max}=%dサンプル）', tau_max)};
plot_titles = {'L=1: 直接波のみ', 'L=10: 10マルチパス'};

for pi = 1:2
    subplot(1,2,pi);
    hold on; grid on; box on;

    cv = corr_example{pi};
    if isempty(cv)
        text(0.5,0.5,'データなし','Units','normalized','HorizontalAlignment','center');
        continue;
    end

    n_show = min(length(cv), 4*tau_max + size(repelem(double(uw(:)),sps),1));
    t_ax   = (0:n_show-1);
    plot(t_ax, abs(cv(1:n_show)), 'b-', 'LineWidth',0.8, 'DisplayName','|相関|');

    % UW長（理想ピーク位置）をマーク
    uw_len = sps * length(uw);
    xline(uw_len, 'r--', 'UW長', 'LabelHorizontalAlignment','left', 'FontSize',8, 'LineWidth',1.2);

    [pv, pi_idx] = max(abs(cv(1:n_show)));
    plot(pi_idx-1, pv, 'r*', 'MarkerSize',12, 'DisplayName', sprintf('ピーク=%.3f', pv));

    xlabel('ラグ [サンプル]');
    ylabel('正規化相関');
    title(plot_titles{pi});
    legend('Location','northeast','FontSize',8);
    ylim([0, max(abs(cv(1:n_show)))*1.15 + eps]);
end
sgtitle('プリアンブル相関波形  L=1 vs L=10  (ETC + AWGN SNR=20dB, seed=37)', ...
    'FontSize',12);
saveas(fig2, fullfile(base_dir,'mp_corr_shape.png'));
fprintf('保存: experiments/mp_corr_shape.png\n');

%% ③ 複合条件（マルチパス＋UAV干渉）
fig3 = figure('Visible','off','Position',[0 0 900 480]);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list_cmb, acc_cmb(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',1.8,'Marker','o','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','left','FontSize',8);
xlabel('SIR [dB]  (ETC電力 / UAV干渉電力)');
ylabel('Accuracy');
title(sprintf('複合条件: マルチパス(L=%d) + UAV干渉  (SNR=%ddB)', L_cmb, snr_db));
legend('Location','southwest'); ylim([0.4 1.05]);
saveas(fig3, fullfile(base_dir,'mp_sir_combined.png'));
fprintf('保存: experiments/mp_sir_combined.png\n');

disp('compare_multipath: OK');


%% ─────────────────────────────────────────────────────────────────
%% ローカル関数
%% ─────────────────────────────────────────────────────────────────

function y = apply_multipath(s, L_paths, tau_max, seed)
% タップ遅延線マルチパスチャネルを適用する
% 直接波（l=0）: tau=0, 固定最大電力
% 反射波（l>=1）: tau=random[1,tau_max], 複素ガウスタップ（レイリー振幅）
% 電力プロファイル: 指数減衰 exp(-l/decay), 全パス正規化
    rng(seed);
    N      = length(s);
    N_out  = N + tau_max;
    buf    = zeros(N_out, 1);

    % 指数減衰電力プロファイル（正規化）
    decay  = L_paths / 3;
    powers = exp(-(0:L_paths-1) / decay);
    powers = powers / sum(powers);

    for l = 1:L_paths
        if l == 1
            tau = 0;                              % 直接波
            h   = sqrt(powers(l));                % 実数・最大電力
        else
            tau = randi([1, tau_max]);            % ランダム遅延
            h   = (randn + 1j*randn)/sqrt(2) * sqrt(powers(l));  % 複素ガウス
        end
        buf(tau+1 : tau+N) = buf(tau+1 : tau+N) + h * s;
    end

    y = buf(1:N);   % 元の長さに切り詰め
end

function y = add_awgn(x, snr_db)
    p_sig   = mean(abs(x).^2);
    p_noise = p_sig / 10^(snr_db/10);
    if isreal(x)
        y = x + sqrt(p_noise)*randn(size(x));
    else
        y = x + sqrt(p_noise/2)*(randn(size(x)) + 1j*randn(size(x)));
    end
end

function y = mix_etc_uav(s_etc, s_uav, sir_lin, snr_db, rs_p, rs_q)
% s_etc はすでにマルチパス後の複素信号でもよい
    s_uav_rs = resample(s_uav, rs_p, rs_q);
    s_etc_n  = s_etc / sqrt(mean(abs(s_etc).^2));
    s_uav_n  = s_uav_rs / sqrt(mean(abs(s_uav_rs).^2));
    L = min(length(s_etc_n), length(s_uav_n));
    mixed = sqrt(sir_lin)*s_etc_n(1:L) + s_uav_n(1:L);
    y = add_awgn(mixed, snr_db);
end

function y = uav_only(s_uav, L_ref, snr_db, rs_p, rs_q)
    s_uav_rs = resample(s_uav, rs_p, rs_q);
    s_uav_n  = s_uav_rs / sqrt(mean(abs(s_uav_rs).^2));
    L = min(length(s_uav_n), L_ref);
    y = add_awgn(s_uav_n(1:L), snr_db);
end

function sc = make_scores(kurt, p_peak)
    k_norm = max(0, min(1, (3-kurt)/2));
    sc = [3-kurt,  p_peak,  0.5*k_norm + 0.5*p_peak];
end

function [acc, f1] = calc_acc_f1(y_true, y_pred)
    tp  = sum(y_true==1 & y_pred==1);
    tn  = sum(y_true==0 & y_pred==0);
    fp  = sum(y_true==0 & y_pred==1);
    fn  = sum(y_true==1 & y_pred==0);
    acc = (tp+tn)/(tp+tn+fp+fn);
    pre = tp/max(tp+fp,1);
    rec = tp/max(tp+fn,1);
    f1  = 2*pre*rec/max(pre+rec,eps);
end

function [peak_corr, peak_pos, detected, corr_abs] = calc_preamble_debug(signal, uw, sps, threshold)
% calc_preamble の内部計算を返す（相関波形可視化用）
    if nargin < 4; threshold = 0.5; end
    ref        = repelem(double(uw(:)), sps);
    sig        = real(signal(:));
    corr_full  = xcorr(sig, ref);
    ref_energy = sum(ref.^2);
    corr_full  = corr_full / (ref_energy + eps);
    zero_lag   = length(sig);
    corr_pos   = corr_full(zero_lag:end);
    [peak_corr, idx] = max(abs(corr_pos));
    peak_pos   = idx - 1;
    detected   = peak_corr >= threshold;
    corr_abs   = abs(corr_pos);
end
