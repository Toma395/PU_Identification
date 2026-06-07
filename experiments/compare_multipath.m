% experiments/compare_multipath.m
% 比較実験: マルチパス環境でのETC検出性能
%
% 改訂: L=1〜10 全整数, チャネル実現数200×信号試行3回=600試行/クラス,
%       副実験SIR -10〜+14dB 2dB刻み, apply_paper_style 適用

base_dir = fileparts(which('compare_multipath'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..'));          % apply_paper_style

%% パラメータ
L_list   = 1:10;    % パス数 全整数（10値）
N_L      = length(L_list);
N_ch     = 200;     % チャネル実現数（チャネル平均用）
N_sig    = 3;       % 信号試行数/チャネル
snr_db   = 20;      % AWGN SNR [dB]
thresh   = [1.0, 0.5, 0.5];

[~, ep]  = gen_ETC('seed', 1);
uw  = ep.uw;
sps = ep.sps;
tau_max = 5 * sps;  % 220 サンプル ≈ 5 μs

RS_P = 11; RS_Q = 5;
sir_list_cmb = -10:2:14;   % dB (13点, 昇順)

mnames = {'Kurtosis','Preamble','Combined'};

%% 主実験: パス数 L vs 検出精度（N_ch × N_sig で平均）
fprintf('=== 主実験: L=1〜%d  N_ch=%d × N_sig=%d =%d試行/クラス ===\n', ...
    L_list(end), N_ch, N_sig, N_ch*N_sig);

acc_mp    = zeros(N_L, 3);
f1_mp     = zeros(N_L, 3);

corr_example = cell(2,1);   % {L=1, L=L_list(end)=10}

for li = 1:N_L
    L_paths = L_list(li);
    n_tot   = N_ch * N_sig;

    scores_pos = zeros(n_tot, 3);
    scores_neg = zeros(n_tot, 3);
    row = 0;

    for ch_i = 1:N_ch
        ch_seed = ch_i * 97 + li * 10007;   % チャネルシード（信号シードと独立）

        for sig_i = 1:N_sig
            row = row + 1;
            sig_seed_etc = ch_i * 1009 + sig_i * 37 + li * 100003;
            sig_seed_uav = sig_seed_etc + 9999983;

            % 正例: ETC over multipath
            [s_etc, ~] = gen_ETC('seed', sig_seed_etc);
            y_mp  = apply_multipath(s_etc, L_paths, tau_max, ch_seed);
            y_pos = add_awgn(y_mp, snr_db);
            k = calc_kurtosis(y_pos);
            p = calc_preamble(y_pos, uw, sps);
            scores_pos(row,:) = make_scores(k,p);

            % 負例: UAV（マルチパスなし）
            [s_uav, ~] = gen_UAV('seed', sig_seed_uav);
            y_neg = add_awgn(s_uav, snr_db);
            k = calc_kurtosis(y_neg);
            p = calc_preamble(y_neg, uw, sps);
            scores_neg(row,:) = make_scores(k,p);

            % 相関波形サンプル保存（ch_i=1, sig_i=1 の最初の1回）
            if ch_i==1 && sig_i==1
                if L_paths == L_list(1)
                    [~,~,~,corr_full] = calc_preamble_debug(y_pos, uw, sps);
                    corr_example{1} = corr_full;
                elseif L_paths == L_list(end)
                    [~,~,~,corr_full] = calc_preamble_debug(y_pos, uw, sps);
                    corr_example{2} = corr_full;
                end
            end
        end
    end

    all_scores = [scores_pos; scores_neg];
    all_labels = [ones(n_tot,1); zeros(n_tot,1)];
    for mi = 1:3
        [acc_mp(li,mi), f1_mp(li,mi)] = calc_acc_f1(all_labels, all_scores(:,mi) >= thresh(mi));
    end
    fprintf('  L=%2d  (ch=%d × sig=%d = %d試行/クラス)  完了\n', L_paths, N_ch, N_sig, n_tot);
end

%% 副実験: マルチパス（L=5） + 干渉UAV  複合条件
L_cmb = 5;
N_sir_cmb = length(sir_list_cmb);
n_tot_cmb = N_ch * N_sig;
fprintf('\n=== 副実験: マルチパス(L=%d) + 干渉UAV  SIR sweep ===\n', L_cmb);

acc_cmb = zeros(N_sir_cmb, 3);

for si = 1:N_sir_cmb
    sir_db_val = sir_list_cmb(si);
    sir_lin    = 10^(sir_db_val/10);
    scores_pos_c = zeros(n_tot_cmb, 3);
    scores_neg_c = zeros(n_tot_cmb, 3);
    row = 0;

    for ch_i = 1:N_ch
        ch_seed = ch_i * 97 + si * 20007;

        for sig_i = 1:N_sig
            row = row + 1;
            seed = ch_i * 1009 + sig_i * 41 + si * 200003;

            [s_etc, ~] = gen_ETC('seed', seed);
            y_mp       = apply_multipath(s_etc, L_cmb, tau_max, ch_seed);
            [s_uav, ~] = gen_UAV('seed', seed + 500013);
            y_pos      = mix_etc_uav(y_mp, s_uav, sir_lin, snr_db, RS_P, RS_Q);
            k = calc_kurtosis(y_pos);
            p = calc_preamble(y_pos, uw, sps);
            scores_pos_c(row,:) = make_scores(k,p);

            [s_uav2, ~] = gen_UAV('seed', seed + 2000017);
            y_neg = uav_only(s_uav2, length(s_etc), snr_db, RS_P, RS_Q);
            k = calc_kurtosis(y_neg);
            p = calc_preamble(y_neg, uw, sps);
            scores_neg_c(row,:) = make_scores(k,p);
        end
    end

    all_sc = [scores_pos_c; scores_neg_c];
    all_lb = [ones(n_tot_cmb,1); zeros(n_tot_cmb,1)];
    for mi = 1:3
        [acc_cmb(si,mi), ~] = calc_acc_f1(all_lb, all_sc(:,mi) >= thresh(mi));
    end
    fprintf('  SIR=%+4ddB  完了\n', sir_db_val);
end

%% サマリ
fprintf('\n=== 主実験: パス数 vs 検出精度 ===\n');
fprintf('%-6s | %-12s %-12s %-12s\n', 'L_paths', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for li = 1:N_L
    fprintf('%-6d | %-12.3f %-12.3f %-12.3f\n', L_list(li), acc_mp(li,:));
end

fprintf('\n=== 副実験: マルチパス(L=%d) + UAV干渉 ===\n', L_cmb);
fprintf('%-8s | %-12s %-12s %-12s\n', 'SIR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_sir_cmb
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

% 重要確認: L=5 + SIR=0dB の複合精度
si0 = find(sir_list_cmb == 0, 1);
if ~isempty(si0)
    fprintf('\n=== 複合条件確認 (L=%d, SIR=0dB) ===\n', L_cmb);
    fprintf('  %-10s: %.3f\n', mnames{1}, acc_cmb(si0,1));
    fprintf('  %-10s: %.3f\n', mnames{2}, acc_cmb(si0,2));
    fprintf('  %-10s: %.3f\n', mnames{3}, acc_cmb(si0,3));
end

%% プロット
lstyles = {'-','--','-.'};
colors  = lines(3);

%% ① パス数 vs 精度・F1（2段）
fig1 = figure('Visible','off','Color','w','Position',[0 0 900 700]);

subplot(2,1,1);
hold on; grid on; box on;
for mi = 1:3
    plot(L_list, acc_mp(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','left','FontSize',12,'HandleVisibility','off');
yline(0.5,'k:','HandleVisibility','off');
xlabel('マルチパス数 L（1=直接波のみ）','FontSize',14);
ylabel('Accuracy','FontSize',14);
title(sprintf('マルチパス vs ETC検出精度  (SNR=%ddB, N_{ch}=%d×N_{sig}=%d)', ...
    snr_db, N_ch, N_sig),'FontSize',16);
legend('Location','southwest','FontSize',12); ylim([0.4 1.05]);
xlim([0.5, L_list(end)+0.5]);
set(gca,'XTick',L_list);

subplot(2,1,2);
hold on; grid on; box on;
for mi = 1:3
    plot(L_list, f1_mp(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','s','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','F1=0.9','LabelHorizontalAlignment','left','FontSize',12,'HandleVisibility','off');
xlabel('マルチパス数 L','FontSize',14);
ylabel('F1 Score','FontSize',14);
title('マルチパス vs F1スコア','FontSize',16);
legend('Location','southwest','FontSize',12); ylim([0.4 1.05]);
xlim([0.5, L_list(end)+0.5]);
set(gca,'XTick',L_list);

apply_paper_style(fig1);
saveas(fig1, fullfile(base_dir,'mp_accuracy.png'));
fprintf('\n保存: experiments/mp_accuracy.png\n');

%% ② 相関波形（L=1 vs L=10, ch_i=1の代表1チャネル）
fig2 = figure('Visible','off','Color','w','Position',[0 0 1100 500]);
plot_titles = {sprintf('L=%d: 直接波のみ', L_list(1)), ...
               sprintf('L=%d: %dマルチパス', L_list(end), L_list(end))};

for pi = 1:2
    subplot(1,2,pi);
    hold on; grid on; box on;

    cv = corr_example{pi};
    if isempty(cv)
        text(0.5,0.5,'データなし','Units','normalized','HorizontalAlignment','center');
        continue;
    end

    n_show = min(length(cv), 4*tau_max + sps*length(uw));
    t_ax   = 0:n_show-1;
    plot(t_ax, abs(cv(1:n_show)), 'b-', 'LineWidth',2.0, 'DisplayName','|相関|');

    uw_len = sps * length(uw);
    xline(uw_len,'r--','UW長','LabelHorizontalAlignment','left','FontSize',12,'LineWidth',1.5,...
        'HandleVisibility','off');
    [pv, pi_idx] = max(abs(cv(1:n_show)));
    plot(pi_idx-1, pv, 'r*', 'MarkerSize',14,'LineWidth',1.5,...
        'DisplayName',sprintf('ピーク=%.3f', pv));

    xlabel('ラグ [サンプル]','FontSize',14);
    ylabel('正規化相関','FontSize',14);
    title(plot_titles{pi},'FontSize',16);
    legend('Location','northeast','FontSize',12);
    ylim([0, max(abs(cv(1:n_show)))*1.15 + eps]);
end
sgtitle(sprintf('プリアンブル相関波形  L=%d vs L=%d  (ETC+AWGN SNR=%ddB)', ...
    L_list(1), L_list(end), snr_db), 'FontSize',15,'Color',[0.15 0.15 0.15]);
apply_paper_style(fig2);
saveas(fig2, fullfile(base_dir,'mp_corr_shape.png'));
fprintf('保存: experiments/mp_corr_shape.png\n');

%% ③ 複合条件（マルチパス＋UAV干渉）
fig3 = figure('Visible','off','Color','w','Position',[0 0 900 480]);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list_cmb, acc_cmb(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
xlabel('SIR [dB]  (ETC電力 / UAV干渉電力)','FontSize',14);
ylabel('Accuracy','FontSize',14);
title(sprintf('複合条件: マルチパス(L=%d) + UAV干渉  (SNR=%ddB)', L_cmb, snr_db),'FontSize',16);
legend('Location','southwest','FontSize',12); ylim([0.4 1.05]);
xlim([sir_list_cmb(1)-1, sir_list_cmb(end)+1]);
apply_paper_style(fig3);
saveas(fig3, fullfile(base_dir,'mp_sir_combined.png'));
fprintf('保存: experiments/mp_sir_combined.png\n');

disp('compare_multipath: OK');


%% ローカル関数
function y = apply_multipath(s, L_paths, tau_max, seed)
    rng(seed);
    N      = length(s);
    N_out  = N + tau_max;
    buf    = zeros(N_out, 1);
    decay  = L_paths / 3;
    powers = exp(-(0:L_paths-1) / decay);
    powers = powers / sum(powers);
    for l = 1:L_paths
        if l == 1
            tau = 0;
            h   = sqrt(powers(l));
        else
            tau = randi([1, tau_max]);
            h   = (randn + 1j*randn)/sqrt(2) * sqrt(powers(l));
        end
        buf(tau+1 : tau+N) = buf(tau+1 : tau+N) + h * s;
    end
    y = buf(1:N);
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
