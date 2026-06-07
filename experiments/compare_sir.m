% experiments/compare_sir.m
% 比較実験: 同一チャネル混信環境でのETC検出性能
%
% 問い1: SIR（ETC電力/UAV干渉電力）を振ってETC検出精度を評価
% 問い2: 干渉UAV数を増やしてETC検出精度の変化を評価
%
% 改訂: SIR -20〜+20dB 1dB刻み(41点), 試行数500, n_uav 1〜15

base_dir = fileparts(which('compare_sir'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..'));          % apply_paper_style

%% パラメータ
sir_list   = -20:1:20;  % 昇順 41点
N_sir      = length(sir_list);
N_trials   = 500;    % 試行数/クラス
snr_awgn   = 20;     % dB (AWGN固定)
thresh     = [1.0, 0.5, 0.5];

n_uav_list = 1:15;   % 問い2: 干渉UAV数
N_nuav     = length(n_uav_list);
sir_q2_db  = 0;      % 問い2固定SIR (0dB/UAV)

mnames = {'Kurtosis','Preamble','Combined'};

[~, ep] = gen_ETC('seed', 1);
uw  = ep.uw;
sps = ep.sps;
RS_P = 11; RS_Q = 5;  % 20MHz → 44MHz

%% 問い1: SIR vs ETC検出精度
fprintf('=== 問い1: SIR vs ETC検出精度 (%d条件 × %d試行/クラス) ===\n', N_sir, N_trials);

scores_q1 = zeros(N_trials*2, 3, N_sir);
labels_q1 = [ones(N_trials,1); zeros(N_trials,1)];
acc_q1    = zeros(N_sir, 3);
f1_q1     = zeros(N_sir, 3);
auc_q1    = zeros(N_sir, 3);

for si = 1:N_sir
    sir_db  = sir_list(si);
    sir_lin = 10^(sir_db/10);

    for ti = 1:N_trials
        seed = ti*37 + si*1000;

        [s_etc, ~] = gen_ETC('seed', seed);
        [s_uav, ~] = gen_UAV('seed', seed + 500013);
        y_pos = mix_etc_uav(s_etc, s_uav, sir_lin, snr_awgn, RS_P, RS_Q);
        k = calc_kurtosis(y_pos);
        p = calc_preamble(y_pos, uw, sps);
        scores_q1(ti,:,si) = make_scores(k,p);

        [s_uav2, ~] = gen_UAV('seed', seed + 1000013);
        y_neg = uav_only(s_uav2, length(s_etc), snr_awgn, RS_P, RS_Q);
        k = calc_kurtosis(y_neg);
        p = calc_preamble(y_neg, uw, sps);
        scores_q1(N_trials+ti,:,si) = make_scores(k,p);
    end

    for mi = 1:3
        [acc_q1(si,mi), f1_q1(si,mi)] = calc_acc_f1(labels_q1, scores_q1(:,mi,si) >= thresh(mi));
    end
    fprintf('  SIR=%+4ddB  完了\n', sir_db);
end

% Q1 AUC（scores_q1 に全スコア保存済み）
for si = 1:N_sir
    for mi = 1:3
        [fpr_a, tpr_a] = roc_curve(labels_q1, scores_q1(:,mi,si));
        auc_q1(si,mi)  = trapz(fpr_a, tpr_a);
    end
end

%% 問い2: 干渉UAV数 vs 検出精度
fprintf('\n=== 問い2: 干渉UAV数 vs 検出精度 (SIR=%.0fdB/UAV) ===\n', sir_q2_db);

acc_q2  = zeros(N_nuav, 3);
auc_q2  = zeros(N_nuav, 3);
sir_eff = sir_q2_db - 10*log10(n_uav_list);   % 実効SIR [dB]
kurt_q2_etc   = zeros(N_nuav, 1);   % 尖度診断用（ETC正例）
kurt_q2_uav   = zeros(N_nuav, 1);   % 尖度診断用（UAV負例）
roc_nuav_list = [8, 15];            % ROC曲線描画用n_uav
scores_roc_q2 = cell(1, length(roc_nuav_list));
labels_roc_q2 = cell(1, length(roc_nuav_list));

for ni = 1:N_nuav
    n_uav   = n_uav_list(ni);
    sir_lin = 10^(sir_q2_db/10);

    scores_q2    = zeros(N_trials*2, 3);
    labels_q2    = [ones(N_trials,1); zeros(N_trials,1)];
    kurt_buf_etc = zeros(N_trials, 1);
    kurt_buf_uav = zeros(N_trials, 1);

    for ti = 1:N_trials
        seed = ti*41 + ni*1000;

        [s_etc, ~] = gen_ETC('seed', seed);
        s_etc_n = s_etc / sqrt(mean(s_etc.^2));
        L = length(s_etc_n);   % ETC全長（11968サンプル）

        % 正例: ETC + n_uav×UAV（ETC全長Lを覆う）
        % ※修正前はUAVが先頭7040サンプルのみ覆い末尾4928が純ETCになっていた
        y_pos = sqrt(sir_lin) * s_etc_n;
        for ku = 1:n_uav
            s_u_long = gen_uav_long(L, seed + ku*300 + 2000000, RS_P, RS_Q);
            s_u_n    = s_u_long / sqrt(mean(abs(s_u_long).^2));
            y_pos    = y_pos + s_u_n;
        end
        y_pos = add_awgn(y_pos, snr_awgn);
        k = calc_kurtosis(y_pos);
        kurt_buf_etc(ti) = k;
        p = calc_preamble(y_pos, uw, sps);
        scores_q2(ti,:) = make_scores(k,p);

        % 負例: n_uav×UAVのみ（全長L、ゼロ埋めなし）
        % ※修正前はzeros(L,1)を基底にしていたため末尾4928がゼロ埋めになっていた
        y_neg = zeros(L, 1);
        for ku = 1:n_uav
            s_u_long = gen_uav_long(L, seed + ku*300 + 5000000, RS_P, RS_Q);
            s_u_n    = s_u_long / sqrt(mean(abs(s_u_long).^2));
            y_neg    = y_neg + s_u_n;
        end
        y_neg = add_awgn(y_neg, snr_awgn);
        k = calc_kurtosis(y_neg);
        kurt_buf_uav(ti) = k;
        p = calc_preamble(y_neg, uw, sps);
        scores_q2(N_trials+ti,:) = make_scores(k,p);
    end

    for mi = 1:3
        [acc_q2(ni,mi), ~] = calc_acc_f1(labels_q2, scores_q2(:,mi) >= thresh(mi));
    end
    kurt_q2_etc(ni) = mean(kurt_buf_etc);
    kurt_q2_uav(ni) = mean(kurt_buf_uav);

    % AUC計算（閾値非依存）
    for mi = 1:3
        [fpr_a, tpr_a] = roc_curve(labels_q2, scores_q2(:,mi));
        auc_q2(ni,mi)  = trapz(fpr_a, tpr_a);
    end

    % ROC描画用スコア保存
    ric = find(roc_nuav_list == n_uav, 1);
    if ~isempty(ric)
        scores_roc_q2{ric} = scores_q2;
        labels_roc_q2{ric} = labels_q2;
    end

    fprintf('  n_UAV=%2d  実効SIR=%5.1fdB  完了\n', n_uav, sir_eff(ni));
end

%% サマリ
fprintf('\n=== 問い1: SIR vs ETC検出精度（Accuracy）===\n');
fprintf('%-8s | %-12s %-12s %-12s\n', 'SIR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_sir
    fprintf('%-8d | %-12.3f %-12.3f %-12.3f\n', sir_list(si), acc_q1(si,:));
end

fprintf('\n=== 問い1: F1スコア ===\n');
fprintf('%-8s | %-12s %-12s %-12s\n', 'SIR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_sir
    fprintf('%-8d | %-12.3f %-12.3f %-12.3f\n', sir_list(si), f1_q1(si,:));
end

fprintf('\n=== 問い2: 干渉UAV数 vs 検出精度 ===\n');
fprintf('%-6s %-10s | %-12s %-12s %-12s\n', 'n_UAV', '実効SIR', mnames{:});
fprintf('%s\n', repmat('-',1,56));
for ni = 1:N_nuav
    fprintf('%-6d %-10.1f | %-12.3f %-12.3f %-12.3f\n', ...
        n_uav_list(ni), sir_eff(ni), acc_q2(ni,:));
end

fprintf('\n=== 問い2: 尖度診断（修正後・UAV干渉がETC全長を覆う）===\n');
fprintf('%-6s %-10s | %-12s %-12s  %-6s\n', 'n_UAV', '実効SIR', 'Kurt(ETC)', 'Kurt(UAV)', '差');
fprintf('%s\n', repmat('-',1,54));
for ni = 1:N_nuav
    fprintf('%-6d %-10.1f | %-12.3f %-12.3f  %+.3f\n', ...
        n_uav_list(ni), sir_eff(ni), kurt_q2_etc(ni), kurt_q2_uav(ni), ...
        kurt_q2_etc(ni) - kurt_q2_uav(ni));
end

fprintf('\n=== 問い1: SIR vs AUC（閾値非依存） ===\n');
fprintf('%-8s | %-12s %-12s %-12s\n', 'SIR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_sir
    fprintf('%-8d | %-12.3f %-12.3f %-12.3f\n', sir_list(si), auc_q1(si,:));
end

fprintf('\n=== 問い2: 干渉UAV数 vs AUC（閾値非依存・Preamble/Combined等価性検証） ===\n');
fprintf('%-6s %-10s | %-12s %-12s %-12s\n', 'n_UAV', '実効SIR', mnames{:});
fprintf('%s\n', repmat('-',1,56));
for ni = 1:N_nuav
    fprintf('%-6d %-10.1f | %-12.3f %-12.3f %-12.3f\n', ...
        n_uav_list(ni), sir_eff(ni), auc_q2(ni,:));
end
fprintf('\n理論予測: k_norm≈0(高n_uav)のとき combined=0.5×preamble → AUC(combined)=AUC(preamble)\n');
fprintf('         combined Accuracy=1.0 は閾値効果であり特徴融合ゲインではない\n\n');

fprintf('\n=== 90%%検出維持の下限SIR（問い1） ===\n');
for mi = 1:3
    idx = find(acc_q1(:,mi) >= 0.90, 1, 'first');  % 昇順なので first = 最低SIR
    if isempty(idx)
        fprintf('  %-10s: 全条件で90%%未達\n', mnames{mi});
    else
        fprintf('  %-10s: SIR >= %+ddB で90%%維持\n', mnames{mi}, sir_list(idx));
    end
end

%% プロット
lstyles = {'-','--','-.'};
colors  = lines(3);

%% ① SIR vs 精度・F1（2段）
fig1 = figure('Visible','off','Color','w','Position',[0 0 900 700]);

subplot(2,1,1);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list, acc_q1(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',5,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
yline(0.5,'k:','HandleVisibility','off');
xlabel('SIR [dB]  (ETC電力 / UAV干渉電力)','FontSize',14);
ylabel('Accuracy','FontSize',14);
title('問い1: SIR vs ETC検出精度  (AWGN SNR=20dB固定)','FontSize',16);
legend('Location','southeast','FontSize',12); ylim([0.4 1.05]);
xlim([sir_list(1)-1, sir_list(end)+1]);

subplot(2,1,2);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list, f1_q1(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','s','MarkerSize',5,'DisplayName',mnames{mi});
end
yline(0.9,'k--','F1=0.9','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
xlabel('SIR [dB]','FontSize',14);
ylabel('F1 Score','FontSize',14);
title('SIR vs F1スコア','FontSize',16);
legend('Location','southeast','FontSize',12); ylim([0.4 1.05]);
xlim([sir_list(1)-1, sir_list(end)+1]);

apply_paper_style(fig1);
saveas(fig1, fullfile(base_dir,'sir_accuracy.png'));
fprintf('\n保存: experiments/sir_accuracy.png\n');

%% ② ROC曲線（問い1: SIR=0dBと-5dB）
roc_sirs = [0, -5];
fig2 = figure('Visible','off','Color','w','Position',[0 0 900 420]);
for pi = 1:2
    si = find(sir_list == roc_sirs(pi), 1);
    subplot(1,2,pi);
    hold on; grid on; box on;
    for mi = 1:3
        [fpr,tpr] = roc_curve(labels_q1, scores_q1(:,mi,si));
        auc = trapz(fpr,tpr);
        plot(fpr, tpr, lstyles{mi}, 'Color',colors(mi,:), 'LineWidth',2.0, ...
            'DisplayName', sprintf('%s (AUC=%.3f)', mnames{mi}, auc));
    end
    plot([0 1],[0 1],'k:','LineWidth',1.0,'HandleVisibility','off');
    xlabel('FPR','FontSize',14); ylabel('TPR','FontSize',14);
    title(sprintf('ROC @ SIR=%+ddB', roc_sirs(pi)),'FontSize',16);
    legend('Location','southeast','FontSize',11);
    xlim([0 1]); ylim([0 1]);
end
apply_paper_style(fig2);
saveas(fig2, fullfile(base_dir,'sir_roc.png'));
fprintf('保存: experiments/sir_roc.png\n');

%% ③ 干渉UAV数 vs 精度（問い2）
fig3 = figure('Visible','off','Color','w','Position',[0 0 900 620]);
main_pos = [0.12, 0.09, 0.78, 0.58];
ax3 = axes(fig3, 'Position', main_pos);
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
for mi = 1:3
    plot(ax3, n_uav_list, acc_q2(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',7,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
xlabel(ax3,'干渉UAV数','FontSize',14);
ylabel(ax3,'Accuracy','FontSize',14);
xlim(ax3, [0.5, 15.5]);
ax3.XTick = 1:15;
ylim(ax3, [0.4 1.05]);

% 上軸に実効SIR表示（2点おきに間引いて重なり防止）
ax3b = axes(fig3,'Position',main_pos, ...
    'XAxisLocation','top','YAxisLocation','right','Color','none');
tick_sel  = 1:2:15;   % 1,3,5,7,9,11,13,15
tick_uav  = n_uav_list(tick_sel);
ax3b.XTick      = tick_uav;
ax3b.XLim       = [0.5, 15.5];
ax3b.XTickLabel = arrayfun(@(x) sprintf('%.0fdB', sir_q2_db - 10*log10(x)), ...
    tick_uav, 'UniformOutput', false);
ax3b.FontSize   = 12;
xlabel(ax3b,'実効SIR [dB]','FontSize',13);
ax3b.YTick = []; ax3b.YColor = 'none';
axes(ax3);
legend(ax3,'Location','northeast','FontSize',12);
apply_paper_style(fig3);
ax3.Position  = main_pos;
ax3b.Position = main_pos;
sgtitle(fig3, sprintf('問い2: 干渉UAV数 vs ETC検出精度  (SIR=%.0fdB/UAV)', sir_q2_db), ...
    'FontSize', 14, 'Color', [0.15 0.15 0.15]);
saveas(fig3, fullfile(base_dir,'nuav_accuracy.png'));
fprintf('保存: experiments/nuav_accuracy.png\n');

%% ④ 問い2 AUC分析: Preamble/Combined等価性 + ROC（閾値非依存）
fig4 = figure('Visible','off','Color','w','Position',[0 0 1200 430]);

% 左: AUC vs n_uav
subplot(1,3,1);
hold on; grid on; box on;
for mi = 1:3
    plot(n_uav_list, auc_q2(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',6,'DisplayName',mnames{mi});
end
yline(0.9,'k--','AUC=0.9','LabelHorizontalAlignment','right','FontSize',11,'HandleVisibility','off');
xlabel('干渉UAV数','FontSize',13);
ylabel('AUC','FontSize',13);
title('AUC vs n_{UAV}（閾値非依存）','FontSize',14);
legend('Location','southwest','FontSize',11);
xlim([0.5 15.5]); ylim([0.4 1.05]);
xticks([1 5 8 10 15]);

% 中・右: ROC at n_uav = roc_nuav_list
roc_titles = {sprintf('ROC @ n_{UAV}=%d', roc_nuav_list(1)), ...
              sprintf('ROC @ n_{UAV}=%d（最悪条件）', roc_nuav_list(2))};
for ri = 1:length(roc_nuav_list)
    subplot(1,3,ri+1);
    hold on; grid on; box on;
    if ~isempty(scores_roc_q2{ri})
        for mi = 1:3
            [fpr_r, tpr_r] = roc_curve(labels_roc_q2{ri}, scores_roc_q2{ri}(:,mi));
            auc_r = trapz(fpr_r, tpr_r);
            plot(fpr_r, tpr_r, lstyles{mi}, 'Color',colors(mi,:), 'LineWidth',2.0, ...
                'DisplayName', sprintf('%s (AUC=%.3f)', mnames{mi}, auc_r));
        end
    end
    plot([0 1],[0 1],'k:','LineWidth',1.0,'HandleVisibility','off');
    xlabel('FPR','FontSize',13); ylabel('TPR','FontSize',13);
    title(roc_titles{ri},'FontSize',13);
    legend('Location','southeast','FontSize',10);
    xlim([0 1]); ylim([0 1]);
end

apply_paper_style(fig4);
sgtitle(fig4, '問い2 AUC分析: Preambleと結合の等価性検証（k_{norm}≈0の条件）', ...
    'FontSize',12,'Color',[0.15 0.15 0.15]);
saveas(fig4, fullfile(base_dir,'nuav_auc_roc.png'));
fprintf('保存: experiments/nuav_auc_roc.png\n');

disp('compare_sir: OK');


%% ローカル関数
function s_long = gen_uav_long(L_target, seed_base, rs_p, rs_q)
    % L_targetサンプルのUAV信号を生成（複数フレームを連結してETC全区間を覆う）
    s_long = zeros(L_target, 1) + 0j;
    offset = 0;
    frame_idx = 0;
    while offset < L_target
        [s_u, ~] = gen_UAV('seed', seed_base + frame_idx * 999983);
        s_u_rs   = resample(s_u, rs_p, rs_q);
        take     = min(length(s_u_rs), L_target - offset);
        s_long(offset+1 : offset+take) = s_u_rs(1:take);
        offset    = offset + take;
        frame_idx = frame_idx + 1;
    end
end

function y = mix_etc_uav(s_etc, s_uav, sir_lin, snr_db, rs_p, rs_q)
    s_uav_rs = resample(s_uav, rs_p, rs_q);
    s_etc_n  = s_etc / sqrt(mean(s_etc.^2));
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

function y = add_awgn(x, snr_db)
    p_sig   = mean(abs(x).^2);
    p_noise = p_sig / 10^(snr_db/10);
    if isreal(x)
        y = x + sqrt(p_noise)*randn(size(x));
    else
        y = x + sqrt(p_noise/2)*(randn(size(x)) + 1j*randn(size(x)));
    end
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

function [fpr_out, tpr_out] = roc_curve(y_true, scores)
    thr = [max(scores)+1e-9; sort(unique(scores),'descend'); min(scores)-1e-9];
    Np  = sum(y_true==1);  Nn = sum(y_true==0);
    fpr_out = zeros(length(thr),1);
    tpr_out = zeros(length(thr),1);
    for i = 1:length(thr)
        yp = scores >= thr(i);
        tpr_out(i) = sum(y_true==1 & yp==1)/max(Np,1);
        fpr_out(i) = sum(y_true==0 & yp==1)/max(Nn,1);
    end
end
