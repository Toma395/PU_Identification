% experiments/compare_sir.m
% 比較実験: 同一チャネル混信環境でのETC検出性能
%
% 問い1: SIR（ETC電力/UAV干渉電力）を振ってETC検出精度を評価
%   正例: ETC + UAV干渉（SIRスケール） → label=1 (ETC含む)
%   負例: UAV干渉のみ                   → label=0 (ETC含まない)
%
% 問い2: 干渉UAV数を増やしてETC検出精度の変化を評価
%   SIR=0dB/UAV固定 → 実効SIR = -10*log10(n_UAV) で劣化
%
% 信号混合の前提:
%   gen_ETC: fs=44MHz（ASK実数）  gen_UAV: fs=20MHz（OFDM複素）
%   → UAVを 20MHz→44MHz に resample(11,5) してから加算する
%   → 混合信号は44MHzの複素信号。calc_kurtosis は real() を使用、
%     calc_preamble も sps=44（44MHzベース）で動作するため整合している。
%
% 出力:
%   experiments/sir_accuracy.png  -- 問い1: SIR vs 精度・F1
%   experiments/sir_roc.png       -- 問い1: ROC曲線
%   experiments/nuav_accuracy.png -- 問い2: 干渉UAV数 vs 精度

base_dir = fileparts(which('compare_sir'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..'));          % apply_paper_style

%% ─────────────────────────────────────────────────────────────────
%% パラメータ
%% ─────────────────────────────────────────────────────────────────
sir_list   = [20, 10, 5, 0, -5, -10, -20];   % dB
N_sir      = length(sir_list);
N_trials   = 200;   % 試行数/クラス (正例×200 + 負例×200 = 400試行/条件)
snr_awgn   = 20;    % dB (AWGN固定: 干渉の純粋な効果を観察)
thresh     = [1.0, 0.5, 0.5];   % kurtosis / preamble / combined

n_uav_list = [1, 2, 3, 5, 10];  % 問い2: 干渉UAV数
N_nuav     = length(n_uav_list);
sir_q2_db  = 0;     % 問い2固定SIR (0dB/UAV)

mnames = {'Kurtosis','Preamble','Combined'};

% テンプレートETC (UW・sps取得)
[~, ep] = gen_ETC('seed', 1);
uw  = ep.uw;
sps = ep.sps;   % = round(44e6/1e6) = 44

% サンプルレート比: ETC=44MHz, UAV=20MHz → 44/20 = 11/5
RS_P = 11; RS_Q = 5;

%% ─────────────────────────────────────────────────────────────────
%% 問い1: SIR vs ETC検出精度
%% ─────────────────────────────────────────────────────────────────
fprintf('=== 問い1: SIR vs ETC検出精度 (%d条件 × %d試行/クラス) ===\n', N_sir, N_trials);

scores_q1 = zeros(N_trials*2, 3, N_sir);
labels_q1 = [ones(N_trials,1); zeros(N_trials,1)];
acc_q1    = zeros(N_sir, 3);
f1_q1     = zeros(N_sir, 3);

for si = 1:N_sir
    sir_db  = sir_list(si);
    sir_lin = 10^(sir_db/10);

    for ti = 1:N_trials
        seed = ti*37 + si*200;

        % ── 正例: ETC + UAV干渉 ──────────────────────────────────
        [s_etc, ~] = gen_ETC('seed', seed);
        [s_uav, ~] = gen_UAV('seed', seed + 500);
        y_pos = mix_etc_uav(s_etc, s_uav, sir_lin, snr_awgn, RS_P, RS_Q);
        k = calc_kurtosis(y_pos);
        p = calc_preamble(y_pos, uw, sps);
        scores_q1(ti,:,si) = make_scores(k,p);

        % ── 負例: UAV干渉のみ ────────────────────────────────────
        [s_uav2, ~] = gen_UAV('seed', seed + 1000);
        y_neg = uav_only(s_uav2, length(s_etc), snr_awgn, RS_P, RS_Q);
        k = calc_kurtosis(y_neg);
        p = calc_preamble(y_neg, uw, sps);
        scores_q1(N_trials+ti,:,si) = make_scores(k,p);
    end

    for mi = 1:3
        [acc_q1(si,mi), f1_q1(si,mi)] = calc_acc_f1(labels_q1, scores_q1(:,mi,si) >= thresh(mi));
    end
    fprintf('  SIR=%4ddB  完了\n', sir_db);
end

%% ─────────────────────────────────────────────────────────────────
%% 問い2: 干渉UAV数 vs 検出精度 (SIR=0dB/UAV)
%% ─────────────────────────────────────────────────────────────────
fprintf('\n=== 問い2: 干渉UAV数 vs 検出精度 (SIR=%.0fdB/UAV) ===\n', sir_q2_db);

acc_q2  = zeros(N_nuav, 3);
sir_eff = sir_q2_db - 10*log10(n_uav_list);   % 実効SIR [dB]

for ni = 1:N_nuav
    n_uav   = n_uav_list(ni);
    sir_lin = 10^(sir_q2_db/10);   % per-UAV SIR = 0dB → ETC電力=1, UAV各1

    scores_q2 = zeros(N_trials*2, 3);
    labels_q2 = [ones(N_trials,1); zeros(N_trials,1)];

    for ti = 1:N_trials
        seed = ti*41 + ni*300;

        % ── 正例: ETC + n_uav × UAV干渉 ─────────────────────────
        [s_etc, ~] = gen_ETC('seed', seed);
        s_etc_n = s_etc / sqrt(mean(s_etc.^2));
        L = length(s_etc_n);
        y_pos = sqrt(sir_lin) * s_etc_n;   % ETC成分（電力=sir_lin=1）

        for ku = 1:n_uav
            [s_u, ~] = gen_UAV('seed', seed + ku*300);
            s_u_rs = resample(s_u, RS_P, RS_Q);           % 20MHz → 44MHz
            s_u_n  = s_u_rs / sqrt(mean(abs(s_u_rs).^2)); % 単位電力
            Lu = min(length(s_u_n), L);
            y_pos(1:Lu) = y_pos(1:Lu) + s_u_n(1:Lu);     % 各UAV電力=1
        end
        y_pos = add_awgn(y_pos, snr_awgn);
        k = calc_kurtosis(y_pos);
        p = calc_preamble(y_pos, uw, sps);
        scores_q2(ti,:) = make_scores(k,p);

        % ── 負例: n_uav × UAV干渉のみ ───────────────────────────
        y_neg = zeros(L, 1);
        for ku = 1:n_uav
            [s_u, ~] = gen_UAV('seed', seed + ku*300 + 5000);
            s_u_rs = resample(s_u, RS_P, RS_Q);
            s_u_n  = s_u_rs / sqrt(mean(abs(s_u_rs).^2));
            Lu = min(length(s_u_n), L);
            y_neg(1:Lu) = y_neg(1:Lu) + s_u_n(1:Lu);
        end
        y_neg = add_awgn(y_neg, snr_awgn);
        k = calc_kurtosis(y_neg);
        p = calc_preamble(y_neg, uw, sps);
        scores_q2(N_trials+ti,:) = make_scores(k,p);
    end

    for mi = 1:3
        [acc_q2(ni,mi), ~] = calc_acc_f1(labels_q2, scores_q2(:,mi) >= thresh(mi));
    end
    fprintf('  n_UAV=%2d  実効SIR=%5.1fdB  完了\n', n_uav, sir_eff(ni));
end

%% ─────────────────────────────────────────────────────────────────
%% 統計サマリ
%% ─────────────────────────────────────────────────────────────────
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

fprintf('\n=== 90%%検出維持の限界SIR（問い1） ===\n');
for mi = 1:3
    idx = find(acc_q1(:,mi) >= 0.90, 1, 'last');
    if isempty(idx)
        fprintf('  %-10s: 全条件で90%%未達\n', mnames{mi});
    else
        fprintf('  %-10s: SIR = %+ddB まで\n', mnames{mi}, sir_list(idx));
    end
end

%% ─────────────────────────────────────────────────────────────────
%% プロット
%% ─────────────────────────────────────────────────────────────────
lstyles = {'-','--','-.'};
colors  = lines(3);

%% ① SIR vs 精度・F1（2段）
fig1 = figure('Visible','off','Color','w','Position',[0 0 900 700]);

subplot(2,1,1);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list, acc_q1(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',8,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','left','FontSize',12,'HandleVisibility','off');
yline(0.5,'k:','50%','LabelHorizontalAlignment','left','FontSize',12,'HandleVisibility','off');
xlabel('SIR [dB]  (ETC電力 / UAV干渉電力)','FontSize',14);
ylabel('Accuracy','FontSize',14);
title('問い1: SIR vs ETC検出精度  (AWGN SNR=20dB固定)','FontSize',16);
legend('Location','southeast','FontSize',12); ylim([0.4 1.05]);
xlim([sir_list(end)-2, sir_list(1)+2]);

subplot(2,1,2);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list, f1_q1(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','s','MarkerSize',8,'DisplayName',mnames{mi});
end
yline(0.9,'k--','F1=0.9','LabelHorizontalAlignment','left','FontSize',12,'HandleVisibility','off');
xlabel('SIR [dB]','FontSize',14);
ylabel('F1 Score','FontSize',14);
title('SIR vs F1スコア','FontSize',16);
legend('Location','southeast','FontSize',12); ylim([0.4 1.05]);
xlim([sir_list(end)-2, sir_list(1)+2]);

apply_paper_style(fig1);
saveas(fig1, fullfile(base_dir,'sir_accuracy.png'));
fprintf('\n保存: experiments/sir_accuracy.png\n');

%% ② ROC曲線（問い1: SIR=0dBと-5dBの2条件）
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
main_pos = [0.12, 0.09, 0.78, 0.58];   % 上マージン 33%確保
ax3 = axes(fig3, 'Position', main_pos);
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
for mi = 1:3
    plot(ax3, n_uav_list, acc_q2(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',8,'DisplayName',mnames{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
xlabel(ax3,'干渉UAV数','FontSize',14);
ylabel(ax3,'Accuracy','FontSize',14);
% title はsgtitleに移動（ax3 title と ax3b 上軸ラベルの重なり回避）
% 上軸に実効SIR表示
ax3b = axes(fig3,'Position',main_pos, ...
    'XAxisLocation','top','YAxisLocation','right','Color','none');
ax3b.XTick      = n_uav_list;
ax3b.XLim       = ax3.XLim;
ax3b.XTickLabel = arrayfun(@(x) sprintf('%.1fdB',x), sir_eff, 'UniformOutput',false);
ax3b.FontSize   = 13;
xlabel(ax3b,'実効SIR [dB]','FontSize',13);
ax3b.YTick = []; ax3b.YColor = 'none';
axes(ax3);
legend(ax3,'Location','northeast','FontSize',12);
apply_paper_style(fig3);
% apply_paper_style 後に位置を強制（FontSize変更による自動レイアウト崩れを防ぐ）
ax3.Position  = main_pos;
ax3b.Position = main_pos;
% sgtitle は figure 最上部（ax3b の上軸ラベルの更に上）に配置される
sgtitle(fig3, sprintf('問い2: 干渉UAV数 vs ETC検出精度  (SIR=%.0fdB/UAV)', sir_q2_db), ...
    'FontSize', 14, 'Color', [0.15 0.15 0.15]);
saveas(fig3, fullfile(base_dir,'nuav_accuracy.png'));
fprintf('保存: experiments/nuav_accuracy.png\n');

disp('compare_sir: OK');


%% ─────────────────────────────────────────────────────────────────
%% ローカル関数
%% ─────────────────────────────────────────────────────────────────

function y = mix_etc_uav(s_etc, s_uav, sir_lin, snr_db, rs_p, rs_q)
% UAVを rs_p:rs_q でリサンプルして同レート化し、SIR比で混合してAWGNを加える
    s_uav_rs = resample(s_uav, rs_p, rs_q);              % 20MHz → 44MHz
    s_etc_n  = s_etc / sqrt(mean(s_etc.^2));             % 単位電力（実数）
    s_uav_n  = s_uav_rs / sqrt(mean(abs(s_uav_rs).^2));  % 単位電力（複素）
    L = min(length(s_etc_n), length(s_uav_n));
    mixed = sqrt(sir_lin)*s_etc_n(1:L) + s_uav_n(1:L);  % ETC電力=sir_lin, UAV電力=1
    y = add_awgn(mixed, snr_db);
end

function y = uav_only(s_uav, L_ref, snr_db, rs_p, rs_q)
% UAV干渉のみ（負例）: リサンプル→単位電力→AWGNを加える
    s_uav_rs = resample(s_uav, rs_p, rs_q);
    s_uav_n  = s_uav_rs / sqrt(mean(abs(s_uav_rs).^2));
    L = min(length(s_uav_n), L_ref);
    y = add_awgn(s_uav_n(1:L), snr_db);
end

function y = add_awgn(x, snr_db)
% 信号全体電力基準でAWGNを加える（複素/実数両対応）
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
