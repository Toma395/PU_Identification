% experiments/compare_methods.m — 比較実験B: 手法 × SNR 識別精度比較
% 3手法（尖度のみ / プリアンブルのみ / 結合）を評価
% SNR: -20〜+20dB 1dB刻み, 試行数500/クラス

base_dir = fileparts(which('compare_methods'));
addpath(fullfile(base_dir, '..', 'signals'));
addpath(fullfile(base_dir, '..', 'features'));
addpath(fullfile(base_dir, '..', 'classifier'));
addpath(fullfile(base_dir, '..'));   % apply_paper_style

%% パラメータ
snr_list  = -20:1:20;    % 41点、遷移域を細かくカバー
N_trials  = 500;         % クラスあたり試行数
N_snr     = length(snr_list);
N_methods = 3;

[~, ep_ref] = gen_ETC('seed', 1);
uw  = ep_ref.uw;
sps = ep_ref.sps;

%% 特徴量収集
scores      = zeros(N_trials*2, N_methods, N_snr);
labels_true = zeros(N_trials*2, N_snr);

fprintf('特徴量計算中 (%d SNR × %d試行 × 2クラス) ...\n', N_snr, N_trials);
for si = 1:N_snr
    snr = snr_list(si);
    for ti = 1:N_trials
        seed_etc = ti * 37  + si * 10000;         % ETC 信号シード
        seed_uav = ti * 41  + si * 10000 + 9999983; % UAV 信号シード（独立）

        [etc, ~] = gen_ETC('seed', seed_etc);
        etc_n = awgn_add(etc, snr, true);
        k = calc_kurtosis(etc_n);
        p = calc_preamble(etc_n, uw, sps);
        scores(ti, :, si)        = make_scores(k, p);
        labels_true(ti, si)      = 1;

        [uav, ~] = gen_UAV('seed', seed_uav);
        uav_n = awgn_add(uav, snr, false);
        k = calc_kurtosis(uav_n);
        p = calc_preamble(uav_n, uw, sps);
        scores(N_trials+ti, :, si)   = make_scores(k, p);
        labels_true(N_trials+ti, si) = 0;
    end
    fprintf('  SNR=%+4ddB 完了\n', snr);
end

%% 精度・F1 集計
thresh_fixed = [1.0, 0.5, 0.5];
acc_mat = zeros(N_snr, N_methods);
f1_mat  = zeros(N_snr, N_methods);
for si = 1:N_snr
    y_true = labels_true(:, si);
    for mi = 1:N_methods
        y_pred = scores(:, mi, si) >= thresh_fixed(mi);
        [acc_mat(si,mi), f1_mat(si,mi)] = calc_acc_f1(y_true, y_pred);
    end
end

%% 結果表示
method_names = {'Kurtosis','Preamble','Combined'};
fprintf('\n=== 識別精度 (Accuracy) ===\n');
fprintf('%-8s | %-10s %-10s %-10s\n', 'SNR[dB]', method_names{:});
fprintf('%s\n', repmat('-',1,44));
for si = 1:N_snr
    fprintf('%-8d | %-10.4f %-10.4f %-10.4f\n', snr_list(si), acc_mat(si,:));
end

fprintf('\n=== 90%%精度を維持できる最低SNR ===\n');
for mi = 1:N_methods
    idx = find(acc_mat(:,mi) >= 0.90, 1, 'first');
    if isempty(idx)
        fprintf('  %-10s: 全SNRで90%%未達\n', method_names{mi});
    else
        fprintf('  %-10s: SNR >= %+d dB で 90%%維持\n', method_names{mi}, snr_list(idx));
    end
end

%% AUC集計
fprintf('\n=== AUC ===\n');
fprintf('%-8s | %-10s %-10s %-10s\n', 'SNR[dB]', method_names{:});
fprintf('%s\n', repmat('-',1,44));
for si = 1:N_snr
    y_true = labels_true(:, si);
    aucs = zeros(1, N_methods);
    for mi = 1:N_methods
        [fpr, tpr] = compute_roc(y_true, scores(:,mi,si));
        aucs(mi) = trapz(fpr, tpr);
    end
    fprintf('%-8d | %-10.4f %-10.4f %-10.4f\n', snr_list(si), aucs);
end

%% プロット① ROC曲線（SNR = 0dB, -5dB）
roc_snr  = [0, -5];
colors   = lines(N_methods);
lstyles  = {'-','--','-.'};

fig1 = figure('Visible','off','Color','w','Position',[0 0 1000 480]);
for pi = 1:2
    si = find(snr_list == roc_snr(pi), 1);
    subplot(1,2,pi);
    hold on; grid on; box on;
    for mi = 1:N_methods
        y_true = labels_true(:, si);
        [fpr, tpr] = compute_roc(y_true, scores(:,mi,si));
        auc = trapz(fpr, tpr);
        plot(fpr, tpr, lstyles{mi}, 'Color',colors(mi,:), 'LineWidth',2.0, ...
            'DisplayName', sprintf('%s (AUC=%.3f)', method_names{mi}, auc));
    end
    plot([0 1],[0 1],'k:','LineWidth',1.0,'HandleVisibility','off');
    xlabel('FPR','FontSize',14); ylabel('TPR','FontSize',14);
    title(sprintf('ROC @ SNR=%+ddB', roc_snr(pi)),'FontSize',15);
    legend('Location','southeast','FontSize',11);
    xlim([0 1]); ylim([0 1]);
end
apply_paper_style(fig1);
saveas(fig1, fullfile(base_dir, 'roc_curves.png'));
fprintf('\nROC曲線保存: experiments/roc_curves.png\n');

%% プロット② 精度 vs SNR
fig2 = figure('Visible','off','Color','w','Position',[0 0 900 500]);
hold on; grid on; box on;
for mi = 1:N_methods
    plot(snr_list, acc_mat(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',5,'DisplayName',method_names{mi});
end
yline(0.9,'k--','90%','LabelHorizontalAlignment','right','FontSize',12,'HandleVisibility','off');
xlabel('SNR [dB]','FontSize',14);
ylabel('Accuracy','FontSize',14);
title(sprintf('識別精度 vs SNR  (AWGN, N_{trials}=%d)', N_trials),'FontSize',16);
legend('Location','southeast','FontSize',12);
ylim([0.4 1.05]); xlim([snr_list(1)-1, snr_list(end)+1]);
apply_paper_style(fig2);
saveas(fig2, fullfile(base_dir, 'accuracy_vs_snr.png'));
fprintf('精度グラフ保存: experiments/accuracy_vs_snr.png\n');

disp('compare_methods: OK');


%% ローカル関数
function sc = make_scores(kurt, p_peak)
    k_norm = max(0, min(1, (3 - kurt) / 2));
    sc = [3 - kurt, p_peak, 0.5*k_norm + 0.5*p_peak];
end

function y = awgn_add(x, snr_db, is_real)
    pwr  = mean(abs(x).^2);
    npwr = pwr / 10^(snr_db/10);
    if is_real
        y = x + sqrt(npwr) * randn(size(x));
    else
        y = x + sqrt(npwr/2) * (randn(size(x)) + 1j*randn(size(x)));
    end
end

function [acc, f1] = calc_acc_f1(y_true, y_pred)
    tp  = sum(y_true==1 & y_pred==1);
    tn  = sum(y_true==0 & y_pred==0);
    fp  = sum(y_true==0 & y_pred==1);
    fn  = sum(y_true==1 & y_pred==0);
    acc = (tp+tn) / (tp+tn+fp+fn);
    pre = tp / max(tp+fp, 1);
    rec = tp / max(tp+fn, 1);
    f1  = 2*pre*rec / max(pre+rec, eps);
end

function [fpr_out, tpr_out] = compute_roc(y_true, scores)
    thresholds = [max(scores)+1e-9; sort(unique(scores),'descend'); min(scores)-1e-9];
    N  = length(thresholds);
    Np = sum(y_true==1);
    Nn = sum(y_true==0);
    fpr_out = zeros(N,1);
    tpr_out = zeros(N,1);
    for i = 1:N
        yp = scores >= thresholds(i);
        tp = sum(y_true==1 & yp==1);
        fp = sum(y_true==0 & yp==1);
        tpr_out(i) = tp / max(Np,1);
        fpr_out(i) = fp / max(Nn,1);
    end
end
