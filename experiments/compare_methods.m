% experiments/compare_methods.m — 比較実験B: 手法 × SNR 識別精度比較
% 3手法（尖度のみ / プリアンブルのみ / 結合）を評価
% 評価指標: 識別精度, F1スコア, ROC曲線（AUC）

base_dir = fileparts(which('compare_methods'));
addpath(fullfile(base_dir, '..', 'signals'));
addpath(fullfile(base_dir, '..', 'features'));
addpath(fullfile(base_dir, '..', 'classifier'));

%% パラメータ
snr_list  = [-10, -5, 0, 5, 10, 15, 20];
N_trials  = 100;    % クラスあたり試行数（ETC×100 + UAV×100 = 200試行/SNR）
N_snr     = length(snr_list);
N_methods = 3;

% テンプレートETC生成（UW・sps取得用）
[~, ep_ref] = gen_ETC('seed', 1);
uw  = ep_ref.uw;
sps = ep_ref.sps;

%% 特徴量収集
% scores: (N_trials*2, N_methods, N_snr)  高いほどETC
scores      = zeros(N_trials*2, N_methods, N_snr);
labels_true = zeros(N_trials*2, N_snr);  % 1=ETC, 0=UAV

fprintf('特徴量計算中 (%d SNR × %d試行 × 2クラス) ...\n', N_snr, N_trials);
for si = 1:N_snr
    snr = snr_list(si);
    for ti = 1:N_trials
        seed = ti * 37;

        % ETC試行 (true label = 1)
        [etc, ~] = gen_ETC('seed', seed);
        etc_n = awgn_add(etc, snr, true);
        k = calc_kurtosis(etc_n);
        p = calc_preamble(etc_n, uw, sps);
        scores(ti, :, si)        = make_scores(k, p);
        labels_true(ti, si)      = 1;

        % UAV試行 (true label = 0)
        [uav, ~] = gen_UAV('seed', seed);
        uav_n = awgn_add(uav, snr, false);
        k = calc_kurtosis(uav_n);
        p = calc_preamble(uav_n, uw, sps);
        scores(N_trials+ti, :, si)   = make_scores(k, p);
        labels_true(N_trials+ti, si) = 0;
    end
    fprintf('  SNR=%4ddB 完了\n', snr);
end

%% 精度・F1 集計
thresh_fixed = [1.0, 0.5, 0.5];  % kurtosis, preamble, combined
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
method_names = {'Kurtosis only', 'Preamble only', 'Combined    '};

fprintf('\n=== 識別精度 (Accuracy) ===\n');
fprintf('%-8s | %-14s %-14s %-14s\n', 'SNR[dB]', method_names{:});
fprintf('%s\n', repmat('-',1,56));
for si = 1:N_snr
    fprintf('%-8d | %-14.4f %-14.4f %-14.4f\n', snr_list(si), acc_mat(si,:));
end

fprintf('\n=== F1スコア ===\n');
fprintf('%-8s | %-14s %-14s %-14s\n', 'SNR[dB]', method_names{:});
fprintf('%s\n', repmat('-',1,56));
for si = 1:N_snr
    fprintf('%-8d | %-14.4f %-14.4f %-14.4f\n', snr_list(si), f1_mat(si,:));
end

% 90%精度を超える最小SNR
fprintf('\n=== 識別精度90%%を超える最小SNR ===\n');
for mi = 1:N_methods
    idx = find(acc_mat(:,mi) >= 0.90, 1, 'first');
    if isempty(idx)
        fprintf('  %-16s: 90%%未達\n', method_names{mi});
    else
        fprintf('  %-16s: %4d dB\n', method_names{mi}, snr_list(idx));
    end
end

%% AUC集計
fprintf('\n=== AUC (ROC面積) ===\n');
fprintf('%-8s | %-14s %-14s %-14s\n', 'SNR[dB]', method_names{:});
fprintf('%s\n', repmat('-',1,56));
for si = 1:N_snr
    y_true = labels_true(:, si);
    aucs = zeros(1, N_methods);
    for mi = 1:N_methods
        [fpr, tpr] = compute_roc(y_true, scores(:,mi,si));
        aucs(mi) = trapz(fpr, tpr);
    end
    fprintf('%-8d | %-14.4f %-14.4f %-14.4f\n', snr_list(si), aucs);
end

%% ROC曲線プロット（SNR = 0dB, -5dB）
roc_snr   = [0, -5];
roc_idx   = arrayfun(@(s) find(snr_list==s,1), roc_snr);
colors    = lines(N_methods);
lstyles   = {'-', '--', '-.'};

fig = figure('Visible', 'off');
for pi = 1:length(roc_idx)
    si     = roc_idx(pi);
    y_true = labels_true(:, si);
    subplot(1, length(roc_idx), pi);
    hold on; grid on; box on;
    for mi = 1:N_methods
        [fpr, tpr] = compute_roc(y_true, scores(:,mi,si));
        auc = trapz(fpr, tpr);
        plot(fpr, tpr, lstyles{mi}, 'Color', colors(mi,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s (AUC=%.3f)', strtrim(method_names{mi}), auc));
    end
    plot([0 1],[0 1],'k:','HandleVisibility','off');
    xlabel('FPR'); ylabel('TPR');
    title(sprintf('ROC @ SNR=%ddB', snr_list(si)));
    legend('Location','southeast','FontSize',7);
    xlim([0 1]); ylim([0 1]);
end
saveas(fig, fullfile(base_dir, 'roc_curves.png'));
fprintf('\nROC曲線保存: experiments/roc_curves.png\n');

%% 精度 vs SNR プロット
fig2 = figure('Visible', 'off');
hold on; grid on; box on;
for mi = 1:N_methods
    plot(snr_list, acc_mat(:,mi), lstyles{mi}, 'Color', colors(mi,:), ...
        'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5, ...
        'DisplayName', strtrim(method_names{mi}));
end
yline(0.9, 'k--', '90%', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
xlabel('SNR [dB]'); ylabel('Accuracy'); title('識別精度 vs SNR');
legend('Location','southeast'); ylim([0 1]);
saveas(fig2, fullfile(base_dir, 'accuracy_vs_snr.png'));
fprintf('精度グラフ保存: experiments/accuracy_vs_snr.png\n');

disp('Phase 3 compare_methods: OK');


%% ローカル関数
function sc = make_scores(kurt, p_peak)
    % 各スコア: 高いほどETC
    k_norm = max(0, min(1, (3 - kurt) / 2));
    sc = [3 - kurt, ...                       % kurtosis
          p_peak,   ...                       % preamble
          0.5*k_norm + 0.5*p_peak];           % combined
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
    % 閾値を高→低で走査 → FPR: 0→1, TPR: 0→1
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
