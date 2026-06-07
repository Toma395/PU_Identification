% experiments/awgn_auc_vs_snr.m
% AWGN SNR sweep AUC分析（compare_methods.m から独立）
%
% compare_methods.m と同一シード・同一パラメータで AWGN sweep スコアを再計算し、
% perfcurve で AUC を算出。閾値非依存の識別能力を評価する。
%
% 出力: experiments/awgn_auc_vs_snr.png
%   左パネル: AUC vs SNR (3手法, AUC=0.9破線)
%   中パネル: ROC @ SNR=0dB
%   右パネル: ROC @ SNR=-10dB

base_dir = fileparts(which('awgn_auc_vs_snr'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..'));   % apply_paper_style

%% パラメータ（compare_methods.m と同一）
snr_list = -20:1:20;
N_snr    = length(snr_list);
N_trials = 500;
mnames   = {'Kurtosis','Preamble','Combined'};

[~, ep_ref] = gen_ETC('seed', 1);
uw  = ep_ref.uw;
sps = ep_ref.sps;

%% スコア計算（compare_methods.m と同一ループ）
labels_true = [ones(N_trials,1); zeros(N_trials,1)];
scores      = zeros(N_trials*2, 3, N_snr);
auc_mat     = zeros(N_snr, 3);

roc_snrs   = [0, -10];      % ROC描画対象 SNR
scores_roc = cell(1, length(roc_snrs));

fprintf('AWGN スコア計算中 (%d SNR点 × %d試行/クラス) ...\n', N_snr, N_trials);
for si = 1:N_snr
    snr = snr_list(si);
    for ti = 1:N_trials
        seed_etc = ti*37  + si*10000;
        seed_uav = ti*41  + si*10000 + 9999983;

        [etc, ~] = gen_ETC('seed', seed_etc);
        etc_n = awgn_add(etc, snr, true);
        k = calc_kurtosis(etc_n);
        p = calc_preamble(etc_n, uw, sps);
        scores(ti,:,si) = make_scores(k, p);

        [uav, ~] = gen_UAV('seed', seed_uav);
        uav_n = awgn_add(uav, snr, false);
        k = calc_kurtosis(uav_n);
        p = calc_preamble(uav_n, uw, sps);
        scores(N_trials+ti,:,si) = make_scores(k, p);
    end

    for mi = 1:3
        [~,~,~, auc_mat(si,mi)] = perfcurve(labels_true, scores(:,mi,si), 1);
    end

    ric = find(roc_snrs == snr, 1);
    if ~isempty(ric)
        scores_roc{ric} = scores(:,:,si);
    end
    fprintf('  SNR=%+4ddB 完了\n', snr);
end

%% サマリ
fprintf('\n=== AWGN AUC vs SNR（perfcurve） ===\n');
fprintf('%-8s | %-12s %-12s %-12s\n', 'SNR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_snr
    fprintf('%-8d | %-12.3f %-12.3f %-12.3f\n', snr_list(si), auc_mat(si,:));
end

fprintf('\nCombined AUC < Preamble AUC の条件:\n');
found = false;
for si = 1:N_snr
    diff_auc = auc_mat(si,2) - auc_mat(si,3);
    if diff_auc > 0.005
        fprintf('  SNR=%+4ddB: Preamble=%.3f > Combined=%.3f  (差=%.3f)\n', ...
            snr_list(si), auc_mat(si,2), auc_mat(si,3), diff_auc);
        found = true;
    end
end
if ~found
    fprintf('  （該当条件なし — Combined AUC が Preamble を下回る帯なし）\n');
end

%% プロット
lstyles = {'-','--','-.'};
colors  = lines(3);

fig = figure('Visible','off','Color','w','Position',[0 0 1400 450]);

% ─── 左: AUC vs SNR ───────────────────────────────────────────
subplot(1,3,1);
hold on; grid on; box on;
for mi = 1:3
    plot(snr_list, auc_mat(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',4,'DisplayName',mnames{mi});
end
yline(0.9,'k--','AUC=0.9','LabelHorizontalAlignment','right','FontSize',11,...
    'HandleVisibility','off');
xlabel('SNR [dB]  (AWGN)','FontSize',13);
ylabel('AUC','FontSize',13);
title('AUC vs SNR（閾値非依存）','FontSize',14);
legend('Location','southeast','FontSize',11);
xlim([snr_list(1)-1  snr_list(end)+1]);
ylim([0.4 1.05]);

% ─── 中・右: ROC @ SNR=0dB / -10dB ──────────────────────────
roc_titles = {sprintf('ROC @ SNR=%+ddB', roc_snrs(1)), ...
              sprintf('ROC @ SNR=%+ddB', roc_snrs(2))};
for ri = 1:2
    subplot(1,3,ri+1);
    hold on; grid on; box on;
    if ~isempty(scores_roc{ri})
        for mi = 1:3
            [fpr_r, tpr_r, ~, auc_r] = perfcurve(labels_true, scores_roc{ri}(:,mi), 1);
            plot(fpr_r, tpr_r, lstyles{mi}, 'Color',colors(mi,:), 'LineWidth',2.0, ...
                'DisplayName', sprintf('%s (AUC=%.3f)', mnames{mi}, auc_r));
        end
    end
    plot([0 1],[0 1],'k:','LineWidth',1.0,'HandleVisibility','off');
    xlabel('FPR','FontSize',13);
    ylabel('TPR','FontSize',13);
    title(roc_titles{ri},'FontSize',13);
    legend('Location','southeast','FontSize',10);
    xlim([0 1]); ylim([0 1]);
end

apply_paper_style(fig);
sgtitle(fig, 'AWGN sweep AUC分析: Preamble vs Combined（閾値非依存）', ...
    'FontSize',12,'Color',[0.15 0.15 0.15]);
saveas(fig, fullfile(base_dir,'awgn_auc_vs_snr.png'));
fprintf('\n保存: experiments/awgn_auc_vs_snr.png\n');
disp('awgn_auc_vs_snr: OK');


%% ─── ローカル関数（compare_methods.m と同一実装） ──────────────
function sc = make_scores(kurt, p_peak)
    k_norm = max(0, min(1, (3-kurt)/2));
    sc = [3-kurt,  p_peak,  0.5*k_norm + 0.5*p_peak];
end

function y = awgn_add(x, snr_db, is_real)
    pwr  = mean(abs(x).^2);
    npwr = pwr / 10^(snr_db/10);
    if is_real
        y = x + sqrt(npwr)*randn(size(x));
    else
        y = x + sqrt(npwr/2)*(randn(size(x)) + 1j*randn(size(x)));
    end
end
