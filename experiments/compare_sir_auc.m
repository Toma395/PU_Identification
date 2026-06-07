% experiments/compare_sir_auc.m
% Q1 AUC vs SIR 分析（compare_sir.m の accuracy 計算を変えずに独立実行）
%
% compare_sir.m と同一シード・同一パラメータで Q1 スコアを再計算し、
% perfcurve で AUC を算出。閾値に依存しない識別能力を評価する。
%
% 出力: experiments/sir_auc_roc.png
%   左パネル: AUC vs SIR (3手法)
%   中パネル: ROC @ SIR=-5dB
%   右パネル: ROC @ SIR=-20dB (最悪条件)

base_dir = fileparts(which('compare_sir_auc'));
addpath(fullfile(base_dir,'..','signals'));
addpath(fullfile(base_dir,'..','features'));
addpath(fullfile(base_dir,'..','classifier'));
addpath(fullfile(base_dir,'..'));   % apply_paper_style

%% パラメータ（compare_sir.m Q1 と同一）
sir_list = -20:1:20;
N_sir    = length(sir_list);
N_trials = 500;
snr_awgn = 20;
RS_P = 11; RS_Q = 5;
mnames   = {'Kurtosis','Preamble','Combined'};

[~, ep] = gen_ETC('seed', 1);
uw = ep.uw;  sps = ep.sps;

%% Q1 スコア計算（compare_sir.m Q1 と同一ループ）
labels_q1  = [ones(N_trials,1); zeros(N_trials,1)];
scores_q1  = zeros(N_trials*2, 3, N_sir);
auc_q1     = zeros(N_sir, 3);

roc_sirs   = [-5, -20];     % ROC 描画対象 SIR
scores_roc = cell(1, length(roc_sirs));   % 保存用

fprintf('Q1 スコア計算中 (%d SIR点 × %d試行/クラス) ...\n', N_sir, N_trials);
for si = 1:N_sir
    sir_db  = sir_list(si);
    sir_lin = 10^(sir_db/10);

    for ti = 1:N_trials
        seed = ti*37 + si*1000;   % compare_sir.m と同一

        [s_etc, ~] = gen_ETC('seed', seed);
        [s_uav,  ~] = gen_UAV('seed', seed + 500013);
        y_pos = mix_etc_uav(s_etc, s_uav, sir_lin, snr_awgn, RS_P, RS_Q);
        k = calc_kurtosis(y_pos);
        p = calc_preamble(y_pos, uw, sps);
        scores_q1(ti,:,si) = make_scores(k, p);

        [s_uav2, ~] = gen_UAV('seed', seed + 1000013);
        y_neg = uav_only(s_uav2, length(s_etc), snr_awgn, RS_P, RS_Q);
        k = calc_kurtosis(y_neg);
        p = calc_preamble(y_neg, uw, sps);
        scores_q1(N_trials+ti,:,si) = make_scores(k, p);
    end

    for mi = 1:3
        [~,~,~, auc_q1(si,mi)] = perfcurve(labels_q1, scores_q1(:,mi,si), 1);
    end

    ric = find(roc_sirs == sir_db, 1);
    if ~isempty(ric)
        scores_roc{ric} = scores_q1(:,:,si);
    end
    fprintf('  SIR=%+4ddB  完了\n', sir_db);
end

%% サマリ
fprintf('\n=== Q1 AUC vs SIR（perfcurve） ===\n');
fprintf('%-8s | %-12s %-12s %-12s\n', 'SIR[dB]', mnames{:});
fprintf('%s\n', repmat('-',1,48));
for si = 1:N_sir
    fprintf('%-8d | %-12.3f %-12.3f %-12.3f\n', sir_list(si), auc_q1(si,:));
end

fprintf('\nCombined AUC < Preamble AUC の条件:\n');
found = false;
for si = 1:N_sir
    diff_auc = auc_q1(si,2) - auc_q1(si,3);
    if diff_auc > 0.005
        fprintf('  SIR=%+4ddB: Preamble=%.3f > Combined=%.3f  (差=%.3f)\n', ...
            sir_list(si), auc_q1(si,2), auc_q1(si,3), diff_auc);
        found = true;
    end
end
if ~found
    fprintf('  （該当条件なし）\n');
end

%% プロット
lstyles = {'-','--','-.'};
colors  = lines(3);

fig = figure('Visible','off','Color','w','Position',[0 0 1400 450]);

% ─── 左: AUC vs SIR ───────────────────────────────────────────
subplot(1,3,1);
hold on; grid on; box on;
for mi = 1:3
    plot(sir_list, auc_q1(:,mi), lstyles{mi}, 'Color',colors(mi,:), ...
        'LineWidth',2.0,'Marker','o','MarkerSize',4,'DisplayName',mnames{mi});
end
yline(0.9,'k--','AUC=0.9','LabelHorizontalAlignment','right','FontSize',11,...
    'HandleVisibility','off');
xlabel('SIR [dB]  (ETC電力 / UAV干渉電力)','FontSize',13);
ylabel('AUC','FontSize',13);
title('AUC vs SIR（閾値非依存）','FontSize',14);
legend('Location','southeast','FontSize',11);
xlim([sir_list(1)-1  sir_list(end)+1]);
ylim([0.4 1.05]);

% ─── 中・右: ROC @ SIR=-5dB / -20dB ─────────────────────────
roc_titles = {sprintf('ROC @ SIR=%+ddB', roc_sirs(1)), ...
              sprintf('ROC @ SIR=%+ddB（最悪条件）', roc_sirs(2))};
for ri = 1:2
    subplot(1,3,ri+1);
    hold on; grid on; box on;
    if ~isempty(scores_roc{ri})
        for mi = 1:3
            [fpr_r, tpr_r, ~, auc_r] = perfcurve(labels_q1, scores_roc{ri}(:,mi), 1);
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
sgtitle(fig, 'SIR sweep AUC分析: Q1 Preamble vs Combined（閾値非依存）', ...
    'FontSize',12,'Color',[0.15 0.15 0.15]);
saveas(fig, fullfile(base_dir,'sir_auc_roc.png'));
fprintf('\n保存: experiments/sir_auc_roc.png\n');
disp('compare_sir_auc: OK');


%% ─── ローカル関数（compare_sir.m と同一実装） ──────────────────
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
