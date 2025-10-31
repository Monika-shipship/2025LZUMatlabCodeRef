% demo_fig5_fig6_parallel.m — (并行加速)
% 注意需要和 run_HH_coupling_single_k.m放在同一文件夹中
clear; clc;  
% 基本 HH 参数
C = 1;
E.Na = 115; E.K = -12; E.L = 10.6;
g.Na = 120; g.K = 36; g.L = 0.3;
% Alpha/Beta 不再需要，因为它们已内联到函数中

% 仿真选项 (打包到 opt 结构体中)
opt.fs      = 2e4;          % 20 kHz → dt=0.05 ms
opt.T       = 3.0e4;        % 30 s = 30000 ms
opt.T_warm  = 3.0e3;        % 丢弃前 3 s 热启动
opt.mu_pre     = 8.4;       % 调整后的值，原论文这里设置的是0，显然不合理
opt.sigma_pre  = 3.0;
opt.mu_post    = 0.0;
opt.sigma_post = 1.0;
opt.seed    = 20251029;
opt.Ifdebug = 1;            % 调试开关

k_list  = 0:0.01:0.20;  % k 扫描
K = numel(k_list);
fprintf('开始并行计算 %d 个 k 值 (T=%.1f s)...\n', K, opt.T/1000);

% --- 并行计算 ---
% (已修复) 使用元胞数组 (cell array) 来收集不同结构体
OUT_cell = cell(1, K); 
parfor ik = 1:K
    % (已修复) 不再传递 Alpha 和 Beta
    OUT_cell{ik} = run_HH_coupling_single_k(k_list(ik), C,E,g,opt);
end

fprintf('计算完成，正在重组数据...\n');

% --- 数据重组 ---
% (已修复) 将元胞数组转换为结构体数组，然后再解包
OUT_array = [OUT_cell{:}]; 

OUT.k                 = [OUT_array.k];
OUT.Pchan_pre_mean    = [OUT_array.Pchan_pre_mean];
OUT.Pchan_post_mean   = [OUT_array.Pchan_post_mean];
OUT.Fr_pre            = [OUT_array.Fr_pre];
OUT.Fr_post           = [OUT_array.Fr_post];
OUT.P_amp_mean        = [OUT_array.P_amp_mean];
OUT.P_diss_post_mean  = [OUT_array.P_diss_post_mean];
OUT.P_net_mean        = [OUT_array.P_net_mean];

%% --- 绘图 (与您的代码完全相同) ---
% 图 5
figure('Name','图5 — 不同 k 下发送/接收神经元通道平均代谢能耗','Color','w');
plot(OUT.k, OUT.Pchan_pre_mean,  '-o','LineWidth',1.5); hold on;
plot(OUT.k, OUT.Pchan_post_mean, '-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('平均代谢能量消耗（通道总和，pJ/ms·cm^{-2}）');
title('图5.（彩色在线）k 变化下发送/接收神经元离子通道的平均代谢能耗（单向耦合）');
legend({'发送神经元（s）','接收神经元（r）'},'Location','northwest','Box','off');
grid on;
axes('Position',[0.60 0.30 0.32 0.32]);
plot(OUT.k, OUT.Fr_pre,  '-o','LineWidth',1.2); hold on;
plot(OUT.k, OUT.Fr_post, '--^','LineWidth',1.2);
xlabel('k'); ylabel('Hz'); title('平均放电频率');
legend({'s','r'},'Location','best','Box','off'); grid on;
% 图 6
figure('Name','图6 — 外介质供能与突触后部耗散（功率）','Color','w');
plot(OUT.k, OUT.P_amp_mean,       '-^','LineWidth',1.5); hold on;
plot(OUT.k, OUT.P_diss_post_mean, '-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('功率（nJ/s·cm^{-2} ≡ pJ/ms·cm^{-2}）');
title('图6.（彩色在线）从细胞外介质供能与突触后部位耗散（单向耦合）');
legend({'从外介质供能（放大器项）','突触后部位耗散'},'Location','northwest','Box','off');
grid on;
axes('Position',[0.60 0.30 0.32 0.32]);
plot(OUT.k, OUT.P_net_mean,'-o','LineWidth',1.2);
xlabel('k'); ylabel('pJ/ms·cm^{-2}'); 
title('连接处净能量导数 \langle P_{supply} + P_{dissip} \rangle'); 
grid on;

% figure('Name',' — 不同 k 下发送/接收神经元通道平均代谢能耗','Color','w');
% plot(OUT.k, OUT.Pchan_pre_mean,  '-','LineWidth',1.5); hold on;
% plot(OUT.k, OUT.Pchan_post_mean, '-','LineWidth',1.5);
% xlabel('突触电导  k  (mS/cm^2)');
% ylabel('平均代谢能量消耗（通道总和，pJ/ms·cm^{-2}）');
% title('k 变化下发送/接收神经元离子通道的平均代谢能耗（单向耦合）');
% legend({'发送神经元（s）','接收神经元（r）'},'Location','northwest','Box','off');
% grid on;
% axes('Position',[0.60 0.30 0.32 0.32]);
% plot(OUT.k, OUT.Fr_pre,  '-','LineWidth',1.2); hold on;
% plot(OUT.k, OUT.Fr_post, '-','LineWidth',1.2);
% xlabel('k'); ylabel('Hz'); title('平均放电频率');
% legend({'s','r'},'Location','best','Box','off'); grid on;
% % 图 6
% figure('Name',' — 外介质供能与突触后部耗散（功率）','Color','w');
% plot(OUT.k, OUT.P_amp_mean,       '-','LineWidth',1.5); hold on;
% plot(OUT.k, OUT.P_diss_post_mean, '-','LineWidth',1.5);
% xlabel('突触电导  k  (mS/cm^2)');
% ylabel('功率（nJ/s·cm^{-2} ≡ pJ/ms·cm^{-2}）');
% title('从细胞外介质供能与突触后部位耗散（单向耦合）');
% legend({'从外介质供能（放大器项）','突触后部位耗散'},'Location','northwest','Box','off');
% grid on;
% axes('Position',[0.60 0.30 0.32 0.32]);
% plot(OUT.k, OUT.P_net_mean,'-','LineWidth',1.2);
% xlabel('k'); ylabel('pJ/ms·cm^{-2}'); 
% title('连接处净能量导数 \langle P_{supply} + P_{dissip} \rangle'); 
% grid on;