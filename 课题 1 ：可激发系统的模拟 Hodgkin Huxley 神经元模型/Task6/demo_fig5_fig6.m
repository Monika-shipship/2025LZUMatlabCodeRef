% demo_fig5_fig6.m — 复现"两个电耦合神经元的能量平衡"小节的图 5 & 图 6（单次长时间仿真）
% 运行前：确保同目录下已有 Task6TwoHH_EnergyCoupled.m

clear; clc;

% 基本 HH 参数
C = 1;
E.Na = 115; E.K = -12; E.L = 10.6;
g.Na = 120; g.K = 36; g.L = 0.3;

% 速率函数（把静息点平移到 0 mV 的常见写法）
Alpha.n = @(u) (0.1 - 0.01.*u) ./ (exp(1 - 0.1.*u) - 1);
Alpha.m = @(u) (2.5 - 0.1.*u) ./ (exp(2.5 - 0.1.*u) - 1);
Alpha.h = @(u) 0.07 .* exp(-u./20);
Beta.n  = @(u) 0.125 .* exp(-u./80);
Beta.m  = @(u) 4 .* exp(-u./18);
Beta.h  = @(u) 1 ./ (exp(3 - 0.1.*u) + 1);

% 仿真选项
% 论文做了很长时间平均（~750 s），先用 30 s 对齐形状，确认后再把 T 拉长
opt.fs      = 2e4;          % 20 kHz → dt=0.05 ms
opt.T       = 3.0e4;        % 30 s = 30000 ms
opt.T_warm  = 3.0e3;        % 丢弃前 3 s 热启动
opt.k_list  = 0:0.05:0.20;  % k 扫描

% 论文设定的随机电流（单位 μA/cm^2）：Pre ~ N(0,3^2)，Post ~ N(0,1^2)
opt.mu_pre     = 8.2;
opt.sigma_pre  = 3.0;
opt.mu_post    = 0.0;
opt.sigma_post = 1.0;

opt.seed    = 20251029;
opt.Ifdebug = 1;            % 打开命令行进度

% 计算
args = namedargs2cell(opt);
OUT = Task6TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,args{:});  % ← 注意这里是 {:}

% 图 5 — 离子通道代谢功率
figure('Name','图5 — 不同 k 下发送/接收神经元通道平均代谢能耗','Color','w');
plot(OUT.k, OUT.Pchan_pre_mean,  '-o','LineWidth',1.5); hold on;
plot(OUT.k, OUT.Pchan_post_mean, '-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('平均代谢能量消耗（通道总和，pJ/ms·cm^{-2}）');
title('图5.（彩色在线）k 变化下发送/接收神经元离子通道的平均代谢能耗（单向耦合）');
legend({'发送神经元（s）','接收神经元（r）'},'Location','northwest','Box','off');
grid on;

% 插图：平均放电频率
axes('Position',[0.60 0.30 0.32 0.32]);
plot(OUT.k, OUT.Fr_pre,  '-','LineWidth',1.2); hold on;
plot(OUT.k, OUT.Fr_post, '--','LineWidth',1.2);
xlabel('k'); ylabel('Hz'); title('平均放电频率');
legend({'s','r'},'Location','best','Box','off'); grid on;

% 图 6 — 突触处能量平衡两部分
figure('Name','图6 — 外介质供能与突触后部耗散（功率）','Color','w');
plot(OUT.k, OUT.P_amp_mean,       '-^','LineWidth',1.5); hold on;
plot(OUT.k, OUT.P_diss_post_mean, '-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('功率（nJ/s·cm^{-2} ≡ pJ/ms·cm^{-2}）');
title('图6.（彩色在线）从细胞外介质供能与突触后部位耗散（单向耦合）');
legend({'从外介质供能（放大器项）','突触后部位耗散'},'Location','northwest','Box','off');
grid on;

% 插图：连接处净能量导数（应为正）
axes('Position',[0.60 0.30 0.32 0.32]);
plot(OUT.k, OUT.P_net_mean,'-','LineWidth',1.2);
xlabel('k'); ylabel('pJ/ms·cm^{-2}'); 
% title('连接处净能量导数 ⟨V_r I_{couple}⟩'); 
title('连接处净能量导数 \langle P_{supply} + P_{dissip} \rangle'); 
grid on;
