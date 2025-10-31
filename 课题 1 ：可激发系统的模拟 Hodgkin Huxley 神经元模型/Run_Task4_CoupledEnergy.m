% 拓展任务：两个电耦合神经元的能量平衡
%
% 目标：
% 1. 模拟两个单向电耦合的 H-H 神经元。
% 2. 神经元1 (发送) 接收均值为 4.5, 方差为 16 的高斯噪声电流。
% 3. 神经元2 (接收) 接收均值为 0,   方差为 16 的高斯噪声电流。
% 4. 仿真不同突触电导 k 值下的系统。
% 5. 绘制 图5：平均代谢能量消耗 vs k (插图：平均发放频率 vs k)
% 6. 绘制 图6：突触后能量供应/耗散 vs k (插图：净能量导数 vs k)

clc; clear; close all;
tic
fprintf('开始执行拓展任务：两个电耦合神经元的能量平衡...\n');

%% 1. 定义 H-H 模型基础参数
% (参数来源于 HHtask2.m，遵循你的风格)

% 电容：
C = 1; % (uF/cm^2)

% 各离子的能斯特电位和对应的电导：
E.Na = 115; E.K = -12; E.L = 10.6; % (mV)
g.Na = 120; g.K = 36; g.L = 0.3;  % (mS/cm^2)

% 对应的 alpha和beta函数： (u = V - V_rest = V - (-65) = V + 65)
% *注意*：你之前的 Alpha/Beta 函数是基于 V_rest=0mV 的。
%         而论文中 V_rest = -65mV。这里的 V 已经是 V-V_rest
%         为了与你的 HHtask2.m 保持一致，我们继续使用 V_rest=0 的定义。
Alpha.n = @(u) (0.1 - 0.01 .* u) ./ (exp(1 - 0.1 .* u) - 1);
Alpha.m = @(u) (2.5 - 0.1 .* u) ./ (exp(2.5 - 0.1 .* u) - 1);
Alpha.h = @(u) 0.07 .* exp(-u ./ 20);
Beta.n = @(u) 0.125 .* exp(-u ./ 80);
Beta.m = @(u) 4 .* exp(-u ./ 18);
Beta.h = @(u) 1 ./ (exp(3 - 0.1 .* u) + 1);

% 计算静息电位时的初值（来自 HHtask2.m）
% 我们假设静息电位为 V_rest = 0 (因为你的alpha/beta函数是这样定义的)
% 求解 fV = 0, fm = 0, fh = 0, fn = 0
syms V0 m0 h0 n0
eqt = [ ...
    (1/C) * (g.Na * m0^3 * h0 * (V0 - E.Na) + g.K * n0^4 * (V0 - E.K) + g.L * (V0 - E.L)) == 0, ...
    Alpha.m(V0) * (1 - m0) - Beta.m(V0) * m0 == 0, ...
    Alpha.h(V0) * (1 - h0) - Beta.h(V0) * h0 == 0, ...
    Alpha.n(V0) * (1 - n0) - Beta.n(V0) * n0 == 0 ...
    ];
VPASol = vpasolve(eqt, [sym('V0'), sym('m0'), sym('h0'), sym('n0')]);
if isempty(VPASol.V0)
    VPASol.V0=0; VPASol.m0=0.05; VPASol.h0=0.6; VPASol.n0=0.3;
    warning('未找到精确静息解，使用近似初值。');
else
    % 确保 VPASol 是 double
    VPASol.V0 = double(VPASol.V0(1));
    VPASol.m0 = double(VPASol.m0(1));
    VPASol.h0 = double(VPASol.h0(1));
    VPASol.n0 = double(VPASol.n0(1));
end
fprintf('静息电位初值：V=%.2f, m=%.2f, h=%.2f, n=%.2f\n', VPASol.V0, VPASol.m0, VPASol.h0, VPASol.n0);


%% 2. 定义仿真参数
fs = 5000;       % 采样频率 (Hz)
T = 2000;        % 仿真总时长 (ms)，需要足够长的时间来平均噪声
Del = 1/fs * 1000; % 时间步长 (ms)，因为H-H方程单位是ms
t_vec = 0:Del:T;   % 时间向量 (ms)
N = length(t_vec); % 总步数

% 突触电导 k
k_vec = linspace(0, 1.0, 51); % 突触电导 k 的范围 (0 到 1.0, 51个点)
N_k = length(k_vec);

%% 3. 定义激励参数 (高斯噪声)
I1_mean = 4.5;   % 神经元1 均值 (uA/cm^2)
I1_var = 16;     % 神经元1 方差
I1_std = sqrt(I1_var);

I2_mean = 0;     % 神经元2 均值 (uA/cm^2)
I2_var = 16;     % 神经元2 方差
I2_std = sqrt(I2_var);

% 预先生成噪声电流向量 (与 parfor_Task4_solver 中的 RK4 兼容)
% 注意：我们需要 (N) 长度的向量。
% 为了在 parfor 中保持一致性，我们在循环外生成一个"种子"或"基础"噪声
% 更好的做法是在 parfor 内部为每次 k 独立生成噪声
% （因为每次 k 都是一次独立的仿真实验）

fprintf('为 %d 个 k 值启动并行仿真...\n', N_k);

% 预分配结果向量
avg_E_meta1 = zeros(1, N_k);
avg_E_meta2 = zeros(1, N_k);
avg_freq1 = zeros(1, N_k);
avg_freq2 = zeros(1, N_k);
avg_E_supply = zeros(1, N_k);
avg_E_dissip = zeros(1, N_k);
avg_E_net = zeros(1, N_k);

% 启动并行池, 模仿 HHtask2.m
p = gcp('nocreate'); if isempty(p), parpool; end

% 并行循环
parfor i = 1:N_k
    k_val = k_vec(i);
    
    % 在 parfor 内部生成独立的噪声流
    % randn('seed', i); % (旧版)
    rng(i); % (新版) 确保每次运行的可复现性和独立性
    
    % 生成高斯白噪声。注意：理论上应使用 sqrt(Del) 缩放，
    % 但这里我们建模为具有特定方差的 *时间序列* 电流，
    % 这更符合"高斯噪声"的描述。
    I_noise1 = I1_mean + I1_std * randn(1, N);
    I_noise2 = I2_mean + I2_std * randn(1, N);

    % 调用求解器
    [E_m1, E_m2, f1, f2, E_s, E_d] = parfor_Task4_solver( ...
        k_val, C, E, g, Alpha, Beta, VPASol, ...
        I_noise1, I_noise2, fs, T ...
        );
    
    % 存储结果
    avg_E_meta1(i) = E_m1;
    avg_E_meta2(i) = E_m2;
    avg_freq1(i) = f1;
    avg_freq2(i) = f2;
    avg_E_supply(i) = E_s;
    avg_E_dissip(i) = E_d;
    
    fprintf('完成 k = %.3f (%d / %d)\n', k_val, i, N_k);
end

% 净能量
avg_E_net = avg_E_supply - avg_E_dissip;

fprintf('仿真完成。\n');
toc

%% 4. 绘制 图5 (代谢能量消耗)

% 创建保存目录
outputDir = 'FigOutput/拓展';
if ~exist(outputDir, 'dir')
   mkdir(outputDir);
end

hFig5 = figure('Name', '图5');
hold on;
plot(k_vec, avg_E_meta1, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(k_vec, avg_E_meta2, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
hold off;
xlabel('突触电导 $k$ $(mS/cm^2)$', 'Interpreter', 'latex');
ylabel('平均代谢能量消耗 $(mW/cm^2)$', 'Interpreter', 'latex'); % 能量单位是 (uA/cm^2) * (mV) = (uW/cm^2)。 1000 uW = 1 mW
title('图5: 发送和接收神经元离子通道的平均代谢能量消耗');
legend('发送神经元 ($V_1$)', '接收神经元 ($V_2$)', 'Interpreter', 'latex');
grid on;
box on;

% 插入图 (Inset)
ax_main = gca;
% 根据主图的位置创建插图
inset_pos = [ax_main.Position(1)+0.55*ax_main.Position(3), ...
             ax_main.Position(2)+0.15*ax_main.Position(4), ...
             0.33 * ax_main.Position(3), ...
             0.33 * ax_main.Position(4)];
axes('Position', inset_pos);
hold on;
plot(k_vec, avg_freq1, 'b--');
plot(k_vec, avg_freq2, 'r--');
hold off;
xlabel('$k$', 'Interpreter', 'latex');
ylabel('平均发放频率 (Hz)');
grid on;
box on;

% 保存图像
saveas(hFig5, fullfile(outputDir, '图5_k值与平均能量消耗.fig'));
saveas(hFig5, fullfile(outputDir, '图5_k值与平均能量消耗.pdf'));
fprintf('图5 已保存到 %s\n', outputDir);

%% 5. 绘制 图6 (突触能量平衡)

hFig6 = figure('Name', '图6');
hold on;
plot(k_vec, avg_E_supply, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(k_vec, avg_E_dissip, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
hold off;
xlabel('突触电导 $k$ $(mS/cm^2)$', 'Interpreter', 'latex');
ylabel('平均能量速率 $(mW/cm^2)$', 'Interpreter', 'latex');
title('图6: 突触后部位的平均能量供应和耗散');
legend('供应 (从细胞外)', '耗散 (向细胞外)');
grid on;
box on;

% 插入图 (Inset)
ax_main2 = gca;
inset_pos2 = [ax_main2.Position(1)+0.15*ax_main2.Position(3), ...
              ax_main2.Position(2)+0.55*ax_main2.Position(4), ...
              0.33 * ax_main2.Position(3), ...
              0.33 * ax_main2.Position(4)];
axes('Position', inset_pos2);
plot(k_vec, avg_E_net, 'g-');
xlabel('$k$', 'Interpreter', 'latex');
ylabel('净能量导数');
grid on;
box on;
ax_inset2 = gca;
ax_inset2.YAxis.Exponent = floor(log10(abs(max(avg_E_net(avg_E_net>0))))); % 自动调整数量级

% 保存图像
saveas(hFig6, fullfile(outputDir, '图6_k值与突触能量平衡.fig'));
saveas(hFig6, fullfile(outputDir, '图6_k值与突触能量平衡.pdf'));
fprintf('图6 已保存到 %s\n', outputDir);

fprintf('所有任务完成。\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  并行计算的 RK4 求解器函数
%  (风格模仿自 HHtask2.m 中的 parfor_Task2I2V_fast)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avg_E_meta1, avg_E_meta2, avg_freq1, avg_freq2, avg_E_supply, avg_E_dissip] ...
    = parfor_Task4_solver(k, C, E, g, Alpha, Beta, VPASol, ...
    I_noise1, I_noise2, fs, T)
    
    % --- 1. 初始化 ---
    Del = 1/fs * 1000; % (ms)
    N = length(I_noise1);

    % 初始化两个神经元的 V, m, h, n 状态向量
    V1 = zeros(1, N); m1 = zeros(1, N); h1 = zeros(1, N); n1 = zeros(1, N);
    V2 = zeros(1, N); m2 = zeros(1, N); h2 = zeros(1, N); n2 = zeros(1, N);

    % 初始化能量(功率)向量
    P_meta1 = zeros(1, N); % 神经元1 的瞬时代谢功率
    P_meta2 = zeros(1, N); % 神经元2 的瞬时代谢功率
    P_supply = zeros(1, N); % 突触供应给 2 的瞬时功率
    P_dissip = zeros(1, N); % 突触从 2 耗散的瞬时功率

    % 设置静息初值
    V1(1) = VPASol.V0; m1(1) = VPASol.m0; h1(1) = VPASol.h0; n1(1) = VPASol.n0;
    V2(1) = VPASol.V0; m2(1) = VPASol.m0; h2(1) = VPASol.h0; n2(1) = VPASol.n0;

    % --- 2. 定义8个微分方程 (f(V, m, h, n)) ---
    
    % 神经元1
    fV1 = @(V1, m1, h1, n1, V2, I_ext, k_syn) (1/C) * (I_ext ...
        - g.Na * m1^3 * h1 * (V1 - E.Na) ...
        - g.K * n1^4 * (V1 - E.K) ...
        - g.L * (V1 - E.L) ...
        - k_syn * (V1 - V2)); % 单向耦合：电流从 1 流向 2
    fm1 = @(V, m) Alpha.m(V) * (1 - m) - Beta.m(V) * m;
    fh1 = @(V, h) Alpha.h(V) * (1 - h) - Beta.h(V) * h;
    fn1 = @(V, n) Alpha.n(V) * (1 - n) - Beta.n(V) * n;

    % 神经元2
    fV2 = @(V2, m2, h2, n2, V1, I_ext, k_syn) (1/C) * (I_ext ...
        - g.Na * m2^3 * h2 * (V2 - E.Na) ...
        - g.K * n2^4 * (V2 - E.K) ...
        - g.L * (V2 - E.L) ...
        + k_syn * (V1 - V2)); % 单向耦合：电流从 1 流向 2
    fm2 = @(V, m) Alpha.m(V) * (1 - m) - Beta.m(V) * m;
    fh2 = @(V, h) Alpha.h(V) * (1 - h) - Beta.h(V) * h;
    fn2 = @(V, n) Alpha.n(V) * (1 - n) - Beta.n(V) * n;

    % --- 3. RK4 求解循环 ---
    for i = 1:N-1
        % 计算 t(i) 时的能量和功率 (基于 V(i), m(i)...)
        % (注意：能量应在更新 V(i+1) 之前计算，或之后计算，保持一致即可)
        % 我们在更新 *之后* 计算（即在 i=1 时计算 V(1) 的能量）
        
        I_Na1 = g.Na * m1(i)^3 * h1(i) * (V1(i) - E.Na);
        I_K1 = g.K * n1(i)^4 * (V1(i) - E.K);
        I_L1 = g.L * (V1(i) - E.L);
        
        I_Na2 = g.Na * m2(i)^3 * h2(i) * (V2(i) - E.Na);
        I_K2 = g.K * n2(i)^4 * (V2(i) - E.K);
        I_L2 = g.L * (V2(i) - E.L);
        
        % 图5：代谢能量 (功率)
        P_meta1(i) = I_Na1 * (V1(i) - E.Na) + I_K1 * (V1(i) - E.K) + I_L1 * (V1(i) - E.L);
        P_meta2(i) = I_Na2 * (V2(i) - E.Na) + I_K2 * (V2(i) - E.K) + I_L2 * (V2(i) - E.L);

        % 图6：突触能量 (功率)
        I_syn = k * (V1(i) - V2(i));
        P_syn_V2 = I_syn * V2(i); % 耦合电流在 V2 上的功率
        
        P_supply(i) = max(P_syn_V2, 0); % 供应
        P_dissip(i) = max(-P_syn_V2, 0); % 耗散
        

        % --- RK4 K1 ---
        I1_t1 = I_noise1(i);
        I2_t1 = I_noise2(i);
        
        k1_V1 = fV1(V1(i), m1(i), h1(i), n1(i), V2(i), I1_t1, k);
        k1_m1 = fm1(V1(i), m1(i));
        k1_h1 = fh1(V1(i), h1(i));
        k1_n1 = fn1(V1(i), n1(i));
        
        k1_V2 = fV2(V2(i), m2(i), h2(i), n2(i), V1(i), I2_t1, k);
        k1_m2 = fm2(V2(i), m2(i));
        k1_h2 = fh2(V2(i), h2(i));
        k1_n2 = fn2(V2(i), n2(i));
        
        % --- RK4 K2 ---
        I1_t2 = (I_noise1(i) + I_noise1(i+1)) / 2;
        I2_t2 = (I_noise2(i) + I_noise2(i+1)) / 2;
        
        k2_V1 = fV1(V1(i) + Del/2 * k1_V1, m1(i) + Del/2 * k1_m1, h1(i) + Del/2 * k1_h1, n1(i) + Del/2 * k1_n1, V2(i) + Del/2 * k1_V2, I1_t2, k);
        k2_m1 = fm1(V1(i) + Del/2 * k1_V1, m1(i) + Del/2 * k1_m1);
        k2_h1 = fh1(V1(i) + Del/2 * k1_V1, h1(i) + Del/2 * k1_h1);
        k2_n1 = fn1(V1(i) + Del/2 * k1_V1, n1(i) + Del/2 * k1_n1);
        
        k2_V2 = fV2(V2(i) + Del/2 * k1_V2, m2(i) + Del/2 * k1_m2, h2(i) + Del/2 * k1_h2, n2(i) + Del/2 * k1_n2, V1(i) + Del/2 * k1_V1, I2_t2, k);
        k2_m2 = fm2(V2(i) + Del/2 * k1_V2, m2(i) + Del/2 * k1_m2);
        k2_h2 = fh2(V2(i) + Del/2 * k1_V2, h2(i) + Del/2 * k1_h2);
        k2_n2 = fn2(V2(i) + Del/2 * k1_V2, n2(i) + Del/2 * k1_n2);
        
        % --- RK4 K3 ---
        % (I1_t2, I2_t2 仍然是中点电流)
        
        k3_V1 = fV1(V1(i) + Del/2 * k2_V1, m1(i) + Del/2 * k2_m1, h1(i) + Del/2 * k2_h1, n1(i) + Del/2 * k2_n1, V2(i) + Del/2 * k2_V2, I1_t2, k);
        k3_m1 = fm1(V1(i) + Del/2 * k2_V1, m1(i) + Del/2 * k2_m1);
        k3_h1 = fh1(V1(i) + Del/2 * k2_V1, h1(i) + Del/2 * k2_h1);
        k3_n1 = fn1(V1(i) + Del/2 * k2_V1, n1(i) + Del/2 * k2_n1);
        
        k3_V2 = fV2(V2(i) + Del/2 * k2_V2, m2(i) + Del/2 * k2_m2, h2(i) + Del/2 * k2_h2, n2(i) + Del/2 * k2_n2, V1(i) + Del/2 * k2_V1, I2_t2, k);
        k3_m2 = fm2(V2(i) + Del/2 * k2_V2, m2(i) + Del/2 * k2_m2);
        k3_h2 = fh2(V2(i) + Del/2 * k2_V2, h2(i) + Del/2 * k2_h2);
        k3_n2 = fn2(V2(i) + Del/2 * k2_V2, n2(i) + Del/2 * k2_n2);

        % --- RK4 K4 ---
        I1_t3 = I_noise1(i+1);
        I2_t3 = I_noise2(i+1);
        
        k4_V1 = fV1(V1(i) + Del * k3_V1, m1(i) + Del * k3_m1, h1(i) + Del * k3_h1, n1(i) + Del * k3_n1, V2(i) + Del * k3_V2, I1_t3, k);
        k4_m1 = fm1(V1(i) + Del * k3_V1, m1(i) + Del * k3_m1);
        k4_h1 = fh1(V1(i) + Del * k3_V1, h1(i) + Del * k3_h1);
        k4_n1 = fn1(V1(i) + Del * k3_V1, n1(i) + Del * k3_n1);
        
        k4_V2 = fV2(V2(i) + Del * k3_V2, m2(i) + Del * k3_m2, h2(i) + Del * k3_h2, n2(i) + Del * k3_n2, V1(i) + Del * k3_V1, I2_t3, k);
        k4_m2 = fm2(V2(i) + Del * k3_V2, m2(i) + Del * k3_m2);
        k4_h2 = fh2(V2(i) + Del * k3_V2, h2(i) + Del * k3_h2);
        k4_n2 = fn2(V2(i) + Del * k3_V2, n2(i) + Del * k3_n2);
        
        % --- 更新 V, m, h, n ---
        V1(i+1) = V1(i) + (Del/6) * (k1_V1 + 2*k2_V1 + 2*k3_V1 + k4_V1);
        m1(i+1) = m1(i) + (Del/6) * (k1_m1 + 2*k2_m1 + 2*k3_m1 + k4_m1);
        h1(i+1) = h1(i) + (Del/6) * (k1_h1 + 2*k2_h1 + 2*k3_h1 + k4_h1);
        n1(i+1) = n1(i) + (Del/6) * (k1_n1 + 2*k2_n1 + 2*k3_n1 + k4_n1);
        
        V2(i+1) = V2(i) + (Del/6) * (k1_V2 + 2*k2_V2 + 2*k3_V2 + k4_V2);
        m2(i+1) = m2(i) + (Del/6) * (k1_m2 + 2*k2_m2 + 2*k3_m2 + k4_m2);
        h2(i+1) = h2(i) + (Del/6) * (k1_h2 + 2*k2_h2 + 2*k3_h2 + k4_h2);
        n2(i+1) = n2(i) + (Del/6) * (k1_n2 + 2*k2_n2 + 2*k3_n2 + k4_n2);
    end

    % --- 4. 计算平均值 ---
    
    % 丢弃前 20% 的数据 (例如 20% * T)，避免瞬态效应
    transient_steps = floor(N * 0.2);
    
    % 计算平均功率 (uA/cm^2 * mV = uW/cm^2)
    % 转换为 mW/cm^2
    avg_E_meta1 = mean(P_meta1(transient_steps:end)) / 1000;
    avg_E_meta2 = mean(P_meta2(transient_steps:end)) / 1000;
    avg_E_supply = mean(P_supply(transient_steps:end)) / 1000;
    avg_E_dissip = mean(P_dissip(transient_steps:end)) / 1000;
    
    % 计算发放频率
    T_stable = T * (1 - 0.2); % 稳定期时长 (ms)
    
    [pks1, ~] = findpeaks(V1(transient_steps:end), 'MinPeakHeight', 20); % 阈值设为 20mV
    [pks2, ~] = findpeaks(V2(transient_steps:end), 'MinPeakHeight', 20);
    
    avg_freq1 = length(pks1) / (T_stable * 1e-3); % 频率 (Hz)
    avg_freq2 = length(pks2) / (T_stable * 1e-3);
    
end