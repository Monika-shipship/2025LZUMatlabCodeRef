% 任务 2：计算直流刺激强度和神经元发放频率之间的关系
clc; clear

tic
% 《Energy and information in Hodgkin-Huxley neurons》A. Moujahid, A. d'Anjou, and F. J. Torrealdea
% 中给出的参数是

% 电容：
C = 1;

% 各离子的能斯特电位和对应的电导：
E.Na = 115; E.K = -12; E.L = 10.6;
g.Na = 120; g.K = 36; g.L = 0.3;

% 对应的 alpha和beta函数：
Alpha.n = @(u) (0.1 - 0.01 .* u) ./ (exp(1 - 0.1 .* u) - 1);
Alpha.m = @(u) (2.5 - 0.1 .* u) ./ (exp(2.5 - 0.1 .* u) - 1);
Alpha.h = @(u) 0.07 .* exp(-u ./ 20);
Beta.n = @(u) 0.125 .* exp(-u ./ 80);
Beta.m = @(u) 4 .* exp(-u ./ 18);
Beta.h = @(u) 1 ./ (exp(3 - 0.1 .* u) + 1);

%%计算静息电位时的初值（在大规模运算时为了提高性能而这样设计）
fV = @(V, m, h, n, t) 1 / C * (-g.L * (V - E.L) - g.Na * m ^ 3 * h * (V - E.Na) - g.K * n ^ 4 * (V - E.K) + I_ext(t));
fm = @(V, m) Alpha.m(V) * (1 - m) - Beta.m(V) * m;
fh = @(V, h) Alpha.h(V) * (1 - h) - Beta.h(V) * h;
fn = @(V, n) Alpha.n(V) * (1 - n) - Beta.n(V) * n;
% 为了求解微分方程，我们还必须获得初始值，在静息时，I_ext肯定是0，V,m,h,n的导数也是零
% 由此可以解方程数值求得初始值，同时也是静息值
% 使用线性方程解出静息时的电位，各离子浓度比，此时I_ext肯定是0，所以得重新命名一个函数fVT来解方程

syms V0 m0 h0 n0
fVT = @(V, m, h, n) 1 / C * (-g.L * (V - E.L) - g.Na * m ^ 3 * h * (V - E.Na) - g.K * n ^ 4 * (V - E.K) + 0);
eqt = [fVT(V0, m0, h0, n0), fm(V0, m0) == 0, fh(V0, h0) == 0, fn(V0, n0) == 0];
VPASol = vpasolve(eqt, [V0 m0 h0 n0]);

% fs= 6000;%采样频率
T = 1000; %总时长(ms)
fs = 13 * T; %采样频率随时间增长而增加，这样可以有效避免解爆炸
t_vec = 0:1 / fs * 1000:T; % 创建时间向量
I_template = (0.1 * T <= t_vec & t_vec <= 0.9 * T); % 预生成电流模板

DelA = 0.02;
Amax = 50;
Avec = 0:DelA:Amax;
N_A = length(Avec);
fPkz = zeros(1, N_A);

toc
Count = 1;

% parfor i = 1:N_A %刺激电流强度

%     A = Avec(i);
%     I_ext_vec = A * I_template; %脉冲刺激向量

%     [V, m, h, n] = Task2I2V_fast(C, E, g, Alpha, Beta, VPASol, ...
%         'fs', fs, ...
%         'T', T, ...
%         'A', A, ...
%         'I_ext', I_ext_vec ...
%     );
%     Threold = max(V) * 2/3;
%     % findpeaks(V.*(V>=Threold));
%     [pks, locs] = findpeaks(V .* (V >= Threold));
%     NPkz = length(pks);
%     fPkz(i) = NPkz / (0.8 * T * 1e-3);
%     % disp(fPkz)
% end

fPkz = parfor_Task2I2V_fast(N_A,Avec,I_template, C, E, g, Alpha, Beta, VPASol, fs, T);

figure;
plot(Avec, fPkz)
xlabel('I_{ext} (\muA/cm^2)'), ylabel('f/Hz')
title('直流刺激强度和神经元发放频率之间的关系 I_{ext}-t')
toc

function fPkz = parfor_Task2I2V_fast(N_A,Avec,I_template, C, E, g, Alpha, Beta, VPASol, fs, T)

    %  启动并行池,没有则自动开
    p = gcp('nocreate'); if isempty(p), parpool; end

    %  预分配结果
    fPkz = zeros(1, N_A);

    % 建立数据通道 并 回调（从 worker 发消息到客户端打印进度）
    dq = parallel.pool.DataQueue;
    done = 0;
    afterEach(dq, @onUpdate);

    % 并行循环,do , 报告 完成1次
    parfor i = 1:N_A
        fprintf('%d / %d \n',i,N_A);
        A = Avec(i);
        I_ext_vec = A * I_template; %脉冲刺激向量

        [V, m, h, n] = Task2I2V_fast(C, E, g, Alpha, Beta, VPASol, ...
            'fs', fs, ...
            'T', T, ...
            'A', A, ...
            'I_ext', I_ext_vec ...
        );
        Threold = max(V) * 2/3;
        % findpeaks(V.*(V>=Threold));
        [pks, locs] = findpeaks(V .* (V >= Threold));
        NPkz = length(pks);
        fPkz(i) = NPkz / (0.8 * T * 1e-3);
        % disp(fPkz)
        send(dq, 1); % 告诉客户端又完成了 1 次
    end

    % 换行收尾
    fprintf('\n完成。\n');

    % 客户端回调 每收到一次 send(dq,1) 就更新一回
    function onUpdate(~)
        done = done + 1;
        % 单行覆盖输出 已完成/总数
        fprintf('\r已完成：%d/%d (%.1f%%)', done, N_A, 100 * done / N_A);
    end

end
