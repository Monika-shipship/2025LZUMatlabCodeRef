% Demo_Fig5_Fig6_Parallel.m — 复现"两个电耦合神经元的能量平衡"的图 5 与图 6（并行版）
% 运行前：确保当前文件夹下已有 Task5TwoHH_EnergyCoupled.m
% 说明：
%   * 并行策略：以 k 为切片维度做 parfor；每个迭代装配一整列，再一次性写入 A(:,ik)
%   * 随机数：为每个段 s 与每个 k 组合设独立子流，保证可复现且与并行顺序无关
%   * 参数口径：对齐论文描述（突触前噪声均值 0，方差 9 => sigma=3；突触后均值 0，sigma=1）

clear; clc;

% 基本 HH 参数（与论文常见口径一致）
C = 1;
E.Na = 115; E.K = -12; E.L = 10.6;
g.Na = 120; g.K = 36; g.L = 0.3;

Alpha.n = @(u) (0.1 - 0.01 .* u) ./ (exp(1 - 0.1 .* u) - 1);
Alpha.m = @(u) (2.5 - 0.1 .* u) ./ (exp(2.5 - 0.1 .* u) - 1);
Alpha.h = @(u) 0.07 .* exp(-u ./ 20);
Beta.n  = @(u) 0.125 .* exp(-u ./ 80);
Beta.m  = @(u) 4 .* exp(-u ./ 18);
Beta.h  = @(u) 1 ./ (exp(3 - 0.1 .* u) + 1);

% 选项（等效总时长 = T * Nseg；先用 5 s × 12 段 = 60 s 调试；要更像论文可把 Nseg 加大）
opt.fs      = 4e4;            % 40 kHz（dt=0.025 ms）
opt.T       = 5e3;            % 每段 5 s（ms 记）
opt.Nseg    = 12;             % 段数；提升到 150 可近似 750 s
opt.T_warm  = 2e3;            % 每段前 2 s 丢弃
opt.k_list  = 0:0.01:0.20;    % k 扫描（mS/cm^2）

% 论文口径：突触前 N(0,3^2)，突触后 N(0,1^2)
opt.mu_pre     = 0.0;
opt.sigma_pre  = 3.0;
opt.mu_post    = 0.0;
opt.sigma_post = 1.0;

opt.seed    = 20251029;
opt.Ifdebug = 0;

% 打开并行池（若未开）
p = gcp('nocreate');
if isempty(p)
    p = parpool("Processes");
end
fprintf('Parallel pool with %d workers.\n', p.NumWorkers);

% 进度显示（DataQueue）
dq = parallel.pool.DataQueue;
done_col = 0;
afterEach(dq, @(msg) fprintf('%s\n', msg));

% 切片维度
S  = opt.Nseg;
K  = numel(opt.k_list);
fprintf('Total segments: %d; k points: %d\n', S, K);

% 结果容器（S × K），列切片写入
Pchan_pre_all  = zeros(S, K);
Pchan_post_all = zeros(S, K);
Psyn_post_all  = zeros(S, K);
Pjunc_all      = zeros(S, K);
Fr_pre_all     = zeros(S, K);
Fr_post_all    = zeros(S, K);

% 为可复现随机流做基数（不同列、不同段得到不同子流）
baseSeed = opt.seed;

% parfor 以 k 为外层切片维度；每次装配一整列
parfor ik = 1:K
    tStart = tic;

    % 本列的局部累加器
    col_Pchan_pre  = zeros(S,1);
    col_Pchan_post = zeros(S,1);
    col_Psyn_post  = zeros(S,1);
    col_Pjunc      = zeros(S,1);
    col_Fr_pre     = zeros(S,1);
    col_Fr_post    = zeros(S,1);

    % 逐段计算（串行，保证每列在同一 worker 内独立）
    for s = 1:S
        opt2 = opt;
        % 每段用独立种子（与 ik、s 绑定，顺序无关，可复现）
        opt2.seed = baseSeed + 10000*ik + s;

        % 调用单段函数（不带 Nseg）
        OUTs = Task5TwoHH_EnergyCoupled(C,E,g,Alpha,Beta, ...
            fs=opt2.fs, T=opt2.T, k_list=opt2.k_list(ik), ...
            mu_pre=opt2.mu_pre, sigma_pre=opt2.sigma_pre, ...
            mu_post=opt2.mu_post, sigma_post=opt2.sigma_post, ...
            seed=opt2.seed, T_warm=opt2.T_warm, Ifdebug=opt2.Ifdebug);

        col_Pchan_pre(s)  = OUTs.Pchan_pre_mean;
        col_Pchan_post(s) = OUTs.Pchan_post_mean;
        col_Psyn_post(s)  = OUTs.Psyn_post_mean;
        col_Pjunc(s)      = OUTs.Pjunc_mean;
        col_Fr_pre(s)     = OUTs.Fr_pre;
        col_Fr_post(s)    = OUTs.Fr_post;
    end

    % 一次性写入这列（满足 parfor "切片变量"规则）
    Pchan_pre_all(:,ik)  = col_Pchan_pre;
    Pchan_post_all(:,ik) = col_Pchan_post;
    Psyn_post_all(:,ik)  = col_Psyn_post;
    Pjunc_all(:,ik)      = col_Pjunc;
    Fr_pre_all(:,ik)     = col_Fr_pre;
    Fr_post_all(:,ik)    = col_Fr_post;

    send(dq, sprintf('  column %d/%d [k=%.3f] done in %.2fs', ...
        ik, K, opt.k_list(ik), toc(tStart)));
end

% 汇总平均（对段平均）
OUT.k               = opt.k_list;
OUT.Pchan_pre_mean  = mean(Pchan_pre_all,  1);
OUT.Pchan_post_mean = mean(Pchan_post_all, 1);
OUT.Psyn_post_mean  = mean(Psyn_post_all,  1);
OUT.Pjunc_mean      = mean(Pjunc_all,      1);
OUT.Fr_pre          = mean(Fr_pre_all,     1);
OUT.Fr_post         = mean(Fr_post_all,    1);

% 图 5
figure('Name','图5 — 不同 k 下发送/接收神经元通道平均代谢能耗','Color','w');
plot(OUT.k, OUT.Pchan_pre_mean, '-o','LineWidth',1.5); hold on;
plot(OUT.k, OUT.Pchan_post_mean,'-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('平均代谢能量消耗（通道总和，pJ/ms·cm^{-2}）');
title('图5.（彩色在线）k 变化下发送/接收神经元离子通道的平均代谢能耗（单向耦合）');
legend({'发送神经元（s）','接收神经元（r）'},'Location','northwest','Box','off'); grid on;

axes('Position',[0.54 0.30 0.32 0.32]);
plot(OUT.k, OUT.Fr_pre, '-','LineWidth',1.2); hold on;
plot(OUT.k, OUT.Fr_post,'--','LineWidth',1.2);
xlabel('k'); ylabel('Hz'); title('平均放电频率');
legend({'s','r'},'Location','best','Box','off'); grid on;

% 图 6
figure('Name','图6 — 外介质供能与突触后部耗散（功率）','Color','w');
plot(OUT.k, OUT.Psyn_post_mean, '-^','LineWidth',1.5); hold on;
plot(OUT.k, OUT.Pchan_post_mean,'-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('功率（nJ/s·cm^{-2} ≡ pJ/ms·cm^{-2}）');
title('图6.（彩色在线）细胞外介质供应（放大器项）与突触后膜通道耗散（单向耦合）');
legend({'从外介质供能（V_s I_{cpl}）','突触后膜通道耗散（Joule）'},'Location','northwest','Box','off');
grid on;

axes('Position',[0.54 0.30 0.32 0.32]);
plot(OUT.k, OUT.Pjunc_mean,'-','LineWidth',1.2);
xlabel('k'); ylabel('pJ/ms·cm^{-2}');
title('连接处净能量导数  \langle k(V_s - V_r)^2 \rangle'); grid on;

% 保存结果（可选）
jsonStr = jsonencode(OUT, 'PrettyPrint', true);
fid = fopen('ProDoubleOUT.txt','w','n','UTF-8'); fprintf(fid,'%s',jsonStr); fclose(fid);
disp('结构体已保存到: ProDoubleOUT.txt');
