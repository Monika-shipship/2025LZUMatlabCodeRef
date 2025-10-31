
% demo_fig5_fig6.m — 复现“两个电耦合神经元的能量平衡”小节的图 5 & 图 6
% 运行前：确保本文件夹下已有 Task5TwoHH_EnergyCoupled.m，且主文件 HodgkinHuxley.m 中的 Alpha/Beta/E/g 一致。

clear; clc;

%  基本 HH 参数
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

% 选项 
opt.fs      = 5e4;           % 50 kHz (dt = 0.02 ms)
opt.T       = 2e4;          %  20000 ms 每段 20 s
opt.Nseg    = 6;             % 分 6 段累计平均  等效 60 s
opt.T_warm  = 5e3;           % 分 6 段，换种子重复 ， 等效 120 s 平均
opt.k_list  = 0:0.05:0.20;   % 突触电导 k 值（mS/cm^2） % 0~0.20 mS/cm^2，步长 0.01

% 外部"噪声"电流（高斯）：
opt.mu_pre     = 5.0;        % 发送端，前神经元均值 推动其稳定放电 让前端平均放电 ≈ 45 Hz（±2Hz ）
opt.sigma_pre  = 1.0;
opt.mu_post    = 0.0;        % 接收端总噪声均值（中间神经元，默认 0） 中间神经元 仅受总噪声 + 耦合
opt.sigma_post = 1.0;

opt.seed    = 20251029;
opt.Ifdebug = 0;

%  计算 
% args = namedargs2cell(opt);
% OUT = Task5TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,args{:});

%% 多段等效长时平均（换随机种子叠加）
% S = opt.Nseg; 
% acc_fields = {'Pchan_pre_mean','Pchan_post_mean','Psyn_post_mean','Pjunc_mean','Fr_pre','Fr_post'};
% OUT = [];% 汇总输出
% for s = 1:S
%     opt2 = opt; 
%     opt2.seed = opt.seed + 97*s;% 每段不同随机种子，改变随机种子以减少统计方差
%     args = namedargs2cell(opt2); % 函数本体不需要 Nseg，将结构体转为 name-value 形参
%     OUTs = Task5TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,args{:});
%     if s==1
%         OUT = OUTs;
%         for f = 1:numel(acc_fields)
%             OUT.(acc_fields{f}) = zeros(size(OUTs.(acc_fields{f})));
%         end
%     end
%     for f = 1:numel(acc_fields)
%         OUT.(acc_fields{f}) = OUT.(acc_fields{f}) + OUTs.(acc_fields{f});
%     end
% end
% for f = 1:numel(acc_fields)
%     OUT.(acc_fields{f}) = OUT.(acc_fields{f}) / S; % 分段平均
% end

%% ===== 计算（含：多段平均 + 命令行进度）=====
S = opt.Nseg;
fields = {'Pchan_pre_mean','Pchan_post_mean','Psyn_post_mean','Pjunc_mean','Fr_pre','Fr_post'};

OUT_accum = [];
fprintf('Total segments: %d\n', S);
tStartAll = tic;

for s = 1:S
    tStartSeg = tic;

    % 每段使用不同随机种子；Task5 接受 Nseg 字段，不影响
    opt2 = opt;  opt2.seed = opt.seed + 97*s;

    % —— 段内计算 ——（含 k 扫描的命令行进度）
    K = numel(opt2.k_list);
    fprintf('  Segment %d/%d | k-sweep (%d points): ', s, S, K);

    args  = namedargs2cell(opt2);
    OUTs  = Task5TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,args{:});    % 真正计算在函数里

    % 在函数内部已完成 k 扫描；这里只负责显示一次"完成"
    fprintf('[done]  (%.1fs)\n', toc(tStartSeg));

    % —— 累计平均 —— 
    if s==1
        OUT = OUTs;  % 复制结构（带 k 轴）
        for f=1:numel(fields), OUT_accum.(fields{f}) = 0; end
    end
    for f=1:numel(fields)
        OUT_accum.(fields{f}) = OUT_accum.(fields{f}) + OUTs.(fields{f});
    end
end

for f=1:numel(fields)
    OUT.(fields{f}) = OUT_accum.(fields{f})/S;
end
fprintf('All segments done in %.1fs.\n\n', toc(tStartAll));

%%  画图

%  图 5（彩色在线） 
figure('Name','图5 — 不同 k 下发送/接收神经元通道平均代谢能耗','Color','w');
plot(OUT.k, OUT.Pchan_pre_mean, '-o','LineWidth',1.5); hold on;
plot(OUT.k, OUT.Pchan_post_mean,'-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('平均代谢能量消耗（通道总和，pJ/ms·cm^{-2}）');
title('图5.（彩色在线）k 变化下发送/接收神经元离子通道的平均代谢能耗（单向耦合）');
legend({'发送神经元（s）','接收神经元（r）'},'Location','northwest','Box','off');
grid on;

% 插图：两个神经元的平均放电频率（Hz）
axes('Position',[0.50 0.30 0.32 0.32]);
plot(OUT.k, OUT.Fr_pre, '-','LineWidth',1.2); hold on;
plot(OUT.k, OUT.Fr_post,'--','LineWidth',1.2);
xlabel('k'); ylabel('Hz'); title('平均放电频率');
legend({'s','r'},'Location','best','Box','off'); grid on;

%  图 6（彩色在线） 
figure('Name','图6 — 从细胞外介质供应与突触后部耗散（功率）','Color','w');
plot(OUT.k, OUT.Psyn_post_mean, '-^','LineWidth',1.5); hold on;
plot(OUT.k, OUT.Pchan_post_mean,'-s','LineWidth',1.5);
xlabel('突触电导  k  (mS/cm^2)');
ylabel('功率（nJ/s·cm^{-2} ≡ pJ/ms·cm^{-2}）');
title('图6.（彩色在线）细胞外介质供应（突触→后）与突触后部耗散（单向耦合）');
legend({'从突触供应到后神经元','突触后膜通道耗散（Joule）'},'Location','northwest','Box','off');
grid on;

% 插图：连接处净能量导数（恒为正）
axes('Position',[0.50 0.30 0.32 0.32]);
plot(OUT.k, OUT.Pjunc_mean,'-','LineWidth',1.2);
xlabel('k'); ylabel('pJ/ms·cm^{-2}'); title('连接处净能量导数  \langle k(V_s - V_r)^2 \rangle');
grid on;


%  将结构体编码为 JSON 字符串
% 'PrettyPrint', true 会添加换行和缩进，使 .txt 文件更易读
jsonStr = jsonencode(OUT, 'PrettyPrint', true);

%  将字符串写入 .txt 文件
filename = 'ProDoubleOUT.txt';
fid = fopen(filename, 'w', 'n', 'UTF-8'); % 使用 UTF-8 编码以支持中文
fprintf(fid, '%s', jsonStr);
fclose(fid);

disp(['结构体已保存到: ', filename]);
