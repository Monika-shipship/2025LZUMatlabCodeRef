
function OUT = Task5TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,options)
    % Task5TwoHH_EnergyCoupled — 两个电耦合 HH 神经元的能量平衡（单向耦合）
    % 说明：
    %   * 模仿 Task1/Task2 的风格，采用 RK4；单位与您现有脚本一致（V：mV，g：mS/cm^2，I：uA/cm^2）。
    %   * 单向电耦合（gap-like）：I_couple = k*(V_pre - V_post) 仅加到接收神经元；发送端不受反馈。
    %   * 通道代谢能耗采用 Hasegawa(2011) 建议的“Joule heat”正定形式：
    %       P_chan = g_Na m^3 h (V-E_Na)^2 + g_K n^4 (V-E_K)^2 + g_L (V-E_L)^2
    %     在 HH 的常用单位下，P_chan 的数值等于 “pJ/ms·cm^{-2}”（亦即 nW/cm^2）。
    %   * “从突触供应的能量”：P_syn_post = V_post * I_couple（功率为正表示向后神经元净输入）；
    %     “突触连接处净能量导数”（插图）：P_junc = k*(V_pre - V_post).^2（恒非负，表示连接电阻上的能量通量/耗散）。
    %
    % 输入（options 结构体关键字段）：
    %   fs         采样频率（Hz），如 100000 对应 0.01 ms 步长
    %   T          总时长（ms）
    %   k_list     耦合电导列表（mS/cm^2）
    %   mu_pre     发送端外部电流均值（uA/cm^2）
    %   sigma_pre  发送端外部电流标准差（uA/cm^2）
    %   mu_post    接收端外部电流均值（uA/cm^2）
    %   sigma_post 接收端外部电流标准差（uA/cm^2）
    %   seed       随机种子
    %   T_warm     统计前丢弃的“热启动”时长（ms），如 100
    %   Ifdebug    调试打印开关（0/1）
    %   Nseg 分多少段平均
    % 输出（OUT 结构体）：
    %   .k                  （1xK）各耦合电导
    %   .Pchan_pre_mean     （1xK）发送端通道平均功率（pJ/ms·cm^{-2}）
    %   .Pchan_post_mean    （1xK）接收端通道平均功率（pJ/ms·cm^{-2}）
    %   .Fr_pre             （1xK）发送端平均放电频率（Hz）
    %   .Fr_post            （1xK）接收端平均放电频率（Hz）
    %   .Psyn_post_mean     （1xK）从突触向接收端供应的平均功率（nJ/s·cm^{-2}）= 与 pJ/ms 等值
    %   .Pjunc_mean         （1xK）连接处净能量导数（k*(ΔV)^2）的平均值（pJ/ms·cm^{-2}）
    %
    % 参考：
    %   Moujahid et al., Phys. Rev. E 83, 031912 (2011)；Hasegawa, arXiv:1106.5862（建议采用包含 Joule 热的功率定义）。

    arguments
        C
        E
        g
        Alpha
        Beta
        options.fs (1,1) double = 1e5
        options.T (1,1) double = 300
        options.k_list double = linspace(0,2.0,9)
        options.mu_pre (1,1) double = 7.5
        options.sigma_pre (1,1) double = 1.0
        options.mu_post (1,1) double = 0.0
        options.sigma_post (1,1) double = 1.0
        options.seed (1,1) double = 2025
        options.T_warm (1,1) double = 100
        options.Ifdebug (1,1) double = 0
        options.Nseg (1,1) double = 1
    end

    fs = options.fs;
    dt = 1000/fs;              % 步长（ms）
    t  = 0:dt:options.T;       % 时间轴（ms）
    N  = numel(t);

    % 预设初始条件（按静息闭式表达）
    V0 = 0; % HH 约定：静息 0 mV
    m0 = Alpha.m(V0)/(Alpha.m(V0)+Beta.m(V0));
    h0 = Alpha.h(V0)/(Alpha.h(V0)+Beta.h(V0));
    n0 = Alpha.n(V0)/(Alpha.n(V0)+Beta.n(V0));

    % 结果容器
    K = numel(options.k_list);
    Pchan_pre_mean  = zeros(1,K);
    Pchan_post_mean = zeros(1,K);
    Fr_pre  = zeros(1,K);
    Fr_post = zeros(1,K);
    Psyn_post_mean = zeros(1,K);
    Pjunc_mean     = zeros(1,K);

    % 随机数
    rng(options.seed);

    % ===== 主循环：不同 k =====
    for ik = 1:K
        k = options.k_list(ik);

        % 状态初始化
        Vpre = zeros(1,N); Vpre(1) = V0;
        mpre = zeros(1,N); mpre(1) = m0;
        hpre = zeros(1,N); hpre(1) = h0;
        npre = zeros(1,N); npre(1) = n0;

        Vpost = zeros(1,N); Vpost(1) = V0;
        mpost = zeros(1,N); mpost(1) = m0;
        hpost = zeros(1,N); hpost(1) = h0;
        npost = zeros(1,N); npost(1) = n0;

        % 外部电流（高斯白噪，步进采样）
        Ipre  = options.mu_pre  + options.sigma_pre  * randn(1,N);
        Ipost = options.mu_post + options.sigma_post * randn(1,N);

        % 记录功率
        Pchan_pre  = zeros(1,N);
        Pchan_post = zeros(1,N);
        Psyn_post  = zeros(1,N);
        Pjunc      = zeros(1,N);

        % ===== RK4 演化 =====
        for i = 1:N-1
            % 耦合电流（只注入到 post）
            I_couple_i = k*(Vpre(i) - Vpost(i));

            % --- 发送端 ---
            [dV1_pre, dm1_pre, dh1_pre, dn1_pre, iNa_pre, iK_pre, iL_pre] = HH_rhs(Vpre(i),  mpre(i),  hpre(i),  npre(i),  Ipre(i),  0, C,E,g,Alpha,Beta);
            % --- 接收端（含耦合） ---
            [dV1_post,dm1_post,dh1_post,dn1_post,iNa_post,iK_post,iL_post] = HH_rhs(Vpost(i),mpost(i),hpost(i),npost(i),Ipost(i),I_couple_i, C,E,g,Alpha,Beta);

            % 中点电流（噪声按样本常数近似；耦合用中点电位估计）
            Vpre_mid  = Vpre(i)  + 0.5*dt*dV1_pre;
            Vpost_mid = Vpost(i) + 0.5*dt*dV1_post;
            I_couple_mid = k*(Vpre_mid - Vpost_mid);

            [dV2_pre, dm2_pre, dh2_pre, dn2_pre, ~,~,~] = HH_rhs(Vpre_mid,  mpre(i)+0.5*dt*dm1_pre,  hpre(i)+0.5*dt*dh1_pre,  npre(i)+0.5*dt*dn1_pre,  Ipre(i),  0, C,E,g,Alpha,Beta);
            [dV2_post,dm2_post,dh2_post,dn2_post,~,~,~] = HH_rhs(Vpost_mid, mpost(i)+0.5*dt*dm1_post, hpost(i)+0.5*dt*dh1_post, npost(i)+0.5*dt*dn1_post, Ipost(i), I_couple_mid, C,E,g,Alpha,Beta);

            Vpre_mid2  = Vpre(i)  + 0.5*dt*dV2_pre;
            Vpost_mid2 = Vpost(i) + 0.5*dt*dV2_post;
            I_couple_mid2 = k*(Vpre_mid2 - Vpost_mid2);

            [dV3_pre, dm3_pre, dh3_pre, dn3_pre, ~,~,~] = HH_rhs(Vpre_mid2,  mpre(i)+0.5*dt*dm2_pre,  hpre(i)+0.5*dt*dh2_pre,  npre(i)+0.5*dt*dn2_pre,  Ipre(i),  0, C,E,g,Alpha,Beta);
            [dV3_post,dm3_post,dh3_post,dn3_post,~,~,~] = HH_rhs(Vpost_mid2, mpost(i)+0.5*dt*dm2_post, hpost(i)+0.5*dt*dh2_post, npost(i)+0.5*dt*dn2_post, Ipost(i), I_couple_mid2, C,E,g,Alpha,Beta);

            Vpre_end  = Vpre(i)  + dt*dV3_pre;
            Vpost_end = Vpost(i) + dt*dV3_post;
            I_couple_end = k*(Vpre_end - Vpost_end);

            [dV4_pre, dm4_pre, dh4_pre, dn4_pre, iNa4_pre, iK4_pre, iL4_pre] = HH_rhs(Vpre_end,  mpre(i)+dt*dm3_pre,  hpre(i)+dt*dh3_pre,  npre(i)+dt*dn3_pre,  Ipre(i),  0, C,E,g,Alpha,Beta);
            [dV4_post,dm4_post,dh4_post,dn4_post,iNa4_post,iK4_post,iL4_post] = HH_rhs(Vpost_end, mpost(i)+dt*dm3_post, hpost(i)+dt*dh3_post, npost(i)+dt*dn3_post, Ipost(i), I_couple_end, C,E,g,Alpha,Beta);

            % 状态更新
            Vpre(i+1)  = Vpre(i)  + dt/6*(dV1_pre + 2*dV2_pre + 2*dV3_pre + dV4_pre);
            mpre(i+1)  = mpre(i)  + dt/6*(dm1_pre + 2*dm2_pre + 2*dm3_pre + dm4_pre);
            hpre(i+1)  = hpre(i)  + dt/6*(dh1_pre + 2*dh2_pre + 2*dh3_pre + dh4_pre);
            npre(i+1)  = npre(i)  + dt/6*(dn1_pre + 2*dn2_pre + 2*dn3_pre + dn4_pre);

            Vpost(i+1) = Vpost(i) + dt/6*(dV1_post + 2*dV2_post + 2*dV3_post + dV4_post);
            mpost(i+1) = mpost(i) + dt/6*(dm1_post + 2*dm2_post + 2*dm3_post + dm4_post);
            hpost(i+1) = hpost(i) + dt/6*(dh1_post + 2*dh2_post + 2*dh3_post + dh4_post);
            npost(i+1) = npost(i) + dt/6*(dn1_post + 2*dn2_post + 2*dn3_post + dn4_post);

            % 实时功率（采用 Joule 形式；单位：pJ/ms·cm^{-2} = nW/cm^2）
            % 用末段估计（更稳定）：
            Pchan_pre(i)  = g.Na*(mpre(i)^3*hpre(i))*(Vpre(i)-E.Na)^2 + g.K*(npre(i)^4)*(Vpre(i)-E.K)^2 + g.L*(Vpre(i)-E.L)^2;
            Pchan_post(i) = g.Na*(mpost(i)^3*hpost(i))*(Vpost(i)-E.Na)^2 + g.K*(npost(i)^4)*(Vpost(i)-E.K)^2 + g.L*(Vpost(i)-E.L)^2;

            % 突触向后神经元的供能功率（可能正或负）：
            I_couple_i = k*(Vpre(i) - Vpost(i));
            % 供能：外介质连接 到 后神经元（取"流入后神经元"为正）
            Psyn_post(i) = Vpost(i) * I_couple_i;%% 以向 post 净输入为正

            % 连接处净能量导数（总为非负）：
            Pjunc(i) = k*(Vpre(i) - Vpost(i))^2;
            % Pjunc(i) = Psyn_post(i) - k*(Vpre(i) - Vpost(i))^2;  % Pjunc(i)= -V_r I_cpl  + ( - g_cpl (ΔV)^2 )
        end

        % 丢弃热启动部分
        idx_warm = round(options.T_warm/dt)+1;
        seg = idx_warm:N;

        Pchan_pre_mean(ik)  = mean(Pchan_pre(seg));
        Pchan_post_mean(ik) = mean(Pchan_post(seg));
        Psyn_post_mean(ik)  = mean(Psyn_post(seg));
        Pjunc_mean(ik)      = mean(Pjunc(seg));

        % 放电频率（阈值 20 mV 上穿计数）
        Fr_pre(ik)  = count_rate(Vpre(seg), fs, 20);
        Fr_post(ik) = count_rate(Vpost(seg), fs, 20);

        if options.Ifdebug
            fprintf('[k=%.3f] 预: P=%.3f  Hz=%.2f | 后: P=%.3f  Psyn=%.3f  Hz=%.2f  (Pjunc=%.3f)\n', ...
                k, Pchan_pre_mean(ik), Fr_pre(ik), Pchan_post_mean(ik), Psyn_post_mean(ik), Fr_post(ik), Pjunc_mean(ik));
        end
    end

    % 输出
    OUT.k               = options.k_list;
    OUT.Pchan_pre_mean  = Pchan_pre_mean;
    OUT.Pchan_post_mean = Pchan_post_mean;
    OUT.Fr_pre          = Fr_pre;
    OUT.Fr_post         = Fr_post;
    OUT.Psyn_post_mean  = Psyn_post_mean;
    OUT.Pjunc_mean      = Pjunc_mean;
end

% ====== 局部函数 ======

function [dV, dm, dh, dn, iNa, iK, iL] = HH_rhs(V,m,h,n, I_ext, I_couple, C,E,g,Alpha,Beta)
    % 通道电流
    iNa = g.Na*(m^3*h)*(V - E.Na);
    iK  = g.K*(n^4)*(V - E.K);
    iL  = g.L*(V - E.L);

    % 膜电位
    dV = ( -iNa - iK - iL + I_ext + I_couple ) / C;

    % Gating
    dm = Alpha.m(V)*(1-m) - Beta.m(V)*m;
    dh = Alpha.h(V)*(1-h) - Beta.h(V)*h;
    dn = Alpha.n(V)*(1-n) - Beta.n(V)*n;
end

function Hz = count_rate(V, fs, thr)
    V = V(:)';
    x  = (V(1:end-1) < thr) & (V(2:end) >= thr);
    nspk = sum(x);
    Tsec = numel(V)/fs;
    Hz = nspk / Tsec;
end
