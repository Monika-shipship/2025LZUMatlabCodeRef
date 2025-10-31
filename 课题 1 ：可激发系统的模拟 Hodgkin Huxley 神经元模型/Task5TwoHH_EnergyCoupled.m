function OUT = Task5TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,opts)
    % Task5TwoHH_EnergyCoupled — 两个电耦合 HH 神经元（单向耦合）单段演化与统计
    % 说明：
    %   * 单向耦合：I_cpl = k*(Vpre - Vpost) 只注入接收端膜方程；发送端不受反馈
    %   * 代谢能耗采用"Joule 热"形式：
    %       P_chan = g_Na m^3 h (V-E_Na)^2 + g_K n^4 (V-E_K)^2 + g_L (V-E_L)^2
    %   * 突触能量项按论文（等式 (8)）的物理解读：
    %       供能（外介质/放大器）P_amp   = V_pre * I_cpl = k * Vpre * (Vpre - Vpost)
    %       连接处耗散（Joule）   P_diss  = - k * (Vpre - Vpost).^2
    %     插图"净能量导数"取正定形式：k*(Vpre - Vpost).^2，以便与文中"为正"表述对应
    %
    % 输入（名称-值）：
    %   fs, T, k_list(标量), mu_pre, sigma_pre, mu_post, sigma_post, seed, T_warm, Ifdebug
    %
    arguments
        C
        E
        g
        Alpha
        Beta
        opts.fs (1,1) double = 1e5
        opts.T (1,1) double = 300
        opts.k_list double = 0.1
        opts.mu_pre (1,1) double = 0
        opts.sigma_pre (1,1) double = 3
        opts.mu_post (1,1) double = 0
        opts.sigma_post (1,1) double = 1
        opts.seed (1,1) double = 2025
        opts.T_warm (1,1) double = 100
        opts.Ifdebug (1,1) double = 0
    end

    fs = opts.fs;
    dt = 1000/fs;                  % ms
    t  = 0:dt:opts.T;              % ms
    N  = numel(t);
    k  = opts.k_list;              % 此函数单次只接收一个 k

    % 初始条件（静息值）
    V0 = 0;
    m0 = Alpha.m(V0)/(Alpha.m(V0)+Beta.m(V0));
    h0 = Alpha.h(V0)/(Alpha.h(V0)+Beta.h(V0));
    n0 = Alpha.n(V0)/(Alpha.n(V0)+Beta.n(V0));

    % 预分配
    Vpre = zeros(1,N); Vpre(1) = V0;
    mpre = zeros(1,N); mpre(1) = m0;
    hpre = zeros(1,N); hpre(1) = h0;
    npre = zeros(1,N); npre(1) = n0;

    Vpost = zeros(1,N); Vpost(1) = V0;
    mpost = zeros(1,N); mpost(1) = m0;
    hpost = zeros(1,N); hpost(1) = h0;
    npost = zeros(1,N); npost(1) = n0;

    % 可复现的随机数
    rng(opts.seed);
    Ipre  = opts.mu_pre  + opts.sigma_pre  * randn(1,N);
    Ipost = opts.mu_post + opts.sigma_post * randn(1,N);

    % 功率记录
    Pchan_pre  = zeros(1,N);
    Pchan_post = zeros(1,N);
    Pamp       = zeros(1,N);   % 放大器供能
    Pjunc_pos  = zeros(1,N);   % k*(ΔV)^2 的正定量，用于插图
    % 连接处耗散项如需单独观察，可按需要外送；此处统计的是放大器供能与通道耗散

    % 演化（经典 RK4）
    for i = 1:N-1
        % step 1
        Ic_i = k*(Vpre(i) - Vpost(i));
        [dV1_pre, dm1_pre, dh1_pre, dn1_pre] = HH_rhs(Vpre(i),  mpre(i),  hpre(i),  npre(i),  Ipre(i),  0,       C,E,g,Alpha,Beta);
        [dV1_pos, dm1_pos, dh1_pos, dn1_pos] = HH_rhs(Vpost(i), mpost(i), hpost(i), npost(i), Ipost(i), Ic_i,    C,E,g,Alpha,Beta);

        % step 2（中点）
        Vpre_m  = Vpre(i)  + 0.5*dt*dV1_pre;
        Vpost_m = Vpost(i) + 0.5*dt*dV1_pos;
        Ic_m    = k*(Vpre_m - Vpost_m);
        [dV2_pre, dm2_pre, dh2_pre, dn2_pre] = HH_rhs(Vpre_m,  mpre(i)+0.5*dt*dm1_pre,  hpre(i)+0.5*dt*dh1_pre,  npre(i)+0.5*dt*dn1_pre,  Ipre(i),  0,    C,E,g,Alpha,Beta);
        [dV2_pos, dm2_pos, dh2_pos, dn2_pos] = HH_rhs(Vpost_m, mpost(i)+0.5*dt*dm1_pos, hpost(i)+0.5*dt*dh1_pos, npost(i)+0.5*dt*dn1_pos, Ipost(i), Ic_m, C,E,g,Alpha,Beta);

        % step 3（中点）
        Vpre_m2  = Vpre(i)  + 0.5*dt*dV2_pre;
        Vpost_m2 = Vpost(i) + 0.5*dt*dV2_pos;
        Ic_m2    = k*(Vpre_m2 - Vpost_m2);
        [dV3_pre, dm3_pre, dh3_pre, dn3_pre] = HH_rhs(Vpre_m2,  mpre(i)+0.5*dt*dm2_pre,  hpre(i)+0.5*dt*dh2_pre,  npre(i)+0.5*dt*dn2_pre,  Ipre(i),  0,     C,E,g,Alpha,Beta);
        [dV3_pos, dm3_pos, dh3_pos, dn3_pos] = HH_rhs(Vpost_m2, mpost(i)+0.5*dt*dm2_pos, hpost(i)+0.5*dt*dh2_pos, npost(i)+0.5*dt*dn2_pos, Ipost(i), Ic_m2, C,E,g,Alpha,Beta);

        % step 4（端点）
        Vpre_e  = Vpre(i)  + dt*dV3_pre;
        Vpost_e = Vpost(i) + dt*dV3_pos;
        Ic_e    = k*(Vpre_e - Vpost_e);
        [dV4_pre, dm4_pre, dh4_pre, dn4_pre] = HH_rhs(Vpre_e,  mpre(i)+dt*dm3_pre,  hpre(i)+dt*dh3_pre,  npre(i)+dt*dn3_pre,  Ipre(i),  0,     C,E,g,Alpha,Beta);
        [dV4_pos, dm4_pos, dh4_pos, dn4_pos] = HH_rhs(Vpost_e, mpost(i)+dt*dm3_pos, hpost(i)+dt*dh3_pos, npost(i)+dt*dn3_pos, Ipost(i), Ic_e, C,E,g,Alpha,Beta);

        % 状态更新
        Vpre(i+1)  = Vpre(i)  + dt/6*(dV1_pre + 2*dV2_pre + 2*dV3_pre + dV4_pre);
        mpre(i+1)  = mpre(i)  + dt/6*(dm1_pre + 2*dm2_pre + 2*dm3_pre + dm4_pre);
        hpre(i+1)  = hpre(i)  + dt/6*(dh1_pre + 2*dh2_pre + 2*dh3_pre + dh4_pre);
        npre(i+1)  = npre(i)  + dt/6*(dn1_pre + 2*dn2_pre + 2*dn3_pre + dn4_pre);

        Vpost(i+1) = Vpost(i) + dt/6*(dV1_pos + 2*dV2_pos + 2*dV3_pos + dV4_pos);
        mpost(i+1) = mpost(i) + dt/6*(dm1_pos + 2*dm2_pos + 2*dm3_pos + dm4_pos);
        hpost(i+1) = hpost(i) + dt/6*(dh1_pos + 2*dh2_pos + 2*dh3_pos + dh4_pos);
        npost(i+1) = npost(i) + dt/6*(dn1_pos + 2*dn2_pos + 2*dn3_pos + dn4_pos);

        % 通道功率（取当前时刻，正定）
        Pchan_pre(i)  = g.Na*(mpre(i)^3*hpre(i))*(Vpre(i)-E.Na)^2 + g.K*(npre(i)^4)*(Vpre(i)-E.K)^2 + g.L*(Vpre(i)-E.L)^2;
        Pchan_post(i) = g.Na*(mpost(i)^3*hpost(i))*(Vpost(i)-E.Na)^2 + g.K*(npost(i)^4)*(Vpost(i)-E.K)^2 + g.L*(Vpost(i)-E.L)^2;

        % 突触能量项（见上）
        Ic  = k*(Vpre(i) - Vpost(i));
        Pamp(i)      =  Vpre(i) * Ic;             % 放大器供能，正负皆可；图 6 取平均
        Pjunc_pos(i) =  k*(Vpre(i) - Vpost(i))^2; % 净能量导数（正定）用于插图
        % 如需连接处耗散可另存：P_diss = -k*(ΔV)^2
    end

    % 丢弃热启动部分
    idx_warm = round(opts.T_warm/dt) + 1;
    seg = idx_warm:N;

    % 汇总
    OUT.k               = k;
    OUT.Pchan_pre_mean  = mean(Pchan_pre(seg));
    OUT.Pchan_post_mean = mean(Pchan_post(seg));
    OUT.Psyn_post_mean  = mean(Pamp(seg));       % "从外介质供能"（放大器项）
    OUT.Pjunc_mean      = mean(Pjunc_pos(seg));  % 连接处净能量导数（正定）

    % 放电频率（阈值 20 mV 上穿计数）
    OUT.Fr_pre  = count_rate(Vpre(seg),  fs, 20);
    OUT.Fr_post = count_rate(Vpost(seg), fs, 20);

    if opts.Ifdebug
        fprintf('[k=%.3f] 预: P=%.3f Hz=%.2f | 后: P=%.3f Psyn=%.3f Hz=%.2f (Pjunc+=%.3f)\n', ...
            k, OUT.Pchan_pre_mean, OUT.Fr_pre, OUT.Pchan_post_mean, OUT.Psyn_post_mean, OUT.Fr_post, OUT.Pjunc_mean);
    end
end

% 局部函数
function [dV, dm, dh, dn] = HH_rhs(V,m,h,n, I_ext, I_cpl, C,E,g,Alpha,Beta)
    iNa = g.Na*(m^3*h)*(V - E.Na);
    iK  = g.K*(n^4)*(V - E.K);
    iL  = g.L*(V - E.L);
    dV  = ( -iNa - iK - iL + I_ext + I_cpl ) / C;
    dm  = Alpha.m(V)*(1-m) - Beta.m(V)*m;
    dh  = Alpha.h(V)*(1-h) - Beta.h(V)*h;
    dn  = Alpha.n(V)*(1-n) - Beta.n(V)*n;
end

function Hz = count_rate(V, fs, thr)
    V = V(:)';
    x  = (V(1:end-1) < thr) & (V(2:end) >= thr);
    nspk = sum(x);
    Tsec = numel(V)/fs;
    Hz = nspk / Tsec;
end
