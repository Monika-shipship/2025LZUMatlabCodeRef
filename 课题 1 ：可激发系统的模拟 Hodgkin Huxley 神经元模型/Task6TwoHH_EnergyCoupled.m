function OUT = Task6TwoHH_EnergyCoupled(C,E,g,Alpha,Beta,options)
% Task5TwoHH_EnergyCoupled — 两个电耦合 HH 神经元的能量平衡（单向耦合）
% 说明：
%   * 单向电耦合：I_couple = k*(V_pre - V_post)，只注入到接收神经元；发送端不受反馈。
%   * 通道代谢能耗采用 "Joule 热"正定形式：
%       P_chan = g_Na m^3 h (V-E_Na)^2 + g_K n^4 (V-E_K)^2 + g_L (V-E_L)^2
%     在 V 用 mV、g 用 mS/cm^2 下，单位为 pJ/ms·cm^{-2}（即 nW/cm^2）。
%   * 突触能量平衡的两部分（与论文图 6 对应）：
%       P_amp       =  k * V_pre * (V_pre - V_post)           % 放大器向突触处供能（正）
%       P_diss_post = -k * (V_pre - V_post).^2                % 突触后部位耗散（负）
%       P_net       =  P_amp + P_diss_post = V_post * I_couple % 连接处净能量导数（正）

    arguments
        C
        E
        g
        Alpha
        Beta
        options.fs (1,1) double = 2e4
        options.T  (1,1) double = 3.0e4
        options.k_list double   = 0:0.01:0.20
        options.mu_pre (1,1) double = 0.0
        options.sigma_pre (1,1) double = 3.0
        options.mu_post (1,1) double = 0.0
        options.sigma_post (1,1) double = 1.0
        options.seed (1,1) double = 2025
        options.T_warm (1,1) double = 3.0e3
        options.Ifdebug (1,1) double = 1
    end

    fs = options.fs;
    dt = 1000/fs;           % 步长（ms）
    t  = 0:dt:options.T;    % 时间轴（ms）
    N  = numel(t);

    % 初始条件（静息 0 mV 的闭式解）
    V0 = 0;
    m0 = Alpha.m(V0)/(Alpha.m(V0)+Beta.m(V0));
    h0 = Alpha.h(V0)/(Alpha.h(V0)+Beta.h(V0));
    n0 = Alpha.n(V0)/(Alpha.n(V0)+Beta.n(V0));

    % 结果容器
    K = numel(options.k_list);
    Pchan_pre_mean   = zeros(1,K);
    Pchan_post_mean  = zeros(1,K);
    Fr_pre           = zeros(1,K);
    Fr_post          = zeros(1,K);
    P_amp_mean       = zeros(1,K);
    P_diss_post_mean = zeros(1,K);
    P_net_mean       = zeros(1,K);

    rng(options.seed);

    if options.Ifdebug
        fprintf('k 扫描：%d 点；总时长 %.3f s；dt=%.3f ms\n', K, options.T/1000, dt);
    end

    % 主循环：对每一个 k 做单次长时间仿真
    for ik = 1:K
        k = options.k_list(ik);

        % 状态
        Vpre  = zeros(1,N); Vpre(1)  = V0;
        mpre  = zeros(1,N); mpre(1)  = m0;
        hpre  = zeros(1,N); hpre(1)  = h0;
        npre  = zeros(1,N); npre(1)  = n0;

        Vpost = zeros(1,N); Vpost(1) = V0;
        mpost = zeros(1,N); mpost(1) = m0;
        hpost = zeros(1,N); hpost(1) = h0;
        npost = zeros(1,N); npost(1) = n0;

        % 高斯输入
        Ipre  = options.mu_pre  + options.sigma_pre  * randn(1,N);
        Ipost = options.mu_post + options.sigma_post * randn(1,N);

        % 即时功率
        Pchan_pre   = zeros(1,N);
        Pchan_post  = zeros(1,N);
        P_amp       = zeros(1,N);
        P_diss_post = zeros(1,N);
        P_net       = zeros(1,N);

        % RK4
        for i = 1:N-1
            I_couple_i = k*(Vpre(i) - Vpost(i));

            % 发送端
            [dV1_pre, dm1_pre, dh1_pre, dn1_pre, ~,~,~] = HH_rhs(Vpre(i), mpre(i), hpre(i), npre(i), Ipre(i), 0, C,E,g,Alpha,Beta);
            % 接收端（含耦合）
            [dV1_post,dm1_post,dh1_post,dn1_post,~,~,~] = HH_rhs(Vpost(i),mpost(i),hpost(i),npost(i),Ipost(i),I_couple_i, C,E,g,Alpha,Beta);

            Vpre_mid   = Vpre(i)  + 0.5*dt*dV1_pre;
            Vpost_mid  = Vpost(i) + 0.5*dt*dV1_post;
            I_couple_m = k*(Vpre_mid - Vpost_mid);

            [dV2_pre, dm2_pre, dh2_pre, dn2_pre, ~,~,~] = HH_rhs(Vpre_mid,  mpre(i)+0.5*dt*dm1_pre,  hpre(i)+0.5*dt*dh1_pre,  npre(i)+0.5*dt*dn1_pre,  Ipre(i),  0, C,E,g,Alpha,Beta);
            [dV2_post,dm2_post,dh2_post,dn2_post,~,~,~] = HH_rhs(Vpost_mid, mpost(i)+0.5*dt*dm1_post, hpost(i)+0.5*dt*dh1_post, npost(i)+0.5*dt*dn1_post, Ipost(i), I_couple_m, C,E,g,Alpha,Beta);

            Vpre_mid2  = Vpre(i)  + 0.5*dt*dV2_pre;
            Vpost_mid2 = Vpost(i) + 0.5*dt*dV2_post;
            I_couple_m2= k*(Vpre_mid2 - Vpost_mid2);

            [dV3_pre, dm3_pre, dh3_pre, dn3_pre, ~,~,~] = HH_rhs(Vpre_mid2,  mpre(i)+0.5*dt*dm2_pre,  hpre(i)+0.5*dt*dh2_pre,  npre(i)+0.5*dt*dn2_pre,  Ipre(i),  0, C,E,g,Alpha,Beta);
            [dV3_post,dm3_post,dh3_post,dn3_post,~,~,~] = HH_rhs(Vpost_mid2, mpost(i)+0.5*dt*dm2_post, hpost(i)+0.5*dt*dh2_post, npost(i)+0.5*dt*dn2_post, Ipost(i), I_couple_m2, C,E,g,Alpha,Beta);

            Vpre_end   = Vpre(i)  + dt*dV3_pre;
            Vpost_end  = Vpost(i) + dt*dV3_post;
            I_couple_e = k*(Vpre_end - Vpost_end);

            [dV4_pre, dm4_pre, dh4_pre, dn4_pre, ~,~,~] = HH_rhs(Vpre_end,  mpre(i)+dt*dm3_pre,  hpre(i)+dt*dh3_pre,  npre(i)+dt*dn3_pre,  Ipre(i),  0, C,E,g,Alpha,Beta);
            [dV4_post,dm4_post,dh4_post,dn4_post,~,~,~] = HH_rhs(Vpost_end, mpost(i)+dt*dm3_post, hpost(i)+dt*dh3_post, npost(i)+dt*dn3_post, Ipost(i), I_couple_e, C,E,g,Alpha,Beta);

            % 状态更新
            Vpre(i+1)  = Vpre(i)  + dt/6*(dV1_pre  + 2*dV2_pre  + 2*dV3_pre  + dV4_pre);
            mpre(i+1)  = mpre(i)  + dt/6*(dm1_pre  + 2*dm2_pre  + 2*dm3_pre  + dm4_pre);
            hpre(i+1)  = hpre(i)  + dt/6*(dh1_pre  + 2*dh2_pre  + 2*dh3_pre  + dh4_pre);
            npre(i+1)  = npre(i)  + dt/6*(dn1_pre  + 2*dn2_pre  + 2*dn3_pre  + dn4_pre);

            Vpost(i+1) = Vpost(i) + dt/6*(dV1_post + 2*dV2_post + 2*dV3_post + dV4_post);
            mpost(i+1) = mpost(i) + dt/6*(dm1_post + 2*dm2_post + 2*dm3_post + dm4_post);
            hpost(i+1) = hpost(i) + dt/6*(dh1_post + 2*dh2_post + 2*dh3_post + dh4_post);
            npost(i+1) = npost(i) + dt/6*(dn1_post + 2*dn2_post + 2*dn3_post + dn4_post);

            % 通道代谢功率（修复了 g.L 的非法表达式）
            Pchan_pre(i)  = g.Na*(mpre(i)^3*hpre(i)) *(Vpre(i) - E.Na)^2 ...
                          + g.K *(npre(i)^4)        *(Vpre(i) - E.K )^2 ...
                          + g.L *                   (Vpre(i) - E.L )^2;

            Pchan_post(i) = g.Na*(mpost(i)^3*hpost(i))*(Vpost(i) - E.Na)^2 ...
                          + g.K *(npost(i)^4)        *(Vpost(i) - E.K )^2 ...
                          + g.L *                    (Vpost(i) - E.L )^2;

            % 突触处能量平衡分解
            % I_cpl           = k*(Vpre(i) - Vpost(i));
            % P_amp(i)        =  k * Vpre(i) * (Vpre(i) - Vpost(i));   % 放大器供能（正）
            % P_diss_post(i)  = -k * (Vpre(i) - Vpost(i))^2;           % 后部位耗散（负）
            % P_net(i)        =  Vpost(i) * I_cpl;                     % 连接处净导数（正）
            % 突触处能量平衡分解 (修正依据 Eq. 8)
            
            I_cpl = k*(Vpre(i) - Vpost(i));

            % 供应项 (Eq. 8 第二项, 放大器供能, 论文图6蓝线)
            P_amp(i)        = k * Vpre(i) * (Vpre(i) - Vpost(i));
            
            % 耗散项 (Eq. 8 第一项, 突触后耗散, 论文图6红线)
            P_diss_post(i)  = k * Vpost(i) * (Vpre(i) - Vpost(i)); 

            % 净能量导数 (Eq. 8 两项之和, 论文图6插图)
            P_net(i)        = P_amp(i) + P_diss_post(i);
        end

        % 丢弃热启动
        idx_warm = max(1, round(options.T_warm/dt)+1);
        seg = idx_warm:N;

        % 平均
        Pchan_pre_mean(ik)   = mean(Pchan_pre(seg));
        Pchan_post_mean(ik)  = mean(Pchan_post(seg));
        P_amp_mean(ik)       = mean(P_amp(seg));
        P_diss_post_mean(ik) = mean(P_diss_post(seg));
        P_net_mean(ik)       = mean(P_net(seg));

        % 频率
        Fr_pre(ik)  = count_rate(Vpre(seg),  fs, 20);
        Fr_post(ik) = count_rate(Vpost(seg), fs, 20);

        % 进度
        if options.Ifdebug
            pct = 100*ik/K;
            fprintf('\r进度：%d/%d (%.1f%%)', ik, K, pct);
            if ik==K, fprintf('\n'); end
        end
    end

    % 输出
    OUT.k                 = options.k_list;
    OUT.Pchan_pre_mean    = Pchan_pre_mean;
    OUT.Pchan_post_mean   = Pchan_post_mean;
    OUT.Fr_pre            = Fr_pre;
    OUT.Fr_post           = Fr_post;
    OUT.P_amp_mean        = P_amp_mean;
    OUT.P_diss_post_mean  = P_diss_post_mean;
    OUT.P_net_mean        = P_net_mean;
end

% ===== 局部函数 =====

function [dV, dm, dh, dn, iNa, iK, iL] = HH_rhs(V,m,h,n, I_ext, I_couple, C,E,g,Alpha,Beta)
    iNa = g.Na*(m^3*h)*(V - E.Na);
    iK  = g.K *(n^4  )*(V - E.K );
    iL  = g.L *(V - E.L);
    dV  = ( -iNa - iK - iL + I_ext + I_couple ) / C;
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
