function OUT_k = run_HH_coupling_single_k(k, C,E,g,options)
% run_HH_coupling_single_k — 运行单个 k 值的 HH 耦合仿真
% 
% (已修复: 移除了 Alpha, Beta 参数)
%
% 输入:
%   k: (double) 单个电导值
%   C,E,g: HH 模型参数
%   options: (struct) 包含 fs, T, T_warm, mu_*, sigma_*, seed 的结构体
%
% 输出:
%   OUT_k: (struct) 包含该 k 值所有平均结果的结构体
%
    % 从 options 结构体中提取参数
    fs = options.fs;
    dt = 1000/fs;
    t  = 0:dt:options.T; % 确保 N 的计算与 options.T 匹配
    N  = numel(t);
    
    % 使用唯一的随机种子
    rng(options.seed + round(k*1e6)); 
    % 初始条件 (Alpha/Beta 函数仅在此处使用一次，然后被丢弃)
    V0 = 0;
    a_m = (2.5 - 0.1*V0) ./ (exp(2.5 - 0.1*V0) - 1);
    b_m = 4 .* exp(-V0./18);
    a_h = 0.07 .* exp(-V0./20);
    b_h = 1 ./ (exp(3 - 0.1*V0) + 1);
    a_n = (0.1 - 0.01*V0) ./ (exp(1 - 0.1*V0) - 1);
    b_n = 0.125 .* exp(-V0./80);
    
    m0 = a_m/(a_m+b_m);
    h0 = a_h/(a_h+b_h);
    n0 = a_n/(a_n+b_n);

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
    % RK4 核心循环
    for i = 1:N-1
        I_couple_i = k*(Vpre(i) - Vpost(i));
        % (已修复: 不再传递 Alpha, Beta)
        % 发送端
        [dV1_pre, dm1_pre, dh1_pre, dn1_pre] = HH_rhs(Vpre(i), mpre(i), hpre(i), npre(i), Ipre(i), 0, C,E,g);
        % 接收端（含耦合）
        [dV1_post,dm1_post,dh1_post,dn1_post] = HH_rhs(Vpost(i),mpost(i),hpost(i),npost(i),Ipost(i),I_couple_i, C,E,g);
        
        Vpre_mid   = Vpre(i)  + 0.5*dt*dV1_pre;
        Vpost_mid  = Vpost(i) + 0.5*dt*dV1_post;
        I_couple_m = k*(Vpre_mid - Vpost_mid);
        [dV2_pre, dm2_pre, dh2_pre, dn2_pre] = HH_rhs(Vpre_mid,  mpre(i)+0.5*dt*dm1_pre,  hpre(i)+0.5*dt*dh1_pre,  npre(i)+0.5*dt*dn1_pre,  Ipre(i),  0, C,E,g);
        [dV2_post,dm2_post,dh2_post,dn2_post] = HH_rhs(Vpost_mid, mpost(i)+0.5*dt*dm1_post, hpost(i)+0.5*dt*dh1_post, npost(i)+0.5*dt*dn1_post, Ipost(i), I_couple_m, C,E,g);
        
        Vpre_mid2  = Vpre(i)  + 0.5*dt*dV2_pre;
        Vpost_mid2 = Vpost(i) + 0.5*dt*dV2_post;
        I_couple_m2= k*(Vpre_mid2 - Vpost_mid2);
        [dV3_pre, dm3_pre, dh3_pre, dn3_pre] = HH_rhs(Vpre_mid2,  mpre(i)+0.5*dt*dm2_pre,  hpre(i)+0.5*dt*dh2_pre,  npre(i)+0.5*dt*dn2_pre,  Ipre(i),  0, C,E,g);
        [dV3_post,dm3_post,dh3_post,dn3_post] = HH_rhs(Vpost_mid2, mpost(i)+0.5*dt*dm2_post, hpost(i)+0.5*dt*dh2_post, npost(i)+0.5*dt*dn2_post, Ipost(i), I_couple_m2, C,E,g);
        
        Vpre_end   = Vpre(i)  + dt*dV3_pre;
        Vpost_end  = Vpost(i) + dt*dV3_post;
        I_couple_e = k*(Vpre_end - Vpost_end);
        [dV4_pre, dm4_pre, dh4_pre, dn4_pre] = HH_rhs(Vpre_end,  mpre(i)+dt*dm3_pre,  hpre(i)+dt*dh3_pre,  npre(i)+dt*dn3_pre,  Ipre(i),  0, C,E,g);
        [dV4_post,dm4_post,dh4_post,dn4_post] = HH_rhs(Vpost_end, mpost(i)+dt*dm3_post, hpost(i)+dt*dh3_post, npost(i)+dt*dn3_post, Ipost(i), I_couple_e, C,E,g);
        
        % 状态更新
        Vpre(i+1)  = Vpre(i)  + dt/6*(dV1_pre  + 2*dV2_pre  + 2*dV3_pre  + dV4_pre);
        mpre(i+1)  = mpre(i)  + dt/6*(dm1_pre  + 2*dm2_pre  + 2*dm3_pre  + dm4_pre);
        hpre(i+1)  = hpre(i)  + dt/6*(dh1_pre  + 2*dh2_pre  + 2*dh3_pre  + dh4_pre);
        npre(i+1)  = npre(i)  + dt/6*(dn1_pre  + 2*dn2_pre  + 2*dn3_pre  + dn4_pre);
        Vpost(i+1) = Vpost(i) + dt/6*(dV1_post + 2*dV2_post + 2*dV3_post + dV4_post);
        mpost(i+1) = mpost(i) + dt/6*(dm1_post + 2*dm2_post + 2*dm3_post + dm4_post);
        hpost(i+1) = hpost(i) + dt/6*(dh1_post + 2*dh2_post + 2*dh3_post + dh4_post);
        npost(i+1) = npost(i) + dt/6*(dn1_post + 2*dn2_post + 2*dn3_post + dn4_post);
        
        % --- 计算功率 (与您的版本相同) ---
        Pchan_pre(i)  = g.Na*(mpre(i)^3*hpre(i)) *(Vpre(i) - E.Na)^2 ...
                      + g.K *(npre(i)^4)        *(Vpre(i) - E.K )^2 ...
                      + g.L * (Vpre(i) - E.L )^2;
        Pchan_post(i) = g.Na*(mpost(i)^3*hpost(i))*(Vpost(i) - E.Na)^2 ...
                      + g.K *(npost(i)^4)        *(Vpost(i) - E.K )^2 ...
                      + g.L * (Vpost(i) - E.L )^2;
        
        I_cpl = k*(Vpre(i) - Vpost(i));
        P_amp(i)        = k * Vpre(i) * (Vpre(i) - Vpost(i));
        P_diss_post(i)  = k * Vpost(i) * (Vpre(i) - Vpost(i)); 
        P_net(i)        = P_amp(i) + P_diss_post(i);
    end
    
    % 丢弃热启动
    idx_warm = max(1, round(options.T_warm/dt)+1);
    seg = idx_warm:N;
    
    % 平均并打包到输出结构体
    OUT_k.k                 = k;
    OUT_k.Pchan_pre_mean    = mean(Pchan_pre(seg));
    OUT_k.Pchan_post_mean   = mean(Pchan_post(seg));
    OUT_k.Fr_pre            = count_rate(Vpre(seg),  fs, 20);
    OUT_k.Fr_post           = count_rate(Vpost(seg), fs, 20);
    OUT_k.P_amp_mean        = mean(P_amp(seg));
    OUT_k.P_diss_post_mean  = mean(P_diss_post(seg));
    OUT_k.P_net_mean        = mean(P_net(seg));
    if options.Ifdebug
        fprintf('k=%.3f 完成. (Fr_pre=%.2f Hz, P_pre=%.2f)\n', k, OUT_k.Fr_pre, OUT_k.Pchan_pre_mean);
    end
end
% ===== 局部函数 (必须包含在文件末尾) =====
function [dV, dm, dh, dn, iNa, iK, iL] = HH_rhs(V,m,h,n, I_ext, I_couple, C,E,g)
    % (已修复: 移除了 Alpha, Beta 参数)
    % (为 codegen 优化：将 Alpha/Beta 内联)
    % 速率函数
    a_m = (2.5 - 0.1*V) ./ (exp(2.5 - 0.1*V) - 1);
    b_m = 4 .* exp(-V./18);
    a_h = 0.07 .* exp(-V./20);
    b_h = 1 ./ (exp(3 - 0.1*V) + 1);
    a_n = (0.1 - 0.01*V) ./ (exp(1 - 0.1*V) - 1);
    b_n = 0.125 .* exp(-V./80);
    % 电流
    iNa = g.Na*(m^3*h)*(V - E.Na);
    iK  = g.K *(n^4  )*(V - E.K );
    iL  = g.L *(V - E.L);
    % 导数
    dV  = ( -iNa - iK - iL + I_ext + I_couple ) / C;
    dm  = a_m.*(1-m) - b_m.*m;
    dh  = a_h.*(1-h) - b_h.*h;
    dn  = a_n.*(1-n) - b_n.*n;
end
function Hz = count_rate(V, fs, thr)
    V = V(:)';
    x  = (V(1:end-1) < thr) & (V(2:end) >= thr);
    nspk = sum(x);
    Tsec = numel(V)/fs;
    Hz = nspk / Tsec;
end