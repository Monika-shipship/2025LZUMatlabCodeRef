%%
%[text] 
%[text] ### 脉冲刺激-使用PPT中的参数和方程
clc;clear;

%[text] %[text:anchor:TMP_6c4e] PPT中给出的参数是
%[text] 电容：$C = 1 \\mu\\text{F/cm}^2$
C=1;
%[text] 各离子的能斯特电位和对应的电导：
%[text] $\\begin{array}{lrr}\n\\hline\\hline\nx & E\_x \\text{ \[mV\]} & g\_x \\text{ \[mS / cm}^2\\text{\]} \\\\\n\\hline\n\\text{Na} & 55 & 40 \\\\\n\\text{K} & -77 & 35 \\\\\n\\text{L} & -65 & 0.3 \\\\\n\\hline\\hline\n\\end{array}$
E.Na=55;E.K=-77;E.L=-65;
g.Na=40;g.K=35;g.L=0.3;
%[text] 对应的 alpha和beta函数：
%[text] $\n\n\n\n\n\n\\begin{array}{lcc}\n\\hline\\hline\nx & \\alpha\_x (u / \\text{ mV}) \\text{ \[ms}^{-1}\\text{\]} & \\beta\_x (u / \\cdot \\text{ mV}) \\text{ \[ms}^{-1}\\text{\]} \\\\\n\\hline\nn & 0.02(u-25) / \[1-e^{-(u-25)/9}\] & -0.002(u-25) / \[1-e^{(u-25)/9}\] \\\\\nm & 0.182(u+35) / \[1-e^{-(u+35)/9}\] & -0.124(u+35) / \[1-e^{(u+35)/9}\] \\\\\nh & 0.25 e^{-(u+90)/12} & 0.25 e^{(u+62)/6} / e^{(u+90)/12} \\\\\n\\hline\\hline\n\\end{array}$
Alpha.n=@(u) 0.02*(u-25)./(1-exp(-(u-25)/9));
Alpha.m=@(u) 0.182*(u+35)./(1-exp(-(u+35)/9));
Alpha.h=@(u) 0.25*exp(-(u+90)/12);
Beta.n=@(u) -0.002*(u-25)./(1-exp((u-25)/9));
Beta.m=@(u) -0.124*(u+35)./(1-exp((u+35)/9));
Beta.h=@(u) 0.25*exp((u+62)/6)./exp((u+90)/12);

%[text] 
fs=5000;%采样频率 (fs=2000 也可以)
T=1;%总时长(s)
t=0:T/fs:T;%时间序列向量
Del=t(2)-t(1);%时间间隔
N=length(t);%时间离散化数量
T_ext=[0.2, 0.4, 0.6, 0.8];%外部电流刺激时间点

% --- 【【【 这 是 唯 一 需 要 修 改 的 地 方 】】】 ---
%
% 1. 定义脉冲宽度
pulse_width = 1e-2; % 10ms
% 2. 定义幅度 (参考 论文2.pdf，使用 10 或 20，而不是 10000)
I_mag = 10; 
% 3. 使用【函数句柄 @rectpuls】来定义 I_ext
%    这一定义方式可以正确处理 RK4 传入的【标量】 t(i)
I_ext = @(u) I_mag * pulstran(u, T_ext, @rectpuls, pulse_width);
%
% --- 【【【 修 改 结 束 】】】 ---

figure('Name','外部电流刺激I_{ext}-t示意图') %[output:60898b3c]
plot(t,I_ext(t)) % 这样绘图仍然是正确的
hold on %[output:60898b3c]
xlabel('t/s'); %[output:60898b3c]
ylabel('外部电流刺激I_{ext}(\muA/cm^{2})') %[output:60898b3c]
title(sprintf('幅度 = %g', I_mag));
%[text] 
%[text] %[text:anchor:TMP_3bb7] PPT中给出的微分方程是
%[text] $&dollar&;&dollar&;\n\\begin{cases}\nC\_m \\frac{dV}{dt} = -g\_L(V-E\_L) - \\bar{g}\_{Na} m^3 h (V-E\_{Na}) - \\bar{g}\_K n^4 (V-E\_K) + I\_{app} \\\\\n\\frac{dm}{dt} = \\alpha\_m(V)(1-m) - \\beta\_m(V)m \\\\\n\\frac{dh}{dt} = \\alpha\_h(V)(1-h) - \\beta\_h(V)h \\\\\n\\frac{dn}{dt} = \\alpha\_n(V)(1-n) - \\beta\_n(V)n\n\\end{cases}\n&dollar&;&dollar&;$
% 初始化，增加速度
V=zeros(1,N);m=zeros(1,N);
h=zeros(1,N);n=zeros(1,N);

% 【 fV 函数现在会调用上面修正后的 I_ext 】
fV=@(V,m,h,n,t_val) 1/C*(-g.L*(V-E.L)-g.Na*m^3*h*(V-E.Na)-g.K*n^4*(V-E.K)+I_ext(t_val));
fm=@(V,m) Alpha.m(V)*(1-m)-Beta.m(V)*m;
fh=@(V,h) Alpha.h(V)*(1-h)-Beta.h(V)*h;
fn=@(V,n) Alpha.n(V)*(1-n)-Beta.n(V)*n;

% --- 初始值：使用【静息电位 -65mV】 ---
V0 = E.L; % -65 mV
m_inf = @(v) Alpha.m(v) ./ (Alpha.m(v) + Beta.m(v));
h_inf = @(v) Alpha.h(v) ./ (Alpha.h(v) + Beta.h(v));
n_inf = @(v) Alpha.n(v) ./ (Alpha.n(v) + Beta.n(v));

V(1) = V0;
m(1) = m_inf(V0);
h(1) = h_inf(V0);
n(1) = n_inf(V0);
% --- 初始值设置结束 ---

for i=1:N-1
    kV1=fV(V(i),m(i),h(i),n(i),t(i));
    km1=fm(V(i),m(i));
    kh1=fh(V(i),h(i));
    kn1=fn(V(i),n(i));

    kV2=fV(V(i)+Del/2*kV1,m(i)+Del/2*km1,h(i)+Del/2*kh1,n(i)+Del/2*kn1,t(i)+Del/2);
    km2=fm(V(i)+Del/2*kV1,m(i)+Del/2*km1);
    kh2=fh(V(i)+Del/2*kV1,h(i)+Del/2*kh1);
    kn2=fn(V(i)+Del/2*kV1,n(i)+Del/2*kn1);

    % kV3=Del*f(x(i)+Del/2,y(i)+1/2*kV2);
    kV3=fV(V(i)+Del/2*kV2,m(i)+Del/2*km2,h(i)+Del/2*kh2,n(i)+Del/2*kn2,t(i)+Del/2);
    km3=fm(V(i)+Del/2*kV2,m(i)+Del/2*km2);
    kh3=fh(V(i)+Del/2*kV2,h(i)+Del/2*kh2);
    kn3=fn(V(i)+Del/2*kV2,n(i)+Del/2*kn2);

    % kV4=Del*f(x(i)+Del,y(i)+kV3);
    kV4=fV(V(i)+Del*kV3,m(i)+Del*km3,h(i)+Del*kh3,n(i)+Del*kn3,t(i)+Del);
    km4=fm(V(i)+Del*kV3,m(i)+Del*km3);
    kh4=fh(V(i)+Del*kV3,h(i)+Del*kh3);
    kn4=fn(V(i)+Del*kV3,n(i)+Del*kn3);

    V(i+1)=V(i)+Del/6*(kV1+2*kV2+2*kV3+kV4);
    m(i+1)=m(i)+Del/6*(km1+2*km2+2*km3+km4);
    h(i+1)=h(i)+Del/6*(kh1+2*kh2+2*kh3+kh4);
    n(i+1)=n(i)+Del/6*(kn1+2*kn2+2*kn3+kn4);
end
figure %[output:5e6c6389]
plot(t,V) %[output:5e6c6389]
xlabel('t/s');
ylabel('膜电位 V (mV)');
title(sprintf('HH 模型 - 脉冲刺激 (幅度 = %g)', I_mag));
%[text] 
%%