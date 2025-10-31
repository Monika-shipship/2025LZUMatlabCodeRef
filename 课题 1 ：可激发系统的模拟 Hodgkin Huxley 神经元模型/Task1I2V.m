function [V,m,h,n]=Task1I2V(C,E,g,Alpha,Beta,options)



    % 输入：
    arguments
        C  %电容，可取1
        E %离子的能斯特电位结构体，如E.Na=115;E.K=-12;E.L=10.6;
        g %离子的电导结构体，如g.Na=120;g.K=36;g.L=0.3;

        Alpha %alpha函数结构体，如
        % Alpha.n=@(u) (0.1 - 0.01.*u) ./ (exp(1 - 0.1.*u) - 1);
        % Alpha.m=@(u) (2.5 - 0.1.*u) ./ (exp(2.5 - 0.1.*u) - 1);
        % Alpha.h=@(u) 0.07 .* exp(-u./20);
        Beta %beta函数结构体，如
        % Beta.n=@(u) 0.125 .* exp(-u./80);
        % Beta.m=@(u) 4 .* exp(-u./18);
        % Beta.h=@(u) 1 ./ (exp(3 - 0.1.*u) + 1);
        % 可选参数：
        options.fs = 4000 %采样频率(Hz)这一项将决定计算的精度
        options.T=200 %总时长(ms) 注意是毫秒
        % options.T_ext
        % options.pluWidth
        options.A =10 %刺激电流强度，单位是(\muA/cm^{2})
        options.I_ext = [] %外部刺激电流，匿名函数，单位是(\muA/cm^{2})
        options.Ifdebug = 1; %是否要详细的中间过程打印
        options.IfPic = 1; %是否要输出图片
    end
    Ifdebug=options.Ifdebug;
    IfPic=options.IfPic;
    %输出：
    % V
    % m
    % h
    % n

    fs= options.fs;%采样频率
    T=options.T;%总时长(ms)
    t=0:T/fs:T;%时间序列向量
    Del=t(2)-t(1);%时间间隔
    N=length(t);%时间离散化数量

    % T_ext=[75,120];%外部电流刺激，重复频率是5hz
    % FunExt=rectpuls(t,0.1);%外部电流刺激，间隔10ms
    % I_ext=@(u) 10*pulstran(u,T_ext,FunExt,fs);%脉冲刺激向量
    % pluWidth=50;%刺激电流宽度

    A=options.A;%刺激电流强度
    % I_ext = @(u) A * pulstran(u, T_ext, @rectpuls, pluWidth);

    %%脉冲刺激向量
    % 检查options.I_ext是否为空，即有没有输入外部的 I_ext，如果输入就用外部输入的，
    % 没输入就用默认的I_ext=@(u) 0+(20<=u & u<= 120).*options.A;
    if isempty(options.I_ext)
        I_ext=@(u) 0+(20<=u & u<= 120).*options.A;
    else
        I_ext=options.I_ext;
    end
    
    if Ifdebug==1

        fprintf("此时使用的参数为：\n C: %.4f \n",C)

        fprintf("g的结构体：\n")
        disp(g)

        fprintf("E的结构体：\n")
        disp(E);

        fprintf("Alpha函数的结构体：\n")
        disp(Alpha)

        fprintf("Beta函数的结构体：\n")
        disp(Beta)

        fprintf("外部刺激电流：\n")
        disp(I_ext)
    end
    if IfPic==1
    fprintf("\n下图是外部电流刺激的时域图\n")
    figure('Name','外部电流刺激I_{ext}-t示意图')
    plot(t,I_ext(t))
    hold on
    xlabel('t/ms');
    ylabel('外部电流刺激I_{ext}(\muA/cm^{2})')
    end

% PPT中给出的微分方程是

    % 初始化，增加速度
    V=zeros(1,N);m=zeros(1,N);
    h=zeros(1,N);n=zeros(1,N);

    fV=@(V,m,h,n,t) 1/C*(-g.L*(V-E.L)-g.Na*m^3*h*(V-E.Na)-g.K*n^4*(V-E.K)+I_ext(t));
    fm=@(V,m) Alpha.m(V)*(1-m)-Beta.m(V)*m;
    fh=@(V,h) Alpha.h(V)*(1-h)-Beta.h(V)*h;
    fn=@(V,n) Alpha.n(V)*(1-n)-Beta.n(V)*n;
    % 为了求解微分方程，我们还必须获得初始值，在静息时，I_ext肯定是0，V,m,h,n的导数也是零
    % 由此可以解方程数值求得初始值，同时也是静息值
    % 使用线性方程解出静息时的电位，各离子浓度比，此时I_ext肯定是0，所以得重新命名一个函数fVT来解方程
    syms V0 m0 h0 n0
    fVT=@(V,m,h,n) 1/C*(-g.L*(V-E.L)-g.Na*m^3*h*(V-E.Na)-g.K*n^4*(V-E.K)+0);
    eqt=[fVT(V0,m0,h0,n0),fm(V0,m0)==0,fh(V0,h0)==0,fn(V0,n0)==0];
    VPASol=vpasolve(eqt,[V0 m0 h0 n0]);
    V(1)=VPASol.V0;
    m(1)=VPASol.m0;
    h(1)=VPASol.h0;
    n(1)=VPASol.n0;

    %%自定义V(1)
    % V(1)=-53;
    % m(1)=Alpha.m(V(1))/(Alpha.m(V(1))+Beta.m(V(1)));
    % h(1)=Alpha.h(V(1))/(Alpha.h(V(1))+Beta.h(V(1)));
    % n(1)=Alpha.n(V(1))/(Alpha.n(V(1))+Beta.n(V(1)));

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
    if IfPic==1

        fprintf("下图是膜电压与时间的关系： \n")
        figure('Name','膜电压与时间的关系')
        plot(t,V)
        xlabel('t/ms')
        ylabel('V/mv')
    end
end