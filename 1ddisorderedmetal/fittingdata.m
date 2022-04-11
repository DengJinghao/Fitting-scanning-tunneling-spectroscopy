tab = readtable('E://STM data1/LHE/20190923 WTE/Bias-Spectroscopy00179.dat');
bias = tab.BiasCalc_V_;didv = tab.LIDemod1X_AVG__A_;current = tab.Current_AVG__A_;
x = bias; y = didv;
x1 = x(~excludedata(x,y,'domain',[0.001 0.05]));y1 = y(~excludedata(x,y,'domain',[0.001 0.05]));
%读取数据，x为bias(V) y 为didv 为了归一化做准备，读取了current作为备用,x1
%y1为范围选定，domain表示对x进行exclude，excludedata返回的是logic0 1，~excludedata用来排除实际的数据
syms a b t t1;
Temp = 4.2.*0.083.*0.001;
fun = @(a,t,t1)(exp(-(a./pi).^0.5.*integral(@(t)(((1-cos(t.*t1))./(t.^1.5.*tanh(t./(2.*Temp))))),0,10,'ArrayValued',true)));
%integral(@(t)(((1-cos(t.*1))./(t.^1.5.*tanh(t./(2.*4.2))))),0,1000,'ArrayValued',true)=14.6433
%integral(@(t)(((1-cos(t.*1))./(t.^1.5.*tanh(t./(2.*4.2))))),0,inf,'ArrayValued',true)=14.7057
%in principle one SHOULD choose inf as the uplimit,for simplified choose 1000
fun2 = @(a,b,x)(b.*2.*Temp.*coth(abs(x)./(2.*Temp)).*integral(@(t1)(fun(a,t,t1).*(sin(abs(x).*t1).*cos((2.*a.*t1).^0.5))./(sinh(pi.*t1.*Temp))),0,100,'ArrayValued',true));
%fun represent the first integral function part,fun2 the second;
%a interaction strength(V);b coefficient;temp temperature(K)
%积分的时候@(t)表示对t积分，fun的（）部分t,a,x要把t积分的放在第一位，用作自变量的系数放在最后一位
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0.00001,-inf],'Upper',[0.1,inf],...
    'Diffminchange',1e-13,'Algorithm','Trust-Region','StartPoint',[0.001,3e-13]...
    ,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10e-20,'TolX',10e-20,'Robust','Bisquare');
%startpoint非常关键，要选取合适的初值，可以先跑一下方程预估初值
%maxfunevals与方程运算精度有关，maxiter是拟合至收敛的最大次数，tolfun tolx是系数和变量的容忍度，algorithm算法
ft = fittype(fun2,'options',fo);
[curve,gof] = fit(x1,y1,ft)%只选取某个范围的数据进行拟合;
hold on %创建新图像时保留旧图像ans
figure(2);
plot(x,y,':','LineWidth',2);%散点图o型
plot(curve);%fit曲线
xlabel('Bias (v)')
ylabel('dI/dV (a.u.)')
axis([-0.05 0.05 0 inf])
hold off
figure(3);
plot(curve,x1,y1,'residuals');