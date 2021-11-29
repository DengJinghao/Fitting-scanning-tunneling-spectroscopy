clear all;
tab = readtable('D:\STM data\LHE\20201120VCl3-NbSe2-He3\Bias-Spectroscopy00671.dat');
bias = tab.BiasCalc_V_;didv = tab.LIDemod1X_A_;current = tab.Current_A_;
x = bias; y = didv;
x = x + 0.00015;%symmetrization
x1 = x(~excludedata(x,y,'domain',[-0.01 0.01]));y1 = y(~excludedata(x,y,'domain',[-0.01 0.01]));
%读取数据，x为bias(V) y 为didv 为了归一化做准备，读取了current作为备用,x1
%y1为范围选定，domain表示对x进行exclude，excludedata返回的是logic0 1，~excludedata用来排除实际的数据
syms nf energy r delta;
%nf = 1; r = 0.1; delta = 5; 
%nf为常数，r为热展宽（mev），delta为超导gap（meV）,delta =1.76kbt
%energy = -20:0.1:20;%meV
%y = [];%density of states
fun = @(nf, r, delta, x)(nf.*real(abs(x + i.*r)./((x + i.*r).^2 - delta.^2).^0.5));
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[5e-14,0.00001,0.0001],'Upper',[1e-11,0.0001,0.002],...
    'Diffminchange',1e-20,'Algorithm','Trust-Region','StartPoint',[1e-13,0.0005,0.001]...
    ,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-50,'TolX',1e-50,'Robust','Bisquare');
%startpoint非常关键，要选取合适的初值，可以先跑一下方程预估初值
%maxfunevals与方程运算精度有关，maxiter是拟合至收敛的最大次数，tolfun tolx是系数和变量的容忍度，algorithm算法
ft = fittype(fun,'options',fo);
[curve,gof] = fit(x1,y1,ft)%只选取某个范围的数据进行拟合;

hold on %创建新图像时保留旧图像ans
figure(3);
plot(x,y,':','LineWidth',2)%散点图o型
plot(curve)%fit曲线
xlabel('Bias (v)')
ylabel('dI/dV (a.u.)')
%axis([-0.05 0.05 0 1e-12])
hold off
figure(4);
plot(curve,x1,y1,'residuals');
