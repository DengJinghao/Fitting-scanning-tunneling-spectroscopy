clear variables;
tab = readtable('D:\STM data\LHE\20201120VCl3-NbSe2-He3\Bias-Spectroscopy00011.dat');
bias = tab.BiasCalc_V_;didv = tab.LIDemod1X_A_;current = tab.Current_A_;
x1 = 1000*bias; y1 = 1e12*didv;
%x = x + 0.00015;%symmetrization
x2 = x1(~excludedata(x1,y1,'domain',[1.1 1.8]));y2 = y1(~excludedata(x1,y1,'domain',[1.1 1.8]));
%读取数据，x为bias(V) y 为didv 为了归一化做准备，读取了current作为备用,x1
%y1为范围选定，domain表示对x进行exclude，excludedata返回的是logic0 1，~excludedata用来排除实际的数据
syms nf energy r delta funbcs Temp ntip x;
%nf = 1; r = 0.1; delta = 5; 
%nf为常数，r为热展宽（mev），delta为超导gap（meV）,delta =1.76kbt for bcs case
%fermi = (1./(exp((energy + bias)./(0.0863.*Temp))+1));
Temp = 0.8;
nf = [];delta = [];energy = [];r = [];x = [];
fun1 = @(x,energy)(0.25.*0.0863.*Temp.*(sech((energy-x)./(2.*0.0863.*Temp))).^2);%derivatation of fermi dirac function
%differmi = matlabFunction(differmi);
%y = integral(@(energy)(fun(Temp,bias,energy).*(nf.*real(abs(energy + 1i.*r)./((energy + 1i.*r).^2 - delta.^2).^0.5))),-1000,1000);
fun = @(nf,r,delta,x)(nf.*integral(@(energy)((fun1(x,energy)).*(real(abs(energy + 1i.*r)./((energy + 1i.*r).^2 - delta.^2).^0.5))),-30,30,'ArrayValued',true));
%积分限略大于bias取值范围，quadgk采用高斯勒让德求解，比integral更强大，sech函数峰型半高宽为10mV，完全宽度20mV
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0.1,0.01,0.1],'Upper',[10000,1,3],...
    'Diffminchange',1e-8,'Algorithm','Trust-Region','StartPoint',[100,0.1,1.2]...
    ,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-15,'TolX',1e-20,'Robust','Bisquare');
%startpoint非常关键，要选取合适的初值，可以先跑一下方程预估初值，其余参数取决于数据
%maxfunevals与方程运算精度有关，maxiter是拟合至收敛的最大次数，tolfun tolx是系数和变量的容忍度，algorithm算法
ft = fittype(fun,'options',fo);
[curve,gof] = fit(x2,y2,ft);%只选取某个范围的数据进行拟合;

hold on %创建新图像时保留旧图像ans
figure(1);
plot(x1,y1,':','LineWidth',2);%散点图o型
plot(curve)%fit曲线
xlabel('Bias (v)')
ylabel('dI/dV (a.u.)')
%axis([-0.05 0.05 0 1e-12])
hold off
figure(2);
plot(curve,x2,y2,'residuals');
