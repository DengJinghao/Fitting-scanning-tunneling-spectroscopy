clear all;
tab = readtable('D:\STM data\LHE\20201120VCl3-NbSe2-He3\Bias-Spectroscopy00671.dat');
bias = tab.BiasCalc_V_;didv = tab.LIDemod1X_A_;current = tab.Current_A_;
x = bias; y = didv;
x = x + 0.00015;%symmetrization
x1 = x(~excludedata(x,y,'domain',[-0.01 0.01]));y1 = y(~excludedata(x,y,'domain',[-0.01 0.01]));
%��ȡ���ݣ�xΪbias(V) y Ϊdidv Ϊ�˹�һ����׼������ȡ��current��Ϊ����,x1
%y1Ϊ��Χѡ����domain��ʾ��x����exclude��excludedata���ص���logic0 1��~excludedata�����ų�ʵ�ʵ�����
syms nf energy r delta;
%nf = 1; r = 0.1; delta = 5; 
%nfΪ������rΪ��չ��mev����deltaΪ����gap��meV��,delta =1.76kbt
%energy = -20:0.1:20;%meV
%y = [];%density of states
fun = @(nf, r, delta, x)(nf.*real(abs(x + i.*r)./((x + i.*r).^2 - delta.^2).^0.5));
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[5e-14,0.00001,0.0001],'Upper',[1e-11,0.0001,0.002],...
    'Diffminchange',1e-20,'Algorithm','Trust-Region','StartPoint',[1e-13,0.0005,0.001]...
    ,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-50,'TolX',1e-50,'Robust','Bisquare');
%startpoint�ǳ��ؼ���Ҫѡȡ���ʵĳ�ֵ����������һ�·���Ԥ����ֵ
%maxfunevals�뷽�����㾫���йأ�maxiter���������������������tolfun tolx��ϵ���ͱ��������̶ȣ�algorithm�㷨
ft = fittype(fun,'options',fo);
[curve,gof] = fit(x1,y1,ft)%ֻѡȡĳ����Χ�����ݽ������;

hold on %������ͼ��ʱ������ͼ��ans
figure(3);
plot(x,y,':','LineWidth',2)%ɢ��ͼo��
plot(curve)%fit����
xlabel('Bias (v)')
ylabel('dI/dV (a.u.)')
%axis([-0.05 0.05 0 1e-12])
hold off
figure(4);
plot(curve,x1,y1,'residuals');
