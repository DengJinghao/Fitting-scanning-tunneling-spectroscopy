clear variables;
tab = readtable('D:\STM data\LHE\20201120VCl3-NbSe2-He3\Bias-Spectroscopy00011.dat');
bias = tab.BiasCalc_V_;didv = tab.LIDemod1X_A_;current = tab.Current_A_;
x1 = 1000*bias; y1 = 1e12*didv;
%x = x + 0.00015;%symmetrization
x2 = x1(~excludedata(x1,y1,'domain',[1.1 1.8]));y2 = y1(~excludedata(x1,y1,'domain',[1.1 1.8]));
%��ȡ���ݣ�xΪbias(V) y Ϊdidv Ϊ�˹�һ����׼������ȡ��current��Ϊ����,x1
%y1Ϊ��Χѡ����domain��ʾ��x����exclude��excludedata���ص���logic0 1��~excludedata�����ų�ʵ�ʵ�����
syms nf energy r delta funbcs Temp ntip x;
%nf = 1; r = 0.1; delta = 5; 
%nfΪ������rΪ��չ��mev����deltaΪ����gap��meV��,delta =1.76kbt for bcs case
%fermi = (1./(exp((energy + bias)./(0.0863.*Temp))+1));
Temp = 0.8;
nf = [];delta = [];energy = [];r = [];x = [];
fun1 = @(x,energy)(0.25.*0.0863.*Temp.*(sech((energy-x)./(2.*0.0863.*Temp))).^2);%derivatation of fermi dirac function
%differmi = matlabFunction(differmi);
%y = integral(@(energy)(fun(Temp,bias,energy).*(nf.*real(abs(energy + 1i.*r)./((energy + 1i.*r).^2 - delta.^2).^0.5))),-1000,1000);
fun = @(nf,r,delta,x)(nf.*integral(@(energy)((fun1(x,energy)).*(real(abs(energy + 1i.*r)./((energy + 1i.*r).^2 - delta.^2).^0.5))),-30,30,'ArrayValued',true));
%�������Դ���biasȡֵ��Χ��quadgk���ø�˹���õ���⣬��integral��ǿ��sech�������Ͱ�߿�Ϊ10mV����ȫ���20mV
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0.1,0.01,0.1],'Upper',[10000,1,3],...
    'Diffminchange',1e-8,'Algorithm','Trust-Region','StartPoint',[100,0.1,1.2]...
    ,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',1e-15,'TolX',1e-20,'Robust','Bisquare');
%startpoint�ǳ��ؼ���Ҫѡȡ���ʵĳ�ֵ����������һ�·���Ԥ����ֵ
%maxfunevals�뷽�����㾫���йأ�maxiter���������������������tolfun tolx��ϵ���ͱ��������̶ȣ�algorithm�㷨
ft = fittype(fun,'options',fo);
[curve,gof] = fit(x2,y2,ft);%ֻѡȡĳ����Χ�����ݽ������;

hold on %������ͼ��ʱ������ͼ��ans
figure(1);
plot(x1,y1,':','LineWidth',2);%ɢ��ͼo��
plot(curve)%fit����
xlabel('Bias (v)')
ylabel('dI/dV (a.u.)')
%axis([-0.05 0.05 0 1e-12])
hold off
figure(2);
plot(curve,x2,y2,'residuals');