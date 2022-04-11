tab = readtable('E://STM data1/LHE/20190923 WTE/Bias-Spectroscopy00179.dat');
bias = tab.BiasCalc_V_;didv = tab.LIDemod1X_AVG__A_;current = tab.Current_AVG__A_;
x = bias; y = didv;
x1 = x(~excludedata(x,y,'domain',[0.001 0.05]));y1 = y(~excludedata(x,y,'domain',[0.001 0.05]));
%��ȡ���ݣ�xΪbias(V) y Ϊdidv Ϊ�˹�һ����׼������ȡ��current��Ϊ����,x1
%y1Ϊ��Χѡ����domain��ʾ��x����exclude��excludedata���ص���logic0 1��~excludedata�����ų�ʵ�ʵ�����
syms a b t t1;
Temp = 4.2.*0.083.*0.001;
fun = @(a,t,t1)(exp(-(a./pi).^0.5.*integral(@(t)(((1-cos(t.*t1))./(t.^1.5.*tanh(t./(2.*Temp))))),0,10,'ArrayValued',true)));
%integral(@(t)(((1-cos(t.*1))./(t.^1.5.*tanh(t./(2.*4.2))))),0,1000,'ArrayValued',true)=14.6433
%integral(@(t)(((1-cos(t.*1))./(t.^1.5.*tanh(t./(2.*4.2))))),0,inf,'ArrayValued',true)=14.7057
%in principle one SHOULD choose inf as the uplimit,for simplified choose 1000
fun2 = @(a,b,x)(b.*2.*Temp.*coth(abs(x)./(2.*Temp)).*integral(@(t1)(fun(a,t,t1).*(sin(abs(x).*t1).*cos((2.*a.*t1).^0.5))./(sinh(pi.*t1.*Temp))),0,100,'ArrayValued',true));
%fun represent the first integral function part,fun2 the second;
%a interaction strength(V);b coefficient;temp temperature(K)
%���ֵ�ʱ��@(t)��ʾ��t���֣�fun�ģ�������t,a,xҪ��t���ֵķ��ڵ�һλ�������Ա�����ϵ���������һλ
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0.00001,-inf],'Upper',[0.1,inf],...
    'Diffminchange',1e-13,'Algorithm','Trust-Region','StartPoint',[0.001,3e-13]...
    ,'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10e-20,'TolX',10e-20,'Robust','Bisquare');
%startpoint�ǳ��ؼ���Ҫѡȡ���ʵĳ�ֵ����������һ�·���Ԥ����ֵ
%maxfunevals�뷽�����㾫���йأ�maxiter���������������������tolfun tolx��ϵ���ͱ��������̶ȣ�algorithm�㷨
ft = fittype(fun2,'options',fo);
[curve,gof] = fit(x1,y1,ft)%ֻѡȡĳ����Χ�����ݽ������;
hold on %������ͼ��ʱ������ͼ��ans
figure(2);
plot(x,y,':','LineWidth',2);%ɢ��ͼo��
plot(curve);%fit����
xlabel('Bias (v)')
ylabel('dI/dV (a.u.)')
axis([-0.05 0.05 0 inf])
hold off
figure(3);
plot(curve,x1,y1,'residuals');