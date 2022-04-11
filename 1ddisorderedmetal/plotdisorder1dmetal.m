clear all;
syms x y a b t1 t;
x = -0.050:0.0002:0.050;
y = [];
a = 0.01;
b = 1e-12;
Temp = 4.2.*0.083.*0.001;
fun = @(a,t,t1)(exp(-(a./pi).^0.5.*integral(@(t)(((1-cos(t.*t1))./(t.^1.5.*tanh(t./(2.*Temp))))),0,10,'ArrayValued',true)));
%integral(@(t)(((1-cos(t.*1))./(t.^1.5.*tanh(t./(2.*4.2))))),0,1000,'ArrayValued',true)=14.6433
%integral(@(t)(((1-cos(t.*1))./(t.^1.5.*tanh(t./(2.*4.2))))),0,inf,'ArrayValued',true)=14.7057
%in principle one SHOULD choose inf as the uplimit,for simplified choose 1000
y = b.*2.*Temp.*coth(abs(x)./(2.*Temp)).*integral(@(t1)(fun(a,t,t1).*(sin(abs(x).*t1).*cos((2.*a.*t1).^0.5))./(sinh(pi.*t1.*Temp))),0,100,'ArrayValued',true);
hold on
figure(4);
plot (x,y,'LineWidth',2)
hold off