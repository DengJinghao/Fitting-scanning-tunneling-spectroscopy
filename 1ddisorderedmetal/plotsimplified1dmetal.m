syms x y a b t1 t;
x = -0.030:0.001:0.030;
y = [];
a = 0.005;
b = 1e-12;
Temp = 4.2.*0.083.*0.001;
fun2 = @(t,a,x)(sin(t.^2).*exp(-2.*t.*((a.*1i)./abs(x)).^(1/2))./t);
%���ֵķ��̲��֣�fun2�Ǽ򻯵ķ��̲��֣�a��ʾ�ŵ��ӹ���ǿ�ȣ�bֻ��һ����ͨ�ķŴ�ϵ����.*����������ʾ�����Ԫ��*Ԫ��
y = b.*real(integral(@(t)fun2(t,a,x),0,10,'ArrayValued',true));
hold on
plot (x,y,'LineWidth',2)
hold off