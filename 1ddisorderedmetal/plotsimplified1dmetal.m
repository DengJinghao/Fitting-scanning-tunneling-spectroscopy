syms x y a b t1 t;
x = -0.030:0.001:0.030;
y = [];
a = 0.005;
b = 1e-12;
Temp = 4.2.*0.083.*0.001;
fun2 = @(t,a,x)(sin(t.^2).*exp(-2.*t.*((a.*1i)./abs(x)).^(1/2))./t);
%积分的方程部分，fun2是简化的方程部分，a标示着电子关联强度，b只是一个普通的放大系数，.*带点的运算表示矩阵的元素*元素
y = b.*real(integral(@(t)fun2(t,a,x),0,10,'ArrayValued',true));
hold on
plot (x,y,'LineWidth',2)
hold off