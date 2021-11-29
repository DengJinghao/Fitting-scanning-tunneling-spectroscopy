clear all;
syms nf energy r delta funbcs bias Temp ntip didv;
nf = 157.4; r = 0.08; delta =  1.353; %nf为常数，r为热展宽（mev），delta为超导gap（meV）
ntip = 2;
Temp = 0.4;
didv = [];
for bias = -5:0.01:5
funbcs = [];%density of states
funbcs = @(nf, r, delta, bias)(nf.*real(abs(bias + i.*r)./((bias + i.*r).^2 - delta.^2).^0.5));
%fermi = (1./(exp((energy + bias)./(0.0863.*Temp))+1));
fun = @(Temp,bias,energy)(0.25.*0.0863.*Temp.*(sech((energy-bias)./(2.*0.0863.*Temp))).^2);%derivatation of fermi dirac function
%differmi = matlabFunction(differmi);
%y = integral(@(energy)(fun(Temp,bias,energy).*(nf.*real(abs(energy + 1i.*r)./((energy + 1i.*r).^2 - delta.^2).^0.5))),-1000,1000);
y = quadgk(@(energy)((fun(Temp,bias,energy)).*(nf.*real(abs(energy + 1i.*r)./((energy + 1i.*r).^2 - delta.^2).^0.5))),-30,30);
%积分限略大于bias取值范围，quadgk采用高斯勒让德求解，比integral更强大，sech函数峰型半高宽为10mV，完全宽度20mV
didv = [didv,y];
end
didv = double(didv);
bias = -5:0.01:5;
plot (bias,didv,'Linewidth',2)