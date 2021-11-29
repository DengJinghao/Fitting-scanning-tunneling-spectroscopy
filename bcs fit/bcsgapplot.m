clear all;
syms nf energy r delta y;
nf = 157.4; r = 0.01; delta =  1.353; %nfΪ������rΪ��չ��mev����deltaΪ����gap��meV��
energy = -5:0.01:5;%meV
y = [];%density of states
y = nf.*real(abs(energy + i.*r)./((energy + i.*r).^2 - delta.^2).^0.5);
hold on
figure(2);
plot (energy,y,'Linewidth',2)
hold off