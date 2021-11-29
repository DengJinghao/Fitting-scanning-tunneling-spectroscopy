clear all;
syms nf energy r delta y;
nf = 157.4; r = 0.01; delta =  1.353; %nf为常数，r为热展宽（mev），delta为超导gap（meV）
energy = -5:0.01:5;%meV
y = [];%density of states
y = nf.*real(abs(energy + i.*r)./((energy + i.*r).^2 - delta.^2).^0.5);
hold on
figure(2);
plot (energy,y,'Linewidth',2)
hold off