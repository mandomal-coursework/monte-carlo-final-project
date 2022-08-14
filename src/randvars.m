x = linspace(600d3,0,1000);



ionp = (0.5)*exp(-x/80d3);
excp = (0.5)*exp(-x/80d3);
elap = 1-(ionp+excp);


figure(3);clf;hold on; grid on;
plot(ionp,x,'linewidth',4)
plot(excp,x,'--','linewidth',4)
plot(elap,x,'linewidth',4)
xlabel 'Probability P'
ylabel 'Height'
title 'Height vs. Interaction Probability'
legend('P(ionize)','P(excite/dissociate)','P(elastic)')
saveeps('probability.png')

%%

r = rand(1,1000);
dte = -log(r)*5;
figure(4)
histogram(abs(dte),100,'normalization','probability')
ylabel 'Probability'
xlabel 'Distance before collision [m]'
title 'Distance before collision'
saveeps('mfp.png')

%%

en1 = rand(1,1000)*50;
en2 = rand(1,1000)*30;

figure(5)
subplot(121)
histogram(en1,100,'normalization','probability')
ylabel 'Probability'
xlabel 'Ionization energy [eV]'
subplot(122)
histogram(en2,100,'normalization','probability')
xlabel 'Excitation/Dissociation energy [eV]'
saveeps('endist.png')