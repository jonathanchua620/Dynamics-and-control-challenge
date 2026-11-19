open('MeasuredSignals-1.mat');

t = MeasuredSignals.t;

x=MeasuredSignals.x;

varphi=MeasuredSignals.varphi;

figure
plot(t, x);

figure
plot(t, varphi);

pppp