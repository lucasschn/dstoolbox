close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')


N = 100;
alpha = [zeros(1,N/2),ones(1,N/2)];
t = linspace(-.1,.1,N);

data = load('naca0012');
data = data.naca0012;

naca0012 = Airfoil('naca0012',0.1);
naca0012.steady = SteadyCurve(data.alpha_st,data.CN_st);
naca0012.steady.computeSlope();

pitching = PitchingMotion('alpha',alpha,'t',t,'V',50);
pitching.plotAlpha()
pitching.computeAttachedFlow(naca0012,'experimental')


figure
plot(t(1:length(pitching.CNI)),pitching.CNI,'DisplayName','C_N^I')
hold on 
plot(t(1:length(pitching.CNC)),pitching.CNC,'DisplayName','C_N^C')
plot(t(1:length(pitching.CNp)),pitching.CNp,'DisplayName','C_N^p')
xlabel('t (s)')
ylabel('C_N')
grid on
legend show
% axes('XAxisLocation','top','YAxisLocation','right','YTick',pitching.S)
% xlabel('S (-)')
title('Indicial response of attached flow coeffs')

