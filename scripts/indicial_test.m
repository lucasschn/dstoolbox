close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')


N = 1001;
if mod(N,2)==0
    alpha = [zeros(1,N/2),ones(1,N/2)];
else
    alpha = [zeros(1,(N-1)/2),ones(1,(N-1)/2+1)];
end
t = linspace(-1e-3,1e-3,N);

data = load('naca0012');
data = data.naca0012;

naca0012 = Airfoil('naca0012',0.1);
naca0012.steady = SteadyCurve(data.alpha_st,data.CN_st);
naca0012.steady.computeSlope();

pitching = PitchingMotion('alpha',alpha,'t',t,'V',50);
pitching.plotAlpha()
pitching.computeAttachedFlow(naca0012,'experimental')


figure
subplot(211)
p1 = plot(pitching.S(1:length(pitching.CNp)),pitching.CNp,'DisplayName','C_N^p','LineWidth',2);
hold on 
p2 = plot(pitching.S(1:length(pitching.CNI)),pitching.CNI,'DisplayName','C_N^I','LineWidth',2);
p3 = plot(pitching.S(1:length(pitching.CNC)),pitching.CNC,'DisplayName','C_N^C','LineWidth',2);

p1.Color(4) = 0.7;
p2.Color(4) = 0.7;
p3.Color(4) = 0.7;
ax = gca;
ax.FontSize = 20;
ylabel('C_N','FontSize',20)
grid on
legend('FontSize',20)
% axes('XAxisLocation','top','YAxisLocation','right','YTick',pitching.S)
% xlabel('S (-)')
title(sprintf('Indicial response of attached flow coeffs ( N = %d)',N))

subplot(212)
plot(pitching.S,pitching.alpha)
xlabel('S')
ylabel('\alpha','FontSize',20)
grid on
ax = gca;
ax.FontSize = 20;


fprintf('4/M = %.2f and phi(0)= %.2f should be equal.\n',4/pitching.M,pitching.CNp(pitching.S==0))

fprintf('2pi/beta = %.4f and phi(inf) = %.4f should be equal.\n',pi/180*2*pi/pitching.beta,pitching.CNC(end))