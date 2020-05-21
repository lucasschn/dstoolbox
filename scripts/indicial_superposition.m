close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')


run('/Users/lucas/src/codes_smarth/labbook.m')

%% Setting up the ramps

% c = 26;
% 
% data = load(loadmat(LB(c).ms,LB(c).mpt),'raw','inert','avg','zero');
% avg = data.avg;
% zero = data.zero;
% raw = data.raw;
% inert = data.inert;
% inert.alpha = raw.alpha(raw.t>=0);
% msname = sprintf('ms%03i',LB(c).ms);
% assignin('base',msname,RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c).U));
% evalin('base',sprintf('%s.setName()',msname))
% ramp = evalin('base',msname);
% ramp.setAlphaDot(LB(c).alphadot) % in degrees
% ramp.setCL(inert.Cl);
% ramp.setCD(inert.Cd);
% ramp.computeAirfoilFrame();
% ramp.isolateRamp();
% 
% ramp.plotAlpha()
% 
% diff(ramp.alpha)
figure
subplot(211)
alpha = [0 0 1 2 3 3];
dalpha = diff(alpha);
t = 0:length(alpha)-1;
plot(t,alpha)
xlabel('t (s)')
ylabel('\alpha (°)')
grid on
subplot(212)
plot(t(2:end),dalpha,'x','MarkerSize',20)
xlabel('t (s)')
ylabel('\Delta\alpha (°)')
grid on

dalpha_lag1 = dalpha(2)*(1-exp(-t));
t1 = t(2:end);
figure
plot(t1,dalpha_lag1(1:end-1))

dalpha_lag2 = dalpha(3)*(1-exp(-t));
t2 = t(3:end);
figure
plot(t2,dalpha_lag2(1:end-2))

dalpha_lag = dalpha_lag1 + dalpha_lag2;
plot(t,dalpha