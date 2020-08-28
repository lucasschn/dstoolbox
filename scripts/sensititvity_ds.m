% This script aims at 
% Author : Lucas Schneeberger
% Date : 21.08.2020

close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,8);

%% Set up the ramp

c = 22;

data = load(loadmat(LB(c).ms,LB(c).mpt),'raw','inert','avg','zero');
raw = data.raw;
inert = data.inert;
inert.alpha = raw.alpha(raw.t>=0);
msname = sprintf('ms%03impt%i',LB(c).ms,LB(c).mpt);
ramp = RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c).U,'alphadot',LB(c).alphadot);
ramp_filt = RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c).U,'alphadot',LB(c).alphadot);
evalin('base',sprintf('ramp.setName(''%s'') ',msname))
evalin('base',sprintf('ramp_filt.setName(''%s'')',msname))

Cl = inert.Cl;
Cd = inert.Cd;
fs = 1/ramp.Ts;
Clf = myFilterTwice(Cl,fs);
Cdf = myFilterTwice(Cd,fs);
ramp.setCL(Clf);
ramp.setCD(Cdf);

ramp.computeAirfoilFrame();
ramp.isolateRamp();
ramp.setPitchRate(airfoil);
% Define stall (convectime must have been set)
ramp.findExpOnset();
%% Run Leishman-Beddoes' model on the ramp
Tp = 1:0.25:9;
CNds = -ones(size(Tp));
CNprime_ds = -ones(size(Tp));

for k=1:length(Tp)
    ramp.BeddoesLeishman(airfoil,Tp(k),0,2.5,ramp.estimateTvl,'experimental');
    [CNds(k),imaxCNLB] = max(ramp.CN_LB); % maximum of the model prediction
    [~,imaxCN] = max(ramp.CN); % maximum of the experimental result
    CNprime_ds(k) = ramp.CNprime(imaxCN);
end

figure 
plot(Tp,CNds,'d','MarkerFaceColor','b','DisplayName','C_{N,ds}')
hold on 
plot(Tp,CNprime_ds,'d','MarkerFaceColor','r','DisplayName','C''_{N,ds}')
yline(1.1457,'--','LineWidth',2,'DisplayName','C_N^{CRIT}')
xlabel('T_p')
ylabel('C_N')
grid on
ax = gca; 
ax.FontSize = 20; 
legend('Location','East')
