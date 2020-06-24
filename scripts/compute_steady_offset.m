% This script is to measure the steady state offset for each pitch rate of
% SH2019 data
% Author : Lucas Schneeberger
% Date : 05.06.2020

close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
addpath('../plot_dir/')
addpath('../src/model/')
addpath('../src/common/')
addpath('../src/lib/')
run('/Users/lucas/src/codes_smarth/labbook.m')
%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load('../static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.CN,14);

%% Setting up the ramps

c = [22,67,26,84,30];

for k=1:length(c)
    if LB(c(k)).ms >= 13 && LB(c(k)).ms < 100
        data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','inert','avg','zero');
        raw = data.raw;
        inert = data.inert;
        inert.alpha = raw.alpha(raw.t>=0);
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        Cl = inert.Cl;
        Cd = inert.Cd;
    else
        data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw');
        raw = data.raw;
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        Cl = raw.Cl;
        Cd = raw.Cd;
    end
    fs = 1/ramp.Ts;
    Cl_fff = myFilter(Cl,fs);
    Cd_fff = myFilter(Cd,fs);
%     ramp.setCL(Cl);
%     ramp.setCD(Cd);
    ramp.setCL(Cl_fff);
    ramp.setCD(Cd_fff);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    ramp.setPitchRate(airfoil);
end

%% Run Beddoes Leishman model on all ramps

r = -ones(size(c));
steady_offset = -ones(size(c));

for k=1:length(c) 
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    evalin('base',sprintf('%s.BeddoesLeishman(airfoil,4.5,4,6,1,''experimental'')',msname))
    evalin('base',sprintf('%s.plotLB(''convectime'')',msname))
    ramp = evalin('base',msname);
    evalin('base',sprintf('saveas(gcf,''../fig/LB_r%03d'',''png'')',round(ramp.r,3)*1000))
    r(k) = evalin('base',sprintf('%s.r',msname));
    steady_offset(k) = evalin('base',sprintf('mean(%s.CN(%s.S > 0.9*%s.CN(end)))-%s.CN_LB(end)',msname,msname,msname,msname));
end

%Compute the static lift offset between experiment and Kirchhoff
kirchhoff_offset = airfoil.steady.CN(end)-kirchhoff(airfoil.steady,30);

% Plot and compare results
figure
plot(r,steady_offset,'.','MarkerSize',20,'Display','steady state offset')
%hold on 
%plot(r,kirchhoff_offset*ones(size(r)),'r-.','DisplayName','static offset due to Kirchhoff')
ax = gca; 
ax.FontSize = 20; 
xlabel('r')
ylabel('\Delta C_{N,\infty}')
grid on 


