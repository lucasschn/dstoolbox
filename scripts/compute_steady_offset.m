% This script is to measure the steady state offset for each pitch rate of
% SH2019 data
% Author : Lucas Schneeberger
% Date : 05.06.2020

close all
clear cll
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

c = [18,14,22,67,26,84,30,89];

for k=1:length(c)
    if LB(c(k)).ms >= 13 && LB(c(k)).ms < 100
        data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','inert','avg','zero');
        raw = data.raw;
        inert = data.inert;
        inert.alpha = raw.alpha(raw.t>=0);
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c(k)).U));
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
    % Butterworth filter
    fc = 35;
    fs = 1/ramp.Ts;
    [b,a] = butter(5,fc/(fs/2));
    Cl_filtered = filter(b,a,Cl);
    Cd_filtered = filter(b,a,Cd);
    ramp.setAlphaDot(LB(c(k)).alphadot) % in degrees
    % Moving average filter
    Cl_ff = movmean(Cl_filtered,30);
    Cd_ff = movmean(Cd_filtered,30);
    % Chebychev type-II filter
    fp = 1/3;
    [b,a] = cheby2(6,20,36*fp/(fs/2));
    Cl_fff = filter(b,a,Cl_ff);
    Cd_fff = filter(b,a,Cd_ff);
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
    evalin('base',sprintf('%s.BeddoesLeishman(airfoil,1.7,6,3,''experimental'')',msname))
    evalin('base',sprintf('%s.plotLB(''convectime'')',msname))
    r(k) = evalin('base',sprintf('%s.r',msname));
    steady_offset(k) = evalin('base',sprintf('%s.CN(end)-%s.CN_LB(end)',msname,msname));
end

%Compute the static lift offset between experiment and Kirchhoff
kirchhoff_offset = airfoil.steady.CN(end)-kirchhoff(airfoil.steady,30);

% Plot and compare results
figure
plot(r,steady_offset,'.','MarkerSize',20,'Display','steady state offset')
hold on 
plot(r,kirchhoff_offset*ones(size(r)),'r-.','DisplayName','static offset due to Kirchhoff')
xlabel('r')
ylabel('steady offset')
grid on 

