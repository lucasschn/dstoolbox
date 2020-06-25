% This script is test Sheng-LB model on SH2019
% Author : Lucas Schneeberger
% Date : 06.06.2020
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
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

%% Setting up the ramps

c = 22;

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
        load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','zero');
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        Cl = raw.Cl-zero.Cl;
        Cd = raw.Cd-zero.Cd;
    end
    fs = 1/ramp.Ts;
    Clf = myFilterTwice(Cl,fs);
    Cdf = myFilterTwice(Cd,fs);
%     ramp.setCL(Cl);
%     ramp.setCD(Cd);
    ramp.setCL(Clf);
    ramp.setCD(Cdf);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    ramp.setPitchRate(airfoil);
end

%% Run Sheng-LB model on all ramps

for k=1:length(c) 
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    evalin('base',sprintf('%s.BLSheng(airfoil,0,1,1,''experimental'')',msname))
    evalin('base',sprintf('%s.plotShengLB(airfoil)',msname))
end

