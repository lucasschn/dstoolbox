% This script is made to test Leishman-Beddoes model on SH2019 with
% user-defined parameters
% Author : Lucas Schneeberger
% Date : 11.06.2020

close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')

run(fullfile('..','labbook.m'))
try
    load('paths','path2fig')
catch
    open setPaths.m
    error('The path to the figure folder has not been set. Please set your paths in setPaths.m, run it and run this script again.')
end

path2static = fullfile('..','static_flatplate');

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(path2static);
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

%% Set up the ramp

c = 71;

data = load(loadmat(LB(c).ms,LB(c).mpt),'raw','zero');
raw = data.raw;
zero = data.zero;
msname = sprintf('ms%03impt%i',LB(c).ms,LB(c).mpt);
ramp = RampUpMotion('alpha',raw.alpha-raw.alpha(1),'t',raw.t,'V',LB(c).U,'alphadot',LB(c).alphadot);
evalin('base',sprintf('ramp.setName(''%s'') ',msname))

Cl = raw.Cl-zero.Cl;
Cd = raw.Cd-zero.Cd;
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
ramp.BeddoesLeishman(airfoil,1,0.65,0.05,3,'experimental')
% saveas(gcf,fullfile(path2fig,'vortex_development_Tp5'),'fig')
