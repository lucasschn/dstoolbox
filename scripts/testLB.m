% This script is made to test Leishman-Beddoes model on SH2019 with
% user-defined parameters
% Author : Lucas Schneeberger
% Date : 11.06.2020

close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')

run(fullfile('..','labbook.m'))
run setPaths

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

ramp = loadRamp(c,true);
ramp.setPitchRate(airfoil);
% Define stall (convectime must have been set)
ramp.findExpOnset();
%% Run Leishman-Beddoes' model on the ramp
ramp.BeddoesLeishman(airfoil,0.5,2,2.5,3,'experimental')
ramp.plotCustom('CN','CN_LB')
saveas(gcf,fullfile(path2fig,'testFig'),'png')
