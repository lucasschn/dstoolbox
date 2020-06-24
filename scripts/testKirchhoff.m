close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')
addpath('../src/model/')
addpath('../src/common/')
addpath('../src/lib/')

airfoil = Airfoil('flatplate',0.15);

static = load('../static_flatplate.mat');
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

airfoil.steady.plotKirchhoff()
airfoil.steady.plotSeparation()

% % Optimal least square: S1=3.3572 S2=2.9679
% plot(airfoil.steady.alpha,kirchhoff(airfoil.steady,airfoil.steady.alpha,[0.1 3.2]),'DisplayName','custom coeffs','LineWidth',2)
% airfoil.steady.computeSeparation()
% airfoil.steady.plotSeparation()
% 
