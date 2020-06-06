close all
clear all
clc

airfoil = Airfoil('flatplate',0.15);

static = load('static_flatplate.mat');
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

airfoil.steady.fitKirchhoff()
airfoil.steady.plotKirchhoff()

% Optimal least square: S1=2.57 S2=6.49
plot(airfoil.steady.alpha,kirchhoff(airfoil.steady,airfoil.steady.alpha,[1.5 4]))
airfoil.steady.computeSeparation()
airfoil.steady.plotSeparation()