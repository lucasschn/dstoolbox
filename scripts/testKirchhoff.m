close all
clear all
clc

airfoil = Airfoil('flatplate',0.15);

static = load('static_flatplate.mat');
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

airfoil.steady.fitKirchhoff()
airfoil.steady.plotKirchhoff()

% Optimal least square: S1=3.3572 S2=2.9679
plot(airfoil.steady.alpha,kirchhoff(airfoil.steady,airfoil.steady.alpha,[1.5 3]),'DisplayName','custom coeffs','LineWidth',2)
airfoil.steady.computeSeparation()
airfoil.steady.plotSeparation()

