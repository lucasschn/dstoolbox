close all
clear all
clc

airfoil = Airfoil('flatplate',0.15);

static = load('static_flatplate.mat');
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

airfoil.steady.fitKirchhoff()
airfoil.steady.plotKirchhoff()

plot(airfoil.steady.alpha,kirchhoff(airfoil.steady,airfoil.steady.alpha,[airfoil.steady.S1 0]))