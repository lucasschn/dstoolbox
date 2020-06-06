close all
clear all
clc
addpath('../src/common')
addpath('../src/lib')
set(0,'DefaultFigureWindowStyle','docked')
load('../static_flatplate.mat')

airfoil = Airfoil('flatplate',0.15);
airfoil.steady = SteadyCurve(alpha,CN,13);

airfoil.steady.fitKirchhoff()
airfoil.steady.plotKirchhoff()
airfoil.steady.plotSeparation()

airfoil.steady.plotViscousRatio()
hold on 
plot(airfoil.steady.fexp,((1+sqrt(airfoil.steady.fexp))/2).^2,'DisplayName','Kirchhoff model')
legend('Location','SouthEast')