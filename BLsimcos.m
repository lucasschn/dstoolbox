close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('data/2008_simcos/matlab/labbook.m')

%% Static data

static_data = csvread('fatma_data.csv');
airfoil = Airfoil('simcos',0.3);
airfoil.steady = SteadyCurve(static_data(:,1),static_data(:,2));

%% Dynamic data
nr=13;

data = load(pressuredata(nr));

pitching = PitchingMotion('alpha',data.alpha,'CL',data.mCl);
pitching.setSinus(airfoil,deg2rad(LB(nr).alpha_0),deg2rad(LB(nr).alpha_1),LB(nr).fosc,LB(nr).FS);
% pitching.computeAirfoilFrame();
% pitching.setCNsteady(airfoil.steady)
% 
% airfoil.steady.fitKirchhoff();
% airfoil.steady.plotKirchhoff();