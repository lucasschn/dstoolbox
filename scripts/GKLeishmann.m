clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')
addpath('plot_dir/')
%% Static data

data = load('naca0012');
data = data.naca0012;

naca0012 = Airfoil('naca0012',0.1);
M = 0.3;
a = 340.3;
nu = 1.48e-5;
Re = M*a*naca0012.c/nu;

% Look at OpenJetCorr if using simcos data
naca0012.steady = SteadyCurve(data.alpha_st,data.CN_st);

% pitching motion
mean_rad = deg2rad(data.mean);
amp_rad = deg2rad(data.amp);
pitching = PitchingMotion('alpha',data.alpha_xp,'CN',data.CN_xp,'V',M*a,'k',0.1);
pitching.setName();
pitching.setSinus(naca0012,mean_rad,amp_rad,2*pi/0.032,-pi/2);


tau1 = 0.01*naca0012.c/pitching.V;
tau2 = 1*naca0012.c/pitching.V;

pitching.GomanKhrabrov(naca0012.steady,tau1,tau2)

pitching.plotGK();