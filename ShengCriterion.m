close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Import test data for a specific airfoil (Low-speed, Ma<0.3)

% alphadot = deg2rad([5,10,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]);
% % 700rpm is 0.5m/s
% V = [0.175,0.25,0.375,0.5,0.6125,0.75];
% 
% Re = V*c/8.9e-7;
% 
% [alphadotgrid,Vgrid] = meshgrid(V,alphadot);
% 
% r = alphadotgrid*c./(2*Vgrid);

calib = load('data/2018_SH_003/Postprocessing/matfiles/loads/ms001mpt001_loads.mat');
cL0 = mean(calib.raw.cL);
cD0 = mean(calib.raw.cD);

path2ramps = 'data/2018_SH_003/Postprocessing/matfiles/loads/';

% avg = load('data/2019_SH/loads/ms013mpt001_loads.mat','avg');
% avg = avg.avg;
% ms013 = RampUpMotion('ms013',avg.alpha);
% ms013.sett(avg.t);
% ms013.V = 0.5; %m/s
% ms013.setCL( avg.Cl);
% ms013.setCD( avg.Cd);
% ms013.computeAirfoilFrame();
% % ms013.isolateRamp();
% ms013.findExpOnset();
% fig1 = ms013.plotCC();

% ideal constructor 
%ms018 = RampUpMotion('ms019'

avg = load([path2ramps 'ms018mpt001_loads.mat'],'avg');
avg = avg.avg;
ms018 = RampUpMotion('ms018',avg.alphaa);
ms018.sett(avg.t);
ms018.V = 0.5; %m/s
ms018.setCL(avg.cL - cL0);
ms018.setCD(avg.cD - cD0);
ms018.computeAirfoilFrame();
ms018.isolateRamp();
% Define stall
ms018.findExpOnset();
fig0 = ms018.plotCC();
ms018.setalphadot(5)

avg = load([path2ramps 'ms019mpt001_loads.mat'],'avg');
avg = avg.avg;
ms019 = RampUpMotion('ms019',avg.alphaa);
ms019.sett(avg.t);
ms019.V = 0.5; %m/s
ms019.setCL(avg.cL - cL0);
ms019.setCD(avg.cD - cD0);
ms019.computeAirfoilFrame();
ms019.isolateRamp();
% Define stall
ms019.findExpOnset();
fig1 = ms019.plotCC();
ms019.setalphadot(10)


avg = load([path2ramps 'ms020mpt001_loads.mat'],'avg');
avg = avg.avg;
ms020 = RampUpMotion('ms020',avg.alphaa);
ms020.sett(avg.t);
ms020.V = 0.5; %m/s
ms020.setCL(avg.cL - cL0);
ms020.setCD(avg.cD - cD0);
ms020.computeAirfoilFrame();
ms020.isolateRamp();
% Define stall
ms020.findExpOnset();
fig2 = ms020.plotCC();
ms020.setalphadot(12.5)

avg = load([path2ramps 'ms021mpt001_loads.mat'],'avg');
avg = avg.avg;
ms021 = RampUpMotion('ms021',avg.alphaa);
ms021.sett(avg.t);
ms021.V = 0.5; %m/s
ms021.setCL(avg.cL - cL0);
ms021.setCD(avg.cD - cD0);
ms021.computeAirfoilFrame();
ms021.isolateRamp();
% Define stall
ms021.findExpOnset();
fig3 = ms021.plotCC();
ms021.setalphadot(15)

avg = load([path2ramps 'ms022mpt001_loads.mat'],'avg');
avg = avg.avg;
ms022 = RampUpMotion('ms022',avg.alphaa);
ms022.sett(avg.t);
ms022.V = 0.5; %m/s
ms022.setCL(avg.cL - cL0);
ms022.setCD(avg.cD - cD0);
ms022.computeAirfoilFrame();
ms022.isolateRamp();
% Define stall
ms022.findExpOnset();
fig4 = ms022.plotCC();
ms022.setalphadot(17.5)

avg = load([path2ramps 'ms023mpt001_loads.mat'],'avg');
avg = avg.avg;
ms023 = RampUpMotion('ms023',avg.alphaa);
ms023.sett(avg.t);
ms023.V = 0.5; %m/s
ms023.setCL(avg.cL - cL0);
ms023.setCD(avg.cD - cD0);
ms023.computeAirfoilFrame();
ms023.isolateRamp();
% Define stall
ms023.findExpOnset();
fig5 = ms023.plotCC();
ms023.setalphadot(20)

avg = load([path2ramps 'ms024mpt001_loads.mat'],'avg');
avg = avg.avg;
ms024 = RampUpMotion('ms024',avg.alphaa);
ms024.sett(avg.t);
ms024.V = 0.5; %m/s
ms024.setCL(avg.cL - cL0);
ms024.setCD(avg.cD - cD0);
ms024.computeAirfoilFrame();
ms024.isolateRamp();
% Define stall
ms024.findExpOnset();
fig6 = ms024.plotCC();
ms024.setalphadot(25)

%% Compute Talpha

airfoil = Airfoil('flatplate',0.15);
% Define alpha_ds0 & compute Talpha
airfoil.Sheng(ms019,ms020,ms021,ms022,ms023,ms024);

figure(fig1)
hold on
plot(ms019.alpha_Shengonset*ones(2,1),fig1.CurrentAxes.YLim,'b--');

figure(fig2)
hold on
plot(ms020.alpha_Shengonset*ones(2,1),fig2.CurrentAxes.YLim,'b--');

figure(fig3)
hold on
plot(ms021.alpha_Shengonset*ones(2,1),fig3.CurrentAxes.YLim,'b--');

figure(fig4)
hold on
plot(ms022.alpha_Shengonset*ones(2,1),fig4.CurrentAxes.YLim,'b--');

figure(fig5)
hold on
plot(ms023.alpha_Shengonset*ones(2,1),fig5.CurrentAxes.YLim,'b--');

figure(fig6)
hold on
plot(ms024.alpha_Shengonset*ones(2,1),fig6.CurrentAxes.YLim,'b--');

%% Import validation data for this airfoil

%% Compute alpha_prime

%% Compute stall onset

%% Save figures