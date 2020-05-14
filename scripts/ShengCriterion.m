close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Calibration
calib = load('data/2018_SH_003/Postprocessing/matfiles/loads/ms014mpt001_loads.mat');
cL0 = mean(calib.raw.cL);
cD0 = mean(calib.raw.cD);

%% Paths to matfiles
path2ramps = 'data/2018_SH_003/Postprocessing/matfiles/loads/';

% avg = load('data/2019_SH/loads/ms013mpt001_loads.mat','avg');

matnames = {'ms018mpt001_loads.mat'};
matnames(end+1) = {'ms019mpt001_loads.mat'};
matnames(end+1) = {'ms020mpt001_loads.mat'};
matnames(end+1) = {'ms021mpt001_loads.mat'};
matnames(end+1) = {'ms022mpt001_loads.mat'};
matnames(end+1) = {'ms023mpt001_loads.mat'};
matnames(end+1) = {'ms024mpt001_loads.mat'};

%% Setting up the ramps
alphadot = [5 10 12.5 15 17.5 20 25];

for k=1:length(matnames)
    path = [path2ramps matnames{k}];
    runExp(path,k,alphadot(k),cL0,cD0)
end

%% Running Sheng experiment

airfoil = Airfoil('flatplate',0.15);
% Define alpha_ds0 & compute Talpha
airfoil.Sheng(ms018,ms019,ms020,ms021,ms022,ms023,ms024);

%% Add Sheng's predicted stall angles to the figures
for k=1:length(matnames)
    path = [path2ramps matnames{k}];
    msname = path(end-20:end-16);
    evalin('base',sprintf('figure(fig%d)',k))
    hold on
    evalin('base',sprintf('plot(%s.alpha_lagonset*ones(2,1),fig%d.CurrentAxes.YLim,''b--'')',msname,k));
end
%% Functions

function runExp(path,k,alphadot,cL0,cD0)
avg = load(path,'avg');
avg = avg.avg;
msname = path(end-20:end-16);
assignin('base',msname,RampUpMotion('alpha',avg.alphaa,'t',avg.t,'V',0.5));
evalin('base',sprintf('%s.setName()',msname))
ramp = evalin('base',msname);
ramp.setAlphaDot(alphadot) % in degrees
ramp.setCL(avg.cL - cL0);
ramp.setCD(avg.cD - cD0);
ramp.computeAirfoilFrame();
ramp.isolateRamp();
% Define stall
ramp.findExpOnset();
evalin('base',sprintf('fig%d = %s.plotCC()',k,msname));
end