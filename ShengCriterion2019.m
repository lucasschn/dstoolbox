close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Setting up the ramps

c = [22,26,84,30];

for k=1:length(c)
    data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','inert','avg','zero');
    avg = data.avg;
    zero = data.zero;
    raw = data.raw;
    inert = data.inert;
    inert.alpha = raw.alpha(raw.t>=0);
    msname = sprintf('ms%03i',LB(c(k)).ms);
    assignin('base',msname,RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c(k)).U));
    evalin('base',sprintf('%s.setName()',msname))
    ramp = evalin('base',msname);
    ramp.setAlphaDot(LB(c(k)).alphadot) % in degrees
    ramp.setCL(inert.Cl);
    ramp.setCD(inert.Cd);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    evalin('base',sprintf('fig%d = %s.plotCC();',k,msname));
end

%% Running Sheng experiment
airfoil = Airfoil('flatplate',0.15);
static = load('static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.Cn,14);
% Define alpha_ds0 & compute Talpha
figs = airfoil.Sheng(ms013,ms014,ms034,ms015);
saveas(gcf,'fig/Sheng/ShengSH2019_dsr.png')
%% Add Sheng's predicted stall angles to the figures
for k=1:length(c)
    msname = sprintf('ms%03i',LB(c(k)).ms);
    evalin('base',sprintf('figure(fig%d)',k))
    hold on
    evalin('base',sprintf('plot(%s.alpha_lagonset*ones(2,1),fig%d.CurrentAxes.YLim,''b--'')',msname,k));
end

figure(figs)