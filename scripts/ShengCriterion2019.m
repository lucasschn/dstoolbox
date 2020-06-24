close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')
addpath('../src/model/')
addpath('../src/common/')
addpath('../src/lib/')

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load('../static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

%% Setting up the ramps

c = [18,14,22,67,26,84,30,89];

for k=1:length(c)    
    if LB(c(k)).ms >= 13 && LB(c(k)).ms < 100
        data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','inert','avg','zero');
        raw = data.raw;
        inert = data.inert;
        inert.alpha = raw.alpha(raw.t>=0);
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        Cl = inert.Cl;
        Cd = inert.Cd;
    else
        load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','zero');
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        Cl = raw.Cl-zero.Cl;
        Cd = raw.Cd-zero.Cd;
    end
    fs = 1/ramp.Ts;
    Cl_fff = myFilter(Cl,fs);
    Cd_fff = myFilter(Cd,fs);
    ramp.setCL(Cl_fff);
    ramp.setCD(Cd_fff);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    ramp.setPitchRate(airfoil);
    evalin('base',sprintf('fig%d = %s.plotCC(''convectime'');',k,msname));
end

%% Running Sheng experiment

% Define alpha_ds0 & compute Talpha
setLinFit(airfoil,ms012mpt1,ms010mpt1,ms013mpt1,ms025mpt1,ms014mpt1,ms034mpt1,ms015mpt1,ms116mpt1);

%saveas(gcf,'../fig/alpha_ds_r','png')
% %% Add Sheng's predicted stall angles to the figures
% for k=1:length(c)
%     msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
%     evalin('base',sprintf('figure(fig%d)',k))
%     hold on
%     evalin('base',sprintf('plot(%s.alpha_lagonset*ones(2,1),fig%d.CurrentAxes.YLim,''b--'')',msname,k));
% end