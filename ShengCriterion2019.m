close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Setting up the ramps

c = [22:33,84:87];
ms = zeros(size(c));

for k=1:length(c)
    ms(k) = LB(c(k)).ms;
end

[ms_unique,index_unique] = unique(ms);

for k=1:length(ms_unique)
    for kk=1:4
        data = load(loadmat(LB(c(index_unique(k)+kk-1)).ms,LB(c(index_unique(k)+kk-1)).mpt),'raw','inert','avg','zero');
        inert = data.inert;
        raw = data.raw;
        tmax = 10;
        alpha = raw.alpha(inert.t<tmax);
        if kk == 1
            Clmean = inert.Cl(inert.t<tmax);
            Cdmean = inert.Cd(inert.t<tmax);
        else
            Clprev = Clmean;
            Cdprev = Cdmean;
            Clmean = (Clprev*(kk-1) + inert.Cl(inert.t<tmax))/kk;
            Cdmean = (Cdprev*(kk-1) + inert.Cd(inert.t<tmax))/kk;
        end
    end
    alpha = raw.alpha(raw.t>0);
    Clavg = inert.avgCl;
    Cdavg = inert.avgCd;
    msname = sprintf('ms%03i',ms_unique(k));
    assignin('base',msname,RampUpMotion('alpha',alpha,'t',inert.t,'V',LB(c(index_unique(k))).U));
    evalin('base',sprintf('%s.setName()',msname))
    ramp = evalin('base',msname);
    ramp.setAlphaDot(LB(c(index_unique(k))).alphadot) % in degrees
    ramp.setCL(Clavg);
    ramp.setCD(Cdavg);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    evalin('base',sprintf('fig%d = %s.plotCC();',k,msname));
end

%% Running Sheng experiment
airfoil = Airfoil('flatplate',0.15);
static = load('static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.Cn,13);
% Define alpha_ds0 & compute Talpha
figs = airfoil.Sheng(ms013,ms034,ms014,ms015);
saveas(gcf,'fig/Sheng/ShengSH2019_dsr.png')
%% Add Sheng's predicted stall angles to the figures
for k=1:length(ms_unique)
    msname = sprintf('ms%03i',LB(c(index_unique(k))).ms);
    evalin('base',sprintf('figure(fig%d)',k))
    hold on
    evalin('base',sprintf('plot(%s.alpha_lagonset*ones(2,1),fig%d.CurrentAxes.YLim,''b--'')',msname,k));
end

figure(figs)