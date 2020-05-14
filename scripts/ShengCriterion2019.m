close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Setting up the ramps

c = [18,14,22,67,26,84,30,89,34,71,38,93,42,75,46,97,50];

for k=1:length(c)
    data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','inert','avg','zero');
    raw = data.raw;
    if LB(c(k)).ms >= 13 && LB(c(k)).ms < 100
        inert = data.inert;
        inert.alpha = raw.alpha(raw.t>=0);
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',inert.alpha,'t',inert.t,'V',LB(c(k)).U));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        fc = 35;
        fs = 1/ramp.Ts;
        [b,a] = butter(5,fc/(fs/2));
        Cl_filtered = filter(b,a,inert.Cl);
        Cd_filtered = filter(b,a,inert.Cd);
    else
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        fc = 35;
        fs = 1/ramp.Ts;
        [b,a] = butter(5,fc/(fs/2));
        Cl_filtered = filter(b,a,raw.Cl);
        Cd_filtered = filter(b,a,raw.Cd);
    end
    ramp.setAlphaDot(LB(c(k)).alphadot) % in degrees
    Cl_ff = movmean(Cl_filtered,30);
    Cd_ff = movmean(Cd_filtered,30);
    fp = 1/3;
    [b,a] = cheby2(6,20,36*fp/(fs/2));
    Cl_fff = filter(b,a,Cl_ff);
    Cd_fff = filter(b,a,Cd_ff);
    ramp.setCL(Cl_fff);
    ramp.setCD(Cd_fff);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    evalin('base',sprintf('fig%d = %s.plotCC();',k,msname));
end


%% Running Sheng experiment
airfoil = Airfoil('flatplate',0.15);
static = load('data/static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.polarforces.CN,13);
% Define alpha_ds0 & compute Talpha
airfoil.Sheng(airfoil,ms012mpt1,ms010mpt1,ms013mpt1,ms025mpt1,ms014mpt1,ms034mpt1,ms015mpt1,ms116mpt1,...
    ms016mpt1,ms026mpt1,ms017mpt1,ms117mpt1,ms018mpt1,ms027mpt1,ms019mpt1,ms118mpt1,ms020mpt1);
saveas(gcf,'fig/Sheng/ShengSH2019_dsr.png')
% %% Add Sheng's predicted stall angles to the figures
% for k=1:length(c)
%     msname = sprintf('ms%03i',LB(c(k)).ms);
%     evalin('base',sprintf('figure(fig%d)',k))
%     hold on
%     evalin('base',sprintf('plot(%s.alpha_lagonset*ones(2,1),fig%d.CurrentAxes.YLim,''b--'')',msname,k));
% end
A = airfoil.A;
B = airfoil.B;
save('expfit_flatplate','A','B')
figure(airfoil.fig)