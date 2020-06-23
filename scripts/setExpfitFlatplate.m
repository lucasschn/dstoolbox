close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')
addpath('../src/model/')
addpath('../src/common/')
addpath('../src/lib/')
if ~strcmp(pwd,'/Users/lucas/Documents/EPFL/PDM/scripts')
    cd('/Users/lucas/Documents/EPFL/PDM/scripts')
end
%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
static = load('../static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

%% Setting up the ramps

c = [18,14,22,67,26,84,30,89,34,71,38,93,42,75,46,97,50];

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
        data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw');
        raw = data.raw;
        msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
        assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot));
        evalin('base',sprintf('%s.setName()',msname))
        ramp = evalin('base',msname);
        Cl = raw.Cl;
        Cd = raw.Cd;
    end
    fs = 1/ramp.Ts;
    Clfff = myFilter(Cl,fs);
    Cdfff = myFilter(Cd,fs);
    ramp.setCL(Clfff);
    ramp.setCD(Cdfff);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    % Define stall
    ramp.findExpOnset();
    ramp.setAlphaDot(LB(c(k)).alphadot) % in degrees
    ramp.setPitchRate(airfoil);
end

%% Running Sheng experiment
% the exponential fit needs to be set
fig = setExpfit(airfoil,ms012mpt1,ms010mpt1,ms013mpt1,ms025mpt1,ms014mpt1,ms034mpt1,ms015mpt1,ms116mpt1,...
    ms016mpt1,ms026mpt1,ms017mpt1,ms117mpt1,ms018mpt1,ms027mpt1,ms019mpt1,ms118mpt1,ms020mpt1);
% what about Talpha and alpha_lag now?
r = -ones(size(c));
Talpha = -ones(size(c));
alpha_lag_ds = -ones(size(c));
for k=1:length(c)
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    evalin('base',sprintf('%s.setName()',msname))
    ramp = evalin('base',msname);
    r(k) = ramp.r;
    [alpha_lag_ds(k),Talpha(k)] = Expfit(airfoil,ramp);
end

alpha_ss = airfoil.steady.alpha_ss;

figure(fig)
subplot(311)
plot(r,alpha_lag_ds,'.','DisplayName','\alpha_{ds} (lagged)','MarkerSize',20)
plot(r,ones(size(r)).*alpha_ss,'--','DisplayName','\alpha_{ss}','LineWidth',2);
legend('FontSize',20,'Location','SouthEast')
% plot Talpha and its fit on the lower graph
subplot(312)
plot(r,Talpha,'.','MarkerSize',20,'DisplayName','T_\alpha')
hold on
xlabel('reduced pitch rate r (-)','FontSize',20);
ylabel('T_\alpha','FontSize',20)
grid on
saveas(gcf,'../fig/Sheng/ShengSH2019_dsr.png')

%% Add Sheng's predicted stall angles to the figures
% for k=1:length(c)
%     msname = sprintf('ms%03i',LB(c(k)).ms);
%     evalin('base',sprintf('figure(fig%d)',k))
%     hold on
%     evalin('base',sprintf('plot(%s.alpha_lagonset*ones(2,1),fig%d.CurrentAxes.YLim,''b--'')',msname,k));
% end


%% Fitting of Talpha(r)
% Talpha = @(x,r) 5*(1-exp(-r/x(1)).*cos(x(2)*r)); % cos without phase so that Talpha(0)=0
% xopt = lsqcurvefit(Talpha,[.04 80],r,obj.Talpha);
% plot(r,Talpha(xopt,r),'DisplayName','fit for T_\alpha')