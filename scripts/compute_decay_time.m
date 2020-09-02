% This script is to measure the decay time from dynamic stall to
% steady-state. Steady-state is defined as when one of the secondary
% vortices drops pass under the finalvalue of LB prediction.
% Author : Lucas Schneeberger
% Date : 03.08.2020

close all
clear all
clc 
set(0,'DefaultFigureWindowStyle','docked')
addpath(fullfile('..','plot_dir'))
addpath(genpath(fullfile('..','src')))
run(fullfile('/Users','lucas','src','codes_smarth','labbook.m'))
%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

%% Setting up the ramps

c = [22,67,26,84,30,34,38];

for k=1:length(c)
    data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','zero');
    raw = data.raw;
    zero = data.zero;
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U));
    evalin('base',sprintf('%s.setName()',msname))
    ramp = evalin('base',msname);
    Cl = raw.Cl - mean(raw.Cl(1:50));
    Cd = raw.Cd - mean(raw.Cd(1:50));
    fs = 1/ramp.Ts;
    Cl_fff = myFilterTwice(Cl,fs);
    Cd_fff = myFilterTwice(Cd,fs);
%     ramp.setCL(Cl);
%     ramp.setCD(Cd);
    ramp.setCL(Cl_fff);
    ramp.setCD(Cd_fff);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();    
    ramp.setPitchRate(airfoil);
    % Define stall
    ramp.findExpOnset();
end

%% Run Beddoes-Leishman on all ramps

r = -ones(size(c));
tc_ds = -ones(size(c));
tc_inf = -ones(size(c));

for k=1:length(c) 
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    evalin('base',sprintf('%s.BeddoesLeishman(airfoil,4.5,4,6,1,''experimental'')',msname))
    evalin('base',sprintf('%s.plotLB(''convectime'')',msname))
    ramp = evalin('base',msname);    
    r(k) = evalin('base',sprintf('%s.r',msname));
    tc_ds(k) = ramp.S(ramp.i_CLonset);
    % defines steady-state as when on of the secondary vortices drops passes
    % under model-predicted value
    i_inf = find(ramp.CN(ramp.i_CLonset:end)<ramp.CN_LB(end),1);
    if ~isempty(i_inf)
        tc_inf(k) = ramp.S(ramp.i_CLonset + i_inf);
    else
        tc_inf(k) = ramp.S(end);
    end
end

delta_tc = tc_inf - tc_ds;

figure
plot(r,delta_tc,'d','MarkerFaceColor','b')
grid on 
xlabel('r')
ylabel('\Delta t_c')

figure
plot(r,tc_ds,'d','MarkerFaceColor','b')
grid on 
xlabel('r')
ylabel('t_{c,ds}')

%% Plot experiments on same figure

figure
for k=1:length(c)
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    ramp = evalin('base',msname);
    plot(ramp.S-tc_ds(k),ramp.CL/ramp.CL(ramp.i_CLonset),'LineWidth',1,'DisplayName',sprintf('%.2f °/s',ramp.alphadot))
    hold on 
end
legend show 
grid on 
xlabel('t_c')
ylabel('C_L')

maxCN = -ones(size(c));

figure
for k=1:length(c)
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    ramp = evalin('base',msname);
    maxCN(k) = max(ramp.CN);
    plot(ramp.S-tc_ds(k),ramp.CN,'LineWidth',1,'DisplayName',sprintf('%.2f °/s',ramp.alphadot))
    hold on 
end
legend show 
grid on 
xlabel('t_c')
ylabel('C_N')

figure 
plot(r,maxCN,'d','MarkerFaceColor','k')
grid on
xlabel('r')
ylabel('max C_N')