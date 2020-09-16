% This script is to measure the time from passing CNcrit defined as the
% limit of CN,ds as r->0 and CN,ds for each r
% Author : Lucas Schneeberger
% Date : 20.08.2020
close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')

run(fullfile('..','labbook.m'))

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

%% Look at the time lag between t_crit and t_ds

c = [18,14,22,67,26,84,30,2,34,71,38,6,42,75,46,10,50,62];

% initialization of the variables and vectors
load('../CNcrit','CNcrit')
% CNcrit = 0.9;
r = -ones(size(c));
tc_crit = -ones(size(c));
tc_ds = -ones(size(c));
ramps = cell(size(c));

for k = 1:length(c)   
    
    ramp = loadRamp(c(k));
    
    % Define reduced pitch rate if necessary ...
    if isempty(ramp.r)
        % compute it with alphadot
        ramp.setPitchRate(airfoil);
    end
    % ... and assign it
    r(k) = ramp.r;
    % Define stall
    ramp.findExpOnset();   
    ramps{k} = ramp;
    i_crit = find(ramp.CN > CNcrit,1);
    tc_crit(k) = ramp.S(i_crit);
    tc_ds(k) = ramp.S(ramp.i_CLonset);
end

%% Plot it as a function of the pitch rate
delta_crit = tc_ds-tc_crit;

t_crit = polyfit(r(r>0.01),delta_crit(r>0.01),0);

figure 
plot(r,delta_crit,'d','MarkerFaceColor','k')
hold on 
plot(r,t_crit*ones(size(r)),'r','LineWidth',2)
grid on 
xlabel('r')
ylabel('\Delta t_{c,crit}')
title(sprintf('std = %.2f',std(delta_crit(delta_crit>0))))

%% Look if CN' at dynamic stall time is constant when using tc_crit as Tp

CNprime_ds = -ones(size(c));

for k=1:length(c)
ramp = ramps{k};
ramp.BeddoesLeishman(airfoil,t_crit,1,0,ramp.estimateTvl,'experimental')
% we are looking at a model value at a point defined experimentally! Does that really make sense?
CNprime_ds(k) = ramp.CNprime(ramp.i_CLonset);
end

%% Plot CN'_ds for all experiments

figure 
plot(r,CNprime_ds,'d','MarkerFaceColor','k')
hold on 
yline(CNcrit,'r--','LineWidth',2)
grid on 
xlabel('r')
ylabel('C''_{N,ds}')