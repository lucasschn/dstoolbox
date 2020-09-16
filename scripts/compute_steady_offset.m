% This script is to measure the steady state offset for each pitch rate of
% SH2019 data
% Author : Lucas Schneeberger
% Date : 05.06.2020

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

%% Setting up the ramps

c = [18,14,22,67,26,84,30,2,34,71,38,6,42,75,46,10,50,62];

%% Run Beddoes Leishman model on all ramps

r = -ones(size(c));
steady_exp = -ones(size(c));
steady_LB = -ones(size(c));
steady_offset = -ones(size(c));

for k=1:length(c) 
    ramp = loadRamp(c(k));
    ramp.setPitchRate(airfoil);
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    ramp.BeddoesLeishman(airfoil,4.5,4,6,1,'experimental')
    saveas(gcf,fullfile('..','fig',sprintf('LB_r%03d',round(ramp.r,3)*1000)),'png')
    r(k) = ramp.r;
    steady_exp(k) = evalin('base',sprintf('mean(%s.CN(%s.S > 0.9*%s.S(end)))',msname,msname,msname));
    steady_LB(k) = evalin('base',sprintf('%s.CN_LB(end)',msname));
end

steady_offset = steady_exp - steady_LB;

%Compute the static lift offset between experiment and Kirchhoff
kirchhoff_offset = airfoil.steady.CN(end)-kirchhoff(airfoil.steady,30);

% Plot and compare results
figure
plot(r,steady_offset,'.','MarkerSize',20,'Display','steady state offset')
%hold on 
%plot(r,kirchhoff_offset*ones(size(r)),'r-.','DisplayName','static offset due to Kirchhoff')
ax = gca; 
ax.FontSize = 20; 
xlabel('r')
ylabel('\Delta C_{N,\infty}')
grid on 


