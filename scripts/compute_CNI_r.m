% This script computes the added mass value for each experiment
% Author : Lucas Schneeberger
% Date : 27.08.2020
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

%% Look at the time lag between t_crit and t_ds

c = [18,14,22,67,26,84,30,2,34,71,38,6,42,75,46,10,50,62];

% initialization of the variables and vectors
r = -ones(size(c));
added_mass = -ones(size(c));
ramps = cell(size(c));

for k = 1:length(c)   
    load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','zero');
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot));
    evalin('base',sprintf('%s.setName()',msname))
    ramp = evalin('base',msname);
    Cl = raw.Cl-zero.Cl;
    Cd = raw.Cd-zero.Cd;
    fs = 1/ramp.Ts;
    Cl_fff = myFilterTwice(Cl,fs);
    Cd_fff = myFilterTwice(Cd,fs);
    ramp.setCL(Cl_fff);
    ramp.setCD(Cd_fff);
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    ramps{k} = ramp;
    % Define reduced pitch rate if necessary ...
    if isempty(ramp.r)
        % compute it with alphadot
        ramp.setPitchRate(airfoil);
    end
    % ... and assign it
    r(k) = ramp.r;
    % Execute LB
    ramp.BeddoesLeishman(airfoil,5.48,1,1,1.7,'experimental')
    added_mass(k) = max(ramp.CNI);
end

%% Plot it as a function of the pitch rate

figure 
plot(r,added_mass,'d','MarkerFaceColor','k','DisplayName','max C_N^I')
grid on 
xlabel('r')
ylabel('C_N')
legend('Location','SouthEast')