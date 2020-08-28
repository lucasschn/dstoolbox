close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath(fullfile('..','src','model')))
run('labbook.m')

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load('../static_flatplate');
airfoil.steady = SteadyCurve(static.alpha,static.CN,8);

%% Define ramps and plot tc_ds(r)

c = [14,22,67,26,84,30,2,34,71,38,6,42,75,46,10,50,62];

% initialization of the vectors
r = -ones(size(c));
Tvl = -ones(size(c));

for k=1:length(c)    
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

    % Define reduced pitch rate if necessary ...
    if isempty(ramp.r)
        % compute it with alphadot
        ramp.setPitchRate(airfoil);
    end
    % ... and assign it
    r(k) = ramp.r;
    % Define experimental stall if necessary ...
    if isempty(ramp.i_CConset)
        ramp.findExpOnset();
    end
    % ... and compute the Tvl estimate according to Boutet 
    Tvl(k) = ramp.estimateTvl;
end

figure
plot(r,Tvl,'d','MarkerFaceColor','k')
xlabel('r')
ylabel('T_{vl}')
grid on