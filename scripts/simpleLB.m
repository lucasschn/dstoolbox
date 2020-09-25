run labbook.m

addpath(fullfile('..','src','model'))
addpath(fullfile('..','src','common'))
addpath(fullfile('..','src','lib'))

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
static = load(fullfile('..','data','static_flatplate.mat'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

% define the case number of the desired experiment (according to labbook)
c = 71;

load(loadmat(LB(c).ms,LB(c).mpt),'raw','zero');

msname = sprintf('ms%03impt%i',LB(c).ms,LB(c).mpt);
ramp = RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c).U,'alphadot',LB(c).alphadot);
ramp.setName()

Cl = raw.Cl-zero.Cl;
Cd = raw.Cd-zero.Cd;
Cl = Cl - mean(Cl(1:50));
Cd = Cd - mean(Cd(1:50));
fs = 1/ramp.Ts;

% filter the load data
Cl_fff = myFilterTwice(Cl,fs);
Cd_fff = myFilterTwice(Cd,fs);
ramp.setCL(Cl_fff);
ramp.setCD(Cd_fff);

ramp.computeAirfoilFrame();
ramp.isolateRamp();

% ramp.BeddoesLeishman(airfoil,Tp,Tf,Tv,Tvl,'experimental')
ramp.BeddoesLeishman(airfoil,1,1,1,1,'experimental')
ramp.plotCustom('CN','CN_LB')