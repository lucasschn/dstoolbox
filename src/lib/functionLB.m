function out = functionLB(input,LBcoeffs)

input=checkOptionNP(input);
LBcoeffs=checkParamNP(LBcoeffs);

% declares an airfoil with associated steady curve
airfoil = Airfoil('airfoil',input.c);
airfoil.steady = SteadyCurve(input.static.alpha,input.static.Cn,input.alpha_ss);

% declares a ramp-up experiment
ramp = RampUpMotion('alpha',input.dyn.alpha,'t',input.dyn.t,'V',input.U,'alphadot',input.alphadot);
ramp.setCL(input.dyn.Cl)
ramp.setCD(input.dyn.Cd)
ramp.computeAirfoilFrame();

% uncomment this if you don't want the whole signal but just the ramp part
%  ramp.isolateRamp();

ramp.BeddoesLeishman(airfoil,LBcoeffs.Tp,LBcoeffs.Tf,LBcoeffs.Tv,LBcoeffs.Tvl,'experimental')
out = struct(ramp);
end

function structout = checkOptionNP(structin)
% check and assign input options

% default values
default.c = 1; % chord length
default.static.alpha = linspace(0,20,41);
default.static.Cn = 2*pi*linspace(0,20,41);
default.alpha_ss = 12;
default.U=1;
default.alphadot=0.2; % in what unit?
default.dyn.t = linspace(0,20,101);
default.dyn.alpha = 0.2*linspace(0,20,101); 
default.dyn.Cl = 2*pi*0.2*linspace(0,20,101);
default.dyn.Cd = zeros(1,101);

champs=fields(default);

% set empty fields to default fields
% tochange=champs(~isfield(structin,champs));
tochange=champs(~cellfun(@(x) isfield(structin,x),champs));

for tochang=tochange'
    structin.(tochang{:}) = default.(tochang{:});
end

structout=structin;
end

function structout = checkParamNP(structin)
% check and assign param options

% default values
default.Tp = 1; % chord length
default.Tf = 1;
default.Tvl = 1;
default.Tv = 1;

champs=fields(default);

% set empty fields to default fields
% tochange=champs(~isfield(structin,champs));
tochange=champs(~cellfun(@(x) isfield(structin,x),champs));

for tochang=tochange'
    structin.(tochang{:}) = default.(tochang{:});
end

structout=structin;
end