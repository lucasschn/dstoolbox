% now testing for simcos data

run labbook_simcos.m 
 
c=13; 

% for rampup motion
% load(loadmat(LB(c).ms,LB(c).mpt),'raw','zero')
static = load(fullfile('..','data','static_flatplate.mat'));

% for pitching motion
data = load(pressuredata(c));

% Drag missing, compute the pressure drag 
CD = zeros(size(data.Cl));
for k=1:length(data.Cl)
    CD(k) = sum(data.Cp(k,:)*data.xk);
end

load(fullfile('..','dynamic_corr'))

n=1000; 

input.c = param.c; 
input.static.alpha = static.alpha;
input.static.Cn = static.CN;
input.alpha_ss = 13;
% input.U=LB(c).U;
input.U = param.U0;
input.alphadot=0.2; % deg/s 
% input.dyn.t = raw.t;
TS = 1/LB(c).FS;
input.dyn.t = 0:TS:(n-1)*TS;
% input.dyn.alpha = raw.alpha; 
input.dyn.alpha = 10*sin(400*input.dyn.t);
%input.dyn.Cl = raw.Cl;
%input.dyn.Cd = raw.Cd;
input.dyn.Cl = rand(size(input.dyn.t));
input.dyn.Cd = zeros(size(input.dyn.t));
% sinus parameters
input.freq = LB(c).fosc;
input.mean_rad = deg2rad(LB(c).alpha_0);
input.amp_rad = deg2rad(LB(c).alpha_1);

% Leishman-Beddoes constants
LBcoeffs.Tp = 3;
LBcoeffs.Tf = 3;
LBcoeffs.Tv = 1;
LBcoeffs.Tvl = 1;

out = functionLB(input,LBcoeffs,'general');

figure, hold all
plot(input.dyn.t,input.dyn.Cl)
plot(input.dyn.t,out.CN_LB)
