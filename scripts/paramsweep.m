% This script is made to run LB model on a broad range of parameters on
% SH2019 data
% Author : Lucas Schneeberger
% Date : 28.08.2020

close all
clear all
clc   
set(0,'DefaultFigureWindowStyle','docked')
addpath('../plot_dir/')
addpath('../src/model/')
addpath('../src/common/')
addpath('../src/lib/')
if ismac  
    run('/Users/lucas/src/codes_smarth/labbook.m')
elseif ispc
    run('labbook')
end
%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

%% Run the parameter sweep
params = -ones([1,4]);

c = 71; % list of tested pitch rates
n = 1e4; % number of samples per pitch rate

res = struct('comment',{});

for k = 1:length(c) % loop over the pitch rates
    data = load(loadmat(LB(c(k)).ms,LB(c(k)).mpt),'raw','inert','avg','zero');
    raw = data.raw;
    zero = data.zero;
    inert.alpha = raw.alpha(raw.t>=0);
    msname = sprintf('ms%03impt%i',LB(c(k)).ms,LB(c(k)).mpt);
    ramp = RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c(k)).U,'alphadot',LB(c(k)).alphadot);
    evalin('base',sprintf('ramp.setName(''%s'') ',msname))

    Cl = raw.Cl-zero.Cl;
    Cd = raw.Cd-zero.Cd;
    fs = 1/ramp.Ts;
    Clf = myFilterTwice(Cl,fs);
    Cdf = myFilterTwice(Cd,fs);
    ramp.setCL(Clf);
    ramp.setCD(Cdf);
    
    ramp.computeAirfoilFrame();
    ramp.isolateRamp();
    ramp.setPitchRate(airfoil);
    % Define stall (convectime must have been set)
    ramp.findExpOnset();     
    
    for kk = 1:n % do n samples per pitch rate
        params(1) = randsample(0:0.1:20,1);
        params(2) = randsample(0:0.1:10,1);
        params(3:4) = randsample(0:0.1:5,2,'true');
        disp(params)
        res(n*(k-1)+kk).r = ramp.r;
        res(n*(k-1)+kk).alphadot = ramp.alphadot;
        res(n*(k-1)+kk).Tp = params(1);
        res(n*(k-1)+kk).Tf = params(2);
        res(n*(k-1)+kk).Tv = params(3);
        res(n*(k-1)+kk).Tvl = params(4);
        
        ramp.BeddoesLeishman(airfoil,params(1),params(2),params(3),params(4),'experimental');
        
        % Postprocessing
        ramp.findPeaks()
        ramp.computeErr()
        
        % Put results in the struct        
        res(n*(k-1)+kk).CN_LB = ramp.CN_LB;
        res(n*(k-1)+kk).maxCN = ramp.maxCN;
        res(n*(k-1)+kk).maxCN_LB = ramp.maxCN_LB;
        res(n*(k-1)+kk).maxCNk = ramp.maxCNk;
        res(n*(k-1)+kk).maxCNf = ramp.maxCNf;
        res(n*(k-1)+kk).maxCNv = ramp.maxCNv;
        res(n*(k-1)+kk).SmaxCN = ramp.SmaxCN;
        res(n*(k-1)+kk).SmaxCN_LB = ramp.SmaxCN_LB;
        res(n*(k-1)+kk).SmaxCNk = ramp.SmaxCNk;
        res(n*(k-1)+kk).SmaxCNf = ramp.SmaxCNf;
        res(n*(k-1)+kk).SmaxCNv = ramp.SmaxCNv;
        res(n*(k-1)+kk).err = ramp.err;
    end
end

if ~isfolder(fullfile('..','data','paramsweep'))
    mkdir(fullfile('..','data','paramsweep'));
end
name = 'res25';
save(fullfile('..','data','paramsweep',name),'res');