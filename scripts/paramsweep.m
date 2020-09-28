% This script is made to run LB model on a broad range of parameters on
% SH2019 data
% Author : Lucas Schneeberger
% Date : 28.08.2020

close all
clear all
clc   
set(0,'DefaultFigureWindowStyle','docked')

run(fullfile('..','labbook'))

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','data','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

%% Run the parameter sweep
params = -ones([1,4]);

c = 71; % list of tested pitch rates
n = 1e4; % number of samples per pitch rate

res = struct('comment',{});

for k = 1:length(c) % loop over the pitch rates    
    ramp = loadRamp(c(k),true);
    ramp.setPitchRate(airfoil);
    % Define stall (convectime must have been set)
    ramp.findExpOnset();     
    
    for kk = 1:n % do n samples per pitch rate
        params(1) = randsample(0:0.1:20,1);
        params(2) = randsample(0:0.1:10,1);
        % does it really make sense to make Tv and Tvl start at 0? Vortex
        % lift will be zero
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
        ramp.computeErrors()
        
        % Put results in the struct
        res(n*(k-1)+kk).CN_LB = ramp.CN_LB;
        res(n*(k-1)+kk).CNk = ramp.CNk;
        res(n*(k-1)+kk).CNf = ramp.CNf;
        res(n*(k-1)+kk).CNv = ramp.CNv;
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
        res(n*(k-1)+kk).errCNk_PeakLoc = ramp.errCNk_PeakLoc;
        res(n*(k-1)+kk).errCNk_PeakHeight = ramp.errCNk_PeakHeight;
        res(n*(k-1)+kk).errCNf_PeakLoc = ramp.errCNf_PeakLoc;
        res(n*(k-1)+kk).errCNf_PeakHeight = ramp.errCNf_PeakHeight;
        res(n*(k-1)+kk).errCNv_PeakLoc = ramp.errCNv_PeakLoc;
        res(n*(k-1)+kk).errCNv_PeakHeight = ramp.errCNv_PeakHeight;
        res(n*(k-1)+kk).errPeakLoc = ramp.errPeakLoc;
        res(n*(k-1)+kk).errPeakHeight = ramp.errPeakHeight;
        res(n*(k-1)+kk).errFirstPeakLoc = ramp.errFirstPeakLoc;
        res(n*(k-1)+kk).errFirstPeakHeight = ramp.errFirstPeakHeight;
        res(n*(k-1)+kk).hasSecondPeak = isempty(ramp.secondPeak);
        res(n*(k-1)+kk).errSecondPeakLoc = ramp.errSecondPeakLoc;
        res(n*(k-1)+kk).errSecondPeakHeight = ramp.errSecondPeakHeight;
        
    end
end

if ~isfolder(fullfile('..','data','paramsweep'))
    mkdir(fullfile('..','data','paramsweep'));
end
name = 'res25uniform4';
save(fullfile('..','data','paramsweep',name),'res','-v7.3');