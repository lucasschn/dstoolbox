% This script has two puroposes:
% 1) Correcting some wrongly determined peaks produced by former versions
% of findPeaks methods. 
% 2) Add error fields to a res structure produced by paramsweep.m. The
% error fields corresponds to error on peak location and height for CNk,
% CNf, CNv and CN_LB curves.
% Author : Lucas Schneeberger
% Date : 05.09.2020
close all
clear all
clc

run(fullfile('..','labbook.m'))

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(path2static);
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);

%% Load the path to the mat-file containing the results of the parameters sweep
try
    % try to find the path2res variable in paths.mat first
    load(fullfile('..','paths.mat'),'path2res')
catch % if path2res was not found in paths.mat, ...
    try % .. maybe setPaths.m has never been run and paths.mat does not exist
        run('setPaths.m')
        load(fullfile('..','paths.mat'),'path2res')
    catch % ... or maybe path2res is ill-defined
        open('setPaths.m')
        error('The path to the paramsweep results mat-file has not been correctly set. Please correct the path in setPaths.m first.') 
    end
end


filename = 'res25uniform4';
load(fullfile(path2res,filename),'res')

%% Create the ramp corresponding to the pitch rate

c = 71;

ramp = loadRamp(c,true);

ramp.setPitchRate(airfoil);
% Define stall (convectime must have been set)
ramp.findExpOnset();

% define the steady state CN for usage in the computation of the error on
% the vortex lift
CNsteady = mean(ramp.CN(end-100:end));

%% Add fields to the res structure

for k=1:length(res)
    %% Wrong values 
    % when no vortex lift
    if res(k).Tv == 0 || res(k).Tvl  == 0 % means there is no vortex lift
        res(k).SmaxCNv = NaN;
        res(k).maxCNv = NaN;
        warning('Peak location has been replaced by NaN.')
    end
    % when no overshoot of steady-state lift
    if res(k).SmaxCNk > 30
        res(k).SmaxCNk = NaN;
        warning('Peak location has been replaced by NaN.')
    end
    %% Define primary peak time range    
    [ramp.maxCN,imaxCN] = max(ramp.CN);
    dCN = diff(ramp.CN(imaxCN:end)); % deltaCN after stall
    i_vortex_end = find(dCN>0.,1) + imaxCN;
    i_ramp_start = find(ramp.S>1e-4,1);
    
    %% Total lift errors
    % Only compute error until the end of primary peak
    res(k).err = mean((res(k).CN_LB(i_ramp_start:i_vortex_end) - ramp.CN(i_ramp_start:i_vortex_end)).^2);
    res(k).errPeakLoc = res(k).SmaxCN_LB - res(k).SmaxCN;
    res(k).errPeakHeight = res(k).maxCN_LB - res(k).maxCN;
    
    %% Kirchhof lift w/o added mass errors
    res(k).errCNk_PeakLoc = res(k).SmaxCNk - res(k).SmaxCN;
    res(k).errCNk_PeakHeight = res(k).maxCNk - res(k).maxCN;
    
    %% Kirchhoff lift w/ added mass errors
    res(k).errCNf_PeakLoc = res(k).SmaxCNf - res(k).SmaxCN;
    res(k).errCNf_PeakHeight = res(k).maxCNf - res(k).maxCN;
    
    %% Vortex lift error
    % Define reference as difference between maxCN and final value of CN
    refCNv = res(k).maxCN - CNsteady; 
    res(k).errCNv_PeakLoc = res(k).SmaxCNv - res(k).SmaxCN; % still compare to primary peak
    res(k).errCNv_PeakHeight = res(k).maxCNv - refCNv;
    
    %% Errors for peaks found with findpeaks
    i_start = find(res(k).CN_LB>=0,1);
    [peaks,peak_indices] = findpeaks(res(k).CN_LB(i_start:end),'MinPeakDistance',150);
    res(k).firstPeak = peaks(1);
    res(k).firstPeakLoc = ramp.S(peak_indices(1));
    if length(peaks)>1
        res(k).secondPeak = peaks(2);
        res(k).secondPeakLoc = ramp.S(peak_indices(2));
    end
    res(k).errFirstPeakLoc = res(k).firstPeakLoc - res(k).SmaxCN;
    res(k).errFirstPeakHeight = res(k).firstPeak - res(k).maxCN;
    if ~isempty(res(k).secondPeak)
        res(k).errSecondPeakLoc = res(k).secondPeakLoc - res(k).SmaxCN;
        res(k).errSecondPeakHeight = res(k).secondPeak - res(k).maxCN;
    end
    %% Add a 0-1 variable that states if a second peak has been detected or not
    res(k).hasSecondPeak = isempty(res(k).secondPeak);
end

save(fullfile(path2res,filename),'res','-v7.3')