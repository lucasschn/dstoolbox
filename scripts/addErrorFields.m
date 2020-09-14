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
if ismac
    run('/Users/lucas/src/codes_smarth/labbook.m')
    path2fig = '../fig';
    path2static = fullfile('..','static_flatplate');
elseif ispc 
    run('labbook')
    path2fig = '\\oscar\macintosh hd\Users\lucas\Documents\EPFL\PDM\fig';
    path2static = '\\oscar\macintosh hd\Users\lucas\Documents\EPFL\PDM\static_flatplate';
end


%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(path2static);
airfoil.steady = SteadyCurve(static.alpha,static.CN,13);


filename = 'res25uniform3';
load(fullfile('..','data','paramsweep',filename))
 
c = 71;

data = load(loadmat(LB(c).ms,LB(c).mpt),'raw','inert','avg','zero');
raw = data.raw;
zero = data.zero;
msname = sprintf('ms%03impt%i',LB(c).ms,LB(c).mpt);
ramp = RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c).U,'alphadot',LB(c).alphadot);
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

% define the steady state CN for usage in the computation of the error on
% the vortex lift
CNsteady = mean(ramp.CN(end-100:end));

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
    res(k).err = (res(k).CN_LB(i_ramp_start:i_vortex_end) - ramp.CN(i_ramp_start:i_vortex_end)).^2;
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
    
end

save(fullfile('..','data','paramsweep',filename),'res')