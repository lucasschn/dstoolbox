%%%
% add short summary of measurements
% when: 15.02.19 + 19.02.19 ...
% who: Sabrina
% where: sharx
% what: shorter version of lb3, only one RE,different alphadot
% 2D2C high speed PIV for some cases; 2 cameras + mirrors; forces for
% all cases; ramp motion flat plate with two sharp edges
% project: smartH
%%%

disp('Labbook has been found')

try
    load('paths','path2fig')
catch
    try 
        load(fullfile('..','paths'),'path2fig')
    catch
        open(fullfile('scripts','setPaths.m'))
        error('The path to the figure folder has not been set. Please set your paths in setPaths.m and run this script again.')
    end
end

path2static = fullfile('..','static_flatplate');
if ismac
    path2smarth = '/Volumes/unfold/smartH/2019_SH';
elseif ispc
    path2smarth = '\\sti1raw.epfl.ch\unfold\smartH\2019_SH';
else
    error('Your platform is not supported.')   
end

%%% where is data - where will data go?
    root.raw = @(nr3) fullfile(path2smarth,sprintf('%i',nr3),'loads'); % Raw data location
    root.data=@(nr3) fullfile(path2smarth,sprintf('%i',nr3));
    root.res = fullfile(path2smarth,'results');
    root.pivmat = fullfile(path2smarth,'postprocessing','matfiles','piv');
    root.loadmat = fullfile(path2smarth,'postprocessing','matfiles','loads');%'\\sti1raw.epfl.ch\unfold\smartH\2019_SH\Postprocessing\matfiles\loads';
    root.fig = fullfile(path2smarth,'postprocessing','matfiles','figurematter');
    root.matlab = fullfile(path2smarth,'matlab');
    root.ftle = fullfile(path2smarth,'postprocessing','ftle');
    root.pod = fullfile(path2smarth,'postprocessing','pod');
    root.pp = fullfile(path2smarth,'postprocessing');

% root.raw = '\\sti1raw.epfl.ch\unfold\smartH\2018_SH_004\20180629\loads'; % Raw data location

%%% shortcuts to frequently used files
raw_force = @(nr3,nr1,nr2) fullfile(root.raw(nr3),sprintf('ms%.3impt%.3i_FORCES*.bin',nr1,nr2)); % Force data
raw_pos = @(nr3,nr1,nr2) fullfile(root.raw(nr3),sprintf('ms%.3impt%.3i_POSITION*.bin',nr1,nr2)); % Position data
raw_par = @(nr3,nr1,nr2) fullfile(root.raw(nr3),sprintf('ms%.3impt%.3i*.txt',nr1,nr2)); % Parameters

pivmat_avg = @(nr1) fullfile(root.pivmat,sprintf('ms%.3i_avg.mat',nr1));
pivmat = @(nr1,nr2) fullfile(root.pivmat,sprintf('ms%.3impt%.3i',nr1,nr2),sprintf('ms%.3impt%.3i.mat',nr1,nr2));
% pivmat_smooth = @(nr1,nr2) fullfile(root.pivmat,sprintf('ms%.3impt%.3i_smooth.mat',nr1,nr2));
loadmat = @(nr1,nr2) fullfile(root.loadmat, sprintf('ms%.3impt%.3i_loads.mat',nr1,nr2));
ftlemat = @(nrpod,nr1,nr2) fullfile(root.ftle,sprintf('ftle_pod%.2i_ms%.3impt%.3i.mat',nrpod,nr1,nr2));
podmatrec = @(nr1,nr2,modes) fullfile(root.pod,sprintf('ms%.3impt%.3i_pod%.2i.mat',nr1,nr2,modes));
podmat = @(nr1,nr2) fullfile(root.pod,sprintf('pod_ms%.3impt%.3i.mat',nr1,nr2));
circmat = @(nr1,nr2) fullfile(root.pivmat,'circulation',sprintf('circ_ms%.3impt%.3i.mat',nr1,nr2));
% vortexmat = fullfile(root.pivmat,tracking);

pressuremat = @(thecase) fullfile(root.mat,sprintf('ms%.3impt%.3i_p.mat',nr1,nr2));
fig = @(thecase,what) fullfile(root.fig,sprintf('ms%.3impt%.3i_%s',thecase,what));
instfig = @(nr1,nr2,what,n) fullfile(root.fig,'plots',sprintf('ms%.3impt%.3i_%s_%.5i.png',nr1,nr2,what,n));
resdat = @(nr1,what) fullfile('C:\Users\sabrinahenne\Documents\SmartH\2019_SH\writing\figures\test\datfiles',sprintf('ms%.3i_%s.dat',nr1,what));
plotmat = @(nr1,what) fullfile('C:\Users\sabrinahenne\Documents\SmartH\2019_SH\postprocessing\matfiles\plotting',sprintf('ms%.3i_%s.mat',nr1,what));
% resdat = @(thecase,what) fullfile(root.fig,'datfiles',sprintf('ms%.3impt%.3i_%s.dat',LB(thecase).ms,LB(thecase).mpt,what));

%%% parameters that are the same for all tested cases
param.c=0.15;       % chorch length in m
param.cmm=15;       % chord length in mm
param.s=0.6;        % span in m
param.c0=0.25;      % position of pitching axis in /c
param.nrowpx=1024;  % number of row of raw image in px
param.ncolpx=1024;  % number of col of raw image in px
param.mu=0.89e-3;   % dynamic viscosity in Pas
param.rho=1000;      % density of water in kg/m^3
%%% some of them can be determined after the measurements

%%% List of MS that should not be saved w/ correlated data. Use for calibration,
%%% static angle sweeps
param.nocorr = [2 24 33 115];

%%% Subdivision into groups of MS that are associated by calibration & Re
grp.ms = { 7:23 25:28 34 116:120};
grp.mpt = {1:4 1:4 1:4 1:4};
grp.mscal = [2 24 33 115];

%%% specific parameters for all testcases
LB = struct('comment' , {});

%%% Set group properties quickly
for i = 1:108
    LB(i).U = 0.5;
    LB(i).Re = param.rho*LB(i).U*param.c/param.mu;
end


%% ramp up motion 0° - 30°, 20s of data

currentindex=1;
LB(currentindex).ms=2;
LB(currentindex).mpt=1;
LB(currentindex).date=20190215;
LB(currentindex).alpha=0; 
LB(currentindex).alphadot=0; %static
LB(currentindex).lcf=1000;
LB(currentindex).comment=' calibration';

currentindex=2; 
LB(currentindex).ms=7;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=3; 
LB(currentindex).ms=7;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=4; 
LB(currentindex).ms=7;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=5; 
LB(currentindex).ms=7;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=6; 
LB(currentindex).ms=8;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=7; 
LB(currentindex).ms=8;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=8; 
LB(currentindex).ms=8;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=9; 
LB(currentindex).ms=8;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=10; 
LB(currentindex).ms=9;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=11; 
LB(currentindex).ms=9;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=12; 
LB(currentindex).ms=9;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=13; 
LB(currentindex).ms=9;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=14; 
LB(currentindex).ms=10;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=5.00;
LB(currentindex).alphatdotdeg=0.75; 
LB(currentindex).alphatdotrad=0.013; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=15; 
LB(currentindex).ms=10;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=5.00;
LB(currentindex).alphatdotdeg=0.75; 
LB(currentindex).alphatdotrad=0.013; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=16; 
LB(currentindex).ms=10;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=5.00;
LB(currentindex).alphatdotdeg=0.75; 
LB(currentindex).alphatdotrad=0.013; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

currentindex=17; 
LB(currentindex).ms=10;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190215; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=5.00;
LB(currentindex).alphatdotdeg=0.75; 
LB(currentindex).alphatdotrad=0.013; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='first test';

%%

currentindex=18; 
LB(currentindex).ms=12;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=2.50;
LB(currentindex).alphatdotdeg=0.38; 
LB(currentindex).alphatdotrad=0.007; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=19; 
LB(currentindex).ms=12;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=2.50;
LB(currentindex).alphatdotdeg=0.38; 
LB(currentindex).alphatdotrad=0.007; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=20; 
LB(currentindex).ms=12;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=2.50;
LB(currentindex).alphatdotdeg=0.38; 
LB(currentindex).alphatdotrad=0.007; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=21; 
LB(currentindex).ms=12;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=2.50;
LB(currentindex).alphatdotdeg=0.38; 
LB(currentindex).alphatdotrad=0.007; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=22; 
LB(currentindex).ms=13;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=7.50;
LB(currentindex).alphatdotdeg=1.13; 
LB(currentindex).alphatdotrad=0.020; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=23; 
LB(currentindex).ms=13;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=7.50;
LB(currentindex).alphatdotdeg=1.13; 
LB(currentindex).alphatdotrad=0.020; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=24; 
LB(currentindex).ms=13;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=7.50;
LB(currentindex).alphatdotdeg=1.13; 
LB(currentindex).alphatdotrad=0.020; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=25; 
LB(currentindex).ms=13;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=7.50;
LB(currentindex).alphatdotdeg=1.13; 
LB(currentindex).alphatdotrad=0.020; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=26; 
LB(currentindex).ms=14;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=12.50;
LB(currentindex).alphatdotdeg=1.88; 
LB(currentindex).alphatdotrad=0.033; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=27; 
LB(currentindex).ms=14;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=12.50;
LB(currentindex).alphatdotdeg=1.88; 
LB(currentindex).alphatdotrad=0.033; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=28; 
LB(currentindex).ms=14;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=12.50;
LB(currentindex).alphatdotdeg=1.88; 
LB(currentindex).alphatdotrad=0.033; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=29; 
LB(currentindex).ms=14;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=12.50;
LB(currentindex).alphatdotdeg=1.88; 
LB(currentindex).alphatdotrad=0.033; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=30; 
LB(currentindex).ms=15;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=31; 
LB(currentindex).ms=15;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=32; 
LB(currentindex).ms=15;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=33; 
LB(currentindex).ms=15;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=34; 
LB(currentindex).ms=16;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=35; 
LB(currentindex).ms=16;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=36; 
LB(currentindex).ms=16;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=37; 
LB(currentindex).ms=16;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=38; 
LB(currentindex).ms=17;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=27.50;
LB(currentindex).alphatdotdeg=4.13; 
LB(currentindex).alphatdotrad=0.072; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=39; 
LB(currentindex).ms=17;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=27.50;
LB(currentindex).alphatdotdeg=4.13; 
LB(currentindex).alphatdotrad=0.072; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=40; 
LB(currentindex).ms=17;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=27.50;
LB(currentindex).alphatdotdeg=4.13; 
LB(currentindex).alphatdotrad=0.072; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=41; 
LB(currentindex).ms=17;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=27.50;
LB(currentindex).alphatdotdeg=4.13; 
LB(currentindex).alphatdotrad=0.072; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=42; 
LB(currentindex).ms=18;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=32.50;
LB(currentindex).alphatdotdeg=4.88; 
LB(currentindex).alphatdotrad=0.085; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=43; 
LB(currentindex).ms=18;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=32.50;
LB(currentindex).alphatdotdeg=4.88; 
LB(currentindex).alphatdotrad=0.085; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=44; 
LB(currentindex).ms=18;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=32.50;
LB(currentindex).alphatdotdeg=4.88; 
LB(currentindex).alphatdotrad=0.085; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=45; 
LB(currentindex).ms=18;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=32.50;
LB(currentindex).alphatdotdeg=4.88; 
LB(currentindex).alphatdotrad=0.085; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=46; 
LB(currentindex).ms=19;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=37.50;
LB(currentindex).alphatdotdeg=5.63; 
LB(currentindex).alphatdotrad=0.098; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=47; 
LB(currentindex).ms=19;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=37.50;
LB(currentindex).alphatdotdeg=5.63; 
LB(currentindex).alphatdotrad=0.098; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=48; 
LB(currentindex).ms=19;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=37.50;
LB(currentindex).alphatdotdeg=5.63; 
LB(currentindex).alphatdotrad=0.098; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=49; 
LB(currentindex).ms=19;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=37.50;
LB(currentindex).alphatdotdeg=5.63; 
LB(currentindex).alphatdotrad=0.098; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=50; 
LB(currentindex).ms=20;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=42.50;
LB(currentindex).alphatdotdeg=6.38; 
LB(currentindex).alphatdotrad=0.111; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=51; 
LB(currentindex).ms=20;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=42.50;
LB(currentindex).alphatdotdeg=6.38; 
LB(currentindex).alphatdotrad=0.111; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=52; 
LB(currentindex).ms=20;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=42.50;
LB(currentindex).alphatdotdeg=6.38; 
LB(currentindex).alphatdotrad=0.111; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=53; 
LB(currentindex).ms=20;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=42.50;
LB(currentindex).alphatdotdeg=6.38; 
LB(currentindex).alphatdotrad=0.111; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=54; 
LB(currentindex).ms=21;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=55; 
LB(currentindex).ms=21;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=56; 
LB(currentindex).ms=21;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=57; 
LB(currentindex).ms=21;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=58; 
LB(currentindex).ms=22;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=47.50;
LB(currentindex).alphatdotdeg=7.13; 
LB(currentindex).alphatdotrad=0.124; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=59; 
LB(currentindex).ms=22;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=47.50;
LB(currentindex).alphatdotdeg=7.13; 
LB(currentindex).alphatdotrad=0.124; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=60; 
LB(currentindex).ms=22;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=47.50;
LB(currentindex).alphatdotdeg=7.13; 
LB(currentindex).alphatdotrad=0.124; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=61; 
LB(currentindex).ms=22;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=47.50;
LB(currentindex).alphatdotdeg=7.13; 
LB(currentindex).alphatdotrad=0.124; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=62; 
LB(currentindex).ms=23;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=50.00;
LB(currentindex).alphatdotdeg=7.50; 
LB(currentindex).alphatdotrad=0.131; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=63; 
LB(currentindex).ms=23;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=50.00;
LB(currentindex).alphatdotdeg=7.50; 
LB(currentindex).alphatdotrad=0.131; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=64; 
LB(currentindex).ms=23;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=50.00;
LB(currentindex).alphatdotdeg=7.50; 
LB(currentindex).alphatdotrad=0.131; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

currentindex=65; 
LB(currentindex).ms=23;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190220; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=50.00;
LB(currentindex).alphatdotdeg=7.50; 
LB(currentindex).alphatdotrad=0.131; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up';

%%

currentindex=66;
LB(currentindex).ms=24;
LB(currentindex).mpt=1;
LB(currentindex).date=20190228;
LB(currentindex).alpha=0; 
LB(currentindex).alphadot=0;
LB(currentindex).lcf=1000;
LB(currentindex).comment=' calibration';

currentindex=67; 
LB(currentindex).ms=25;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=10.00;
LB(currentindex).alphatdotdeg=1.50; 
LB(currentindex).alphatdotrad=0.026; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=68; 
LB(currentindex).ms=25;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=10.00;
LB(currentindex).alphatdotdeg=1.50; 
LB(currentindex).alphatdotrad=0.026; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=69; 
LB(currentindex).ms=25;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=10.00;
LB(currentindex).alphatdotdeg=1.50; 
LB(currentindex).alphatdotrad=0.026; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=70; 
LB(currentindex).ms=25;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=10.00;
LB(currentindex).alphatdotdeg=1.50; 
LB(currentindex).alphatdotrad=0.026; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=71; 
LB(currentindex).ms=26;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=25.00;
LB(currentindex).alphatdotdeg=3.75; 
LB(currentindex).alphatdotrad=0.065; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=72; 
LB(currentindex).ms=26;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=25.00;
LB(currentindex).alphatdotdeg=3.75; 
LB(currentindex).alphatdotrad=0.065; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=73; 
LB(currentindex).ms=26;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=25.00;
LB(currentindex).alphatdotdeg=3.75; 
LB(currentindex).alphatdotrad=0.065; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=74; 
LB(currentindex).ms=26;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=25.00;
LB(currentindex).alphatdotdeg=3.75; 
LB(currentindex).alphatdotrad=0.065; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=75; 
LB(currentindex).ms=27;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=35.00;
LB(currentindex).alphatdotdeg=5.25; 
LB(currentindex).alphatdotrad=0.092; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=76; 
LB(currentindex).ms=27;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=35.00;
LB(currentindex).alphatdotdeg=5.25; 
LB(currentindex).alphatdotrad=0.092; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=77; 
LB(currentindex).ms=27;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=35.00;
LB(currentindex).alphatdotdeg=5.25; 
LB(currentindex).alphatdotrad=0.092; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=78; 
LB(currentindex).ms=27;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=35.00;
LB(currentindex).alphatdotdeg=5.25; 
LB(currentindex).alphatdotrad=0.092; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=79; 
LB(currentindex).ms=28;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=80; 
LB(currentindex).ms=28;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=81; 
LB(currentindex).ms=28;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=82; 
LB(currentindex).ms=28;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190228; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=45.00;
LB(currentindex).alphatdotdeg=6.75; 
LB(currentindex).alphatdotrad=0.118; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

%% originally 50 cases, here 4

currentindex=83;
LB(currentindex).ms=33;
LB(currentindex).mpt=1;
LB(currentindex).date=20190311;
LB(currentindex).alpha=0; 
LB(currentindex).alphadot=0;
LB(currentindex).lcf=1000;
LB(currentindex).comment=' calibration';

currentindex=84; 
LB(currentindex).ms=34;
LB(currentindex).mpt=1; 
LB(currentindex).date=20190311; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=15.00;
LB(currentindex).alphatdotdeg=2.25; 
LB(currentindex).alphatdotrad=0.039; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=85; 
LB(currentindex).ms=34;
LB(currentindex).mpt=2; 
LB(currentindex).date=20190311; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=15.00;
LB(currentindex).alphatdotdeg=2.25; 
LB(currentindex).alphatdotrad=0.039; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=86; 
LB(currentindex).ms=34;
LB(currentindex).mpt=3; 
LB(currentindex).date=20190311; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=15.00;
LB(currentindex).alphatdotdeg=2.25; 
LB(currentindex).alphatdotrad=0.039; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=87; 
LB(currentindex).ms=34;
LB(currentindex).mpt=4; 
LB(currentindex).date=20190311; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=15.00;
LB(currentindex).alphatdotdeg=2.25; 
LB(currentindex).alphatdotrad=0.039; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

%%

currentindex=88;
LB(currentindex).ms=115;
LB(currentindex).mpt=1;
LB(currentindex).date=20200204;
LB(currentindex).alpha=0; 
LB(currentindex).alphadot=0;
LB(currentindex).lcf=1000;
LB(currentindex).comment=' calibration';

currentindex=89; 
LB(currentindex).ms=116;
LB(currentindex).mpt=1; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=90; 
LB(currentindex).ms=116;
LB(currentindex).mpt=2; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=91; 
LB(currentindex).ms=116;
LB(currentindex).mpt=3; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=92; 
LB(currentindex).ms=116;
LB(currentindex).mpt=4; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=20.00;
LB(currentindex).alphatdotdeg=3.00; 
LB(currentindex).alphatdotrad=0.052; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=93; 
LB(currentindex).ms=117;
LB(currentindex).mpt=1; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=94; 
LB(currentindex).ms=117;
LB(currentindex).mpt=2; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=95; 
LB(currentindex).ms=117;
LB(currentindex).mpt=3; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=96; 
LB(currentindex).ms=117;
LB(currentindex).mpt=4; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=30.00;
LB(currentindex).alphatdotdeg=4.50; 
LB(currentindex).alphatdotrad=0.079; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=97; 
LB(currentindex).ms=118;
LB(currentindex).mpt=1; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=98; 
LB(currentindex).ms=118;
LB(currentindex).mpt=2; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=99; 
LB(currentindex).ms=118;
LB(currentindex).mpt=3; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=100; 
LB(currentindex).ms=118;
LB(currentindex).mpt=4; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=40.00;
LB(currentindex).alphatdotdeg=6.00; 
LB(currentindex).alphatdotrad=0.105; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=101; 
LB(currentindex).ms=119;
LB(currentindex).mpt=1; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=102; 
LB(currentindex).ms=119;
LB(currentindex).mpt=2; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=103; 
LB(currentindex).ms=119;
LB(currentindex).mpt=3; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=104; 
LB(currentindex).ms=119;
LB(currentindex).mpt=4; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=17.50;
LB(currentindex).alphatdotdeg=2.63; 
LB(currentindex).alphatdotrad=0.046; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=105; 
LB(currentindex).ms=120;
LB(currentindex).mpt=1; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=106; 
LB(currentindex).ms=120;
LB(currentindex).mpt=2; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=107; 
LB(currentindex).ms=120;
LB(currentindex).mpt=3; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

currentindex=108; 
LB(currentindex).ms=120;
LB(currentindex).mpt=4; 
LB(currentindex).date=20200204; 
LB(currentindex).alpha=30; 
LB(currentindex).alphadot=22.50;
LB(currentindex).alphatdotdeg=3.38; 
LB(currentindex).alphatdotrad=0.059; 
LB(currentindex).PIV=1; 
LB(currentindex).f_acqu=500; 
LB(currentindex).N=10000; 
LB(currentindex).dt=2000; 
LB(currentindex).lcf=1000;
LB(currentindex).comment='ramp-up + hsPIV';

