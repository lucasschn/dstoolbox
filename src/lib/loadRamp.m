function ramp = loadRamp(c,filtered)
% loads ramp experiment number c. filtered argument is a boolean that specifies if the
% signals should be filtered (default) or not.
disp(pwd)

try
    run(fullfile('..','labbook.m'))
catch
    error('Labbook was not found at any paths.')
end

if length(c) > 1
    error('loadRamp can only manage one experiment at a time. Verify your case number size.')
end

try
    load(loadmat(LB(c).ms,LB(c).mpt),'raw','zero');
catch IOError    
    error('Matlab could not read the experimental data. Are you sure you are connected to the raw server?')
end

msname = sprintf('ms%03impt%i',LB(c).ms,LB(c).mpt);
assignin('base',msname,RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c).U,'alphadot',LB(c).alphadot));
evalin('base',sprintf('%s.setName()',msname))
ramp = evalin('base',msname);
Cl = raw.Cl-zero.Cl;
Cd = raw.Cd-zero.Cd;
Cl = Cl - mean(Cl(1:50));
Cd = Cd - mean(Cd(1:50));
fs = 1/ramp.Ts;
if nargin > 1 && filtered
    Cl_fff = myFilterTwice(Cl,fs);
    Cd_fff = myFilterTwice(Cd,fs);
    ramp.setCL(Cl_fff);
    ramp.setCD(Cd_fff);
else
    ramp.setCL(Cl);
    ramp.setCD(Cd);
end
ramp.computeAirfoilFrame();
ramp.isolateRamp();
end