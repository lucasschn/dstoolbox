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

c = [18,14,22,67,26,84,30,2,34,71,38,6,42,75,46,10,50,62];

% initialization of the vectors
r = -ones(size(c));
tc_ds = -ones(size(c));
tc_ss = -ones(size(c));
CN_ds = -ones(size(c));

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
    % ... and assign it
    tc_ss(k) = interp1(ramp.analpha,ramp.S,airfoil.steady.alpha_ss);   
    if ramp.alpha_CLonset > airfoil.steady.alpha_ss        
        tc_ds(k) = ramp.S(ramp.i_CLonset);
        CN_ds(k) = ramp.CN(ramp.i_CLonset);        
    else
        fit_error('Both CC- and CL-based stall angles are smaller than static stall angle.')
    end
end

figure
plot(r,tc_ds-tc_ss,'d','MarkerFaceColor','k')
xlabel('r')
ylabel('t_{c,ds}-t_{c,ss}')
grid on

pr = polyfit(r(r>airfoil.r0),CN_ds(r>airfoil.r0),0);
pl = lsqcurvefit(@(x,xdata) x*(xdata-airfoil.r0)+pr,20,r(r<airfoil.r0),CN_ds(r<airfoil.r0));

figure
plot(r,CN_ds,'d','MarkerFaceColor','k')
hold on
plot(sort([0 r airfoil.r0]),[pl*([0 r(r<airfoil.r0)]-airfoil.r0)+pr pr*ones(1,length(r(r>airfoil.r0))+1)],'r','LineWidth',1)
xlabel('r')
ylabel('C_{N,ds}')
grid on
title('C_{N,ds} = max(C_N)')
saveas(gcf,'../fig/Cds0','png')
CNcrit = pl*(0-airfoil.r0)+pr;
save('../CNcrit','CNcrit')
fprintf('The limit of CNds as r->0 is %1.3f. \n',CNcrit)