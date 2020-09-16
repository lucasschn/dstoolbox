close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run(fullfile('..','labbook.m'))

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,8);

%% Define ramps and plot tc_ds(r)

c = [18,14,22,67,26,84,30,2,34,71,38,6,42,75,46,10,50,62];

% initialization of the vectors
r = -ones(size(c));
tc_ds = -ones(size(c));
tc_ss = -ones(size(c));
CN_ds = -ones(size(c));

for k=1:length(c)    
    ramp = loadRamp(c(k));
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
    i_start = find(ramp.analpha == 0,1,'last');
    i_end = find(ramp.analpha >= 30,1);
    tc_ss(k) = interp1(ramp.analpha(i_start:i_end),ramp.S(i_start:i_end),airfoil.steady.alpha_ss);   
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