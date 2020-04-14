close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('data/2008_simcos/matlab/labbook.m')

%% flow parameters
V = 50; %m/s
Re = 9.2e5; % chord as characteristic length
M = 0.14; 


%% Static data

% static_data = csvread('fatma_data.csv');
% static_data = csvread('T1_Re3.070_M0.00_N5.0_XtrTop 0%.csv',10);

airfoil = Airfoil('OA209',0.3);
%airfoil.steady = SteadyCurve(static_data(:,1),static_data(:,2));
load('static_corr')
airfoil.steady = SteadyCurve(mA_corr,mCl_corr);
%% Dynamic data
nr=13;

data = load(pressuredata(nr));


% Drag missing, compute the pressure drag 
CD = zeros(size(data.Cl));
for k=1:length(data.Cl)
    CD(k) = sum(data.Cp(k,:)*data.xk);
end

load('dynamic_corr')
% assuming CL = CN
pitching = PitchingMotion('alpha',Alpha_corr,'CN',Cl_corr.*cosd(Alpha_corr),'k',LB(nr).k,'freq',LB(nr).fosc,'V',50);
pitching.setSinus(airfoil,deg2rad(LB(nr).alpha_0),deg2rad(LB(nr).alpha_1),LB(nr).fosc,LB(nr).FS);
pitching.setName('simcos')
pitching.setCNsteady(airfoil.steady)

% set static slope and zero lift
airfoil.steady.computeSlope();
airfoil.steady.CN0 = interp1(airfoil.steady.alpha,airfoil.steady.CN,0,'linear','extrap');

airfoil.steady.fitKirchhoff();
airfoil.steady.plotKirchhoff(); 


%% Define Beddoes-Leishman model
Tp = 3;
Tf = 1;
Tv = 1;
pitching.BeddoesLeishman(airfoil,Tp,Tf,Tv);

%% Plot results
figure
plot(pitching.alpha(1:length(pitching.CN_LB)),pitching.CN_LB,'DisplayName','C_{N,LB}','LineWidth',2)
hold on
plot(pitching.alpha,pitching.CN,'DisplayName','C_{N,xp}','LineWidth',2)
plot(pitching.alpha,pitching.alpha*2*pi*pi/180,'r--','DisplayName','2\pi\alpha')
legend('Location','SouthWest')
xlabel('\alpha (°)')
ylabel('C_N')
grid on