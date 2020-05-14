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
airfoil.steady = SteadyCurve(mA,mCl_corr);
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
pitching = PitchingMotion('alpha',Alpha,'CN',Cl_corr.*cosd(Alpha),'k',LB(nr).k,'freq',LB(nr).fosc,'V',50,'Ts',1/LB(nr).FS);
pitching.setSinus(airfoil,deg2rad(LB(nr).alpha_0),deg2rad(LB(nr).alpha_1),LB(nr).fosc*2*pi);
pitching.setName('simcos')
pitching.setCNsteady(airfoil.steady)

% set static slope and zero lift
airfoil.steady.computeSlope();
airfoil.steady.setCN0(interp1(airfoil.steady.alpha,airfoil.steady.CN,0,'linear','extrap'));
airfoil.steady.fitKirchhoff();
airfoil.steady.plotKirchhoff(); 


%% Define Beddoes-Leishman model
Tp = 1.7;
Tf = 3;
Tv = 6;
pitching.BeddoesLeishman(airfoil,Tp,Tf,Tv,'analytical');

%% Plot results
figure
plot(pitching.alpha,pitching.CN,'DisplayName','C_{N,xp}','Color','b','LineWidth',2)
hold on
plot(pitching.alpha(1:length(pitching.CN_LB)),pitching.CN_LB,'DisplayName','C_{N,LB}','Color','r','LineWidth',2)
plot(pitching.alpha,pitching.alpha*2*pi*pi/180,'r--','DisplayName','2\pi\alpha')
legend('Location','SouthWest')
xlabel('\alpha (°)')
ylabel('C_N')
grid on