close all
clear all 
clc 
set(0,'DefaultFigureWindowStyle','docked')

M = 0.12;
V = M*343;

%% NACA 23012
naca23012 = Airfoil('NACA 23012',0.55);

stall23012 = readmatrix('data/sheng_naca23012.csv');
CNalpha23012 = readmatrix('data/CN_alpha.csv');
CNalphalag23012 = readmatrix('data/CN_alpha_lag.csv');

r = stall23012(:,1);
alpha_ds = stall23012(:,2);

for k=1:length(r)
    assignin('base',sprintf('ms%d',k),RampUpMotion('r',r(k),'V',V,'alpha_CConset',alpha_ds(k)));
    eval(sprintf('ms%d.setName()',k)) % defines the name property from the name of the instance
    eval(sprintf('ms%d.setalphadot(%3.4f)',k,r(k)*2*V/naca23012.c))
end

% Define alpha_ds0 & compute Talpha
naca23012.Sheng(ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10,ms11,ms12,ms13,ms14,ms15,ms16,ms17,ms18,ms19,ms20,ms21,ms22,ms23,ms24,ms25,ms26,ms27);

r = 0.025;
alphadot = r*2*V/naca23012.c;
alpha = CNalpha23012(:,1);
CN_alpha = CNalpha23012(:,2);
alpha_lag = CNalphalag23012(:,1);
CN_alpha_lag = CNalphalag23012(:,2);
t = alpha/alphadot;
% alpha_lag_alpha = interp1(CN_alpha_lag,alpha_lag,CN_alpha);

fig1 = figure;
plot(alpha,CN_alpha,'DisplayName','C_N(\alpha)')
hold on 
plot(alpha_lag,CN_alpha_lag,'DisplayName','C_N(\alpha'')')
xlabel('\alpha (°)')
ylabel('C_N')
grid on

% fig2 = figure; 
% plot(t,alpha,'DisplayName','\alpha')
% hold on
% plot(t,alpha_lag,'DisplayName','\alpha''')
% xlabel('t (s)')
% ylabel('\alpha (°)')
% legend show
% grid on

%% NACA 0012
naca0012 = Airfoil('NACA 0012',0.55); 

data0012 = readmatrix('data/sheng_naca0012.csv');

r = data0012(:,1);
alpha_ds = data0012(:,2);

for k=1:length(r)
    assignin('base',sprintf('ms%d',k),RampUpMotion('r',r(k),'V',V,'alpha_CConset',alpha_ds(k),'alphadot',r(k)*2*V/naca23012.c));
    eval(sprintf('ms%d.setName()',k)) % defines the name property from the name of the instance
    eval(sprintf('ms%d.setalphadot(%3.4f)',k,r(k)*2*V/naca23012.c))
end

% Define alpha_ds0 & compute Talpha
naca0012.Sheng(ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10,ms11,ms12,ms13,ms14,ms15,ms16,ms17,ms18,ms19,ms20,ms21,ms22,ms23);
% How to automate this argument passing? Sheng should be accepting a vector
% of ramps instead?
