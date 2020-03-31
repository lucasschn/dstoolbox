close all
clear all 
clc 


M = 0.12;
V = M*343;

%% NACA 23012
naca23012 = Airfoil('NACA 23012',0.55);

data23012 = readmatrix('data/sheng_naca23012.csv');

r = data23012(:,1);
alpha_ds = data23012(:,2);

for k=1:length(r)
    assignin('base',sprintf('ms%d',k),RampUpMotion(sprintf('ms%d',k)));
    evalin('base',sprintf('ms%d.r=%.3e;',k,r(k)));
    evalin('base',sprintf('ms%d.V=%.2e;',k,V));
    evalin('base',sprintf('ms%d.alpha_CConset=%.2f;',k,alpha_ds(k)));
    evalin('base',sprintf('ms%d.setalphadot(%.2f);',k,r(k)*2*V/naca23012.c));
end


% Define alpha_ds0 & compute Talpha
naca23012.Sheng(ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10,ms11,ms12,ms13,ms14,ms15,ms16,ms17,ms18,ms19,ms20,ms21,ms22,ms23,ms24,ms25,ms26,ms27);

%% NACA 0012
naca0012 = Airfoil('NACA 0012',0.55); 

data0012 = readmatrix('data/sheng_naca0012.csv');

r = data0012(:,1);
alpha_ds = data0012(:,2);

for k=1:length(r)
    assignin('base',sprintf('ms%d',k),RampUpMotion(sprintf('ms%d',k)));
    evalin('base',sprintf('ms%d.r=%.3e;',k,r(k)));
    evalin('base',sprintf('ms%d.V=%.2e;',k,V));
    evalin('base',sprintf('ms%d.alpha_CConset=%.2f;',k,alpha_ds(k)));
    evalin('base',sprintf('ms%d.setalphadot(%.2f);',k,r(k)*2*V/naca0012.c));
end

% Define alpha_ds0 & compute Talpha
naca0012.Sheng(ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10,ms11,ms12,ms13,ms14,ms15,ms16,ms17,ms18,ms19,ms20,ms21,ms22,ms23);
% How to automate this argument passing?
