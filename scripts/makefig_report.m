% This script is made to save all figures needed for the report
% Author : Lucas Schneeberger
% Date : 06.06.2020
close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')

run(fullfile('..','labbook.m'))

%% Defines paths 

%% Define the airfoil and the associated steady curve

airfoil = Airfoil('flatplate',0.15);
airfoil.r0 = 0.04;
static = load(fullfile('..','static_flatplate'));
airfoil.steady = SteadyCurve(static.alpha,static.CN,13.5);

%% Setting up the non-filterd ramp

ramp = loadRamp(22,false);
ramp.setPitchRate(airfoil);
ramp.findExpOnset();

%% Setting up the filtered ramps

c = [22,67,75];

for k=1:length(c)
    assignin('base',sprintf('ramp%d_filt',k),loadRamp(c(k),true));
    evalin('base',sprintf('ramp%d_filt.setPitchRate(airfoil)',k));
    evalin('base',sprintf('ramp%d_filt.findExpOnset()',k));
end

%% Make figures of the plain experimental data

ramp.plotAlpha('normal')
% saveas(gcf,'ramp_exmaple','png')
ramp.plotCN('convectime')
% saveas(gcf,'../fig/CN','png')
ramp.plotCC('convectime')
% saveas(gcf,'../fig/CC','png')

for k=1:length(c)
    evalin('base',sprintf('ramp%d_filt.plotCN(''convectime'')',k))
    % saveas(gcf,'../fig/CN_filt','png')
    evalin('base',sprintf('ramp%d_filt.plotCC(''convectime'');',k))
    % saveas(gcf,'../fig/CC_filt','png')
end

%% Test LB

Tp = 3.6; Tf = 1; Tv = 2; Tvl = 0.5;

save_LB = false;

ramp1_filt.BeddoesLeishman(airfoil,Tp,Tf,Tv,Tvl,'experimental')
ramp1_filt.plotLB('convectime')
if save_LB
    saveas(gcf,fullfile('..','fig','CN_LB_r020'),'png')
end

CN1_LB = ramp1_filt.CN_LB;

ramp2_filt.BeddoesLeishman(airfoil,Tp,Tf,Tv,Tvl,'experimental')
ramp2_filt.plotLB('convectime')
if save_LB
    saveas(gcf,fullfile('..','fig','CN_LB_r026'),'png')
end

CN2_LB = ramp2_filt.CN_LB;

ramp3_filt.BeddoesLeishman(airfoil,Tp,Tf,Tv,Tvl,'experimental')
ramp3_filt.plotLB('convectime')
if save_LB
    saveas(gcf,fullfile('..','fig','CN_LB_r092'),'png')
end

CN3_LB = ramp3_filt.CN_LB;

ramp1_filt.plotSeparation(airfoil,'convectime',1)

fprintf('LB: r=.020, r=.026, r=0.092 \n')
fprintf('%.1f%%, %.1f%%, %.1f%% \n',max(CN1_LB)/max(ramp1_filt.CN)*100-100,max(CN2_LB)/max(ramp2_filt.CN)*100-100,max(CN3_LB)/max(ramp3_filt.CN)*100-100) 
%% Test Sheng-LB

Tf = 0; Tv = 0.5; Tvl = 3;

save_shengLB = false;

ramp1_filt.BLSheng(airfoil,Tf,Tv,Tvl,'experimental')
ramp1_filt.plotShengLB(airfoil)
if save_shengLB
    saveas(gcf,'../fig/CN_ShengLB_r020','png')
end
CN1_ShengLB = ramp1_filt.CN_LB;

Tf = 0; Tv = 0.5; Tvl = 2; % Tvl higher than 3 makes no difference

ramp2_filt.BLSheng(airfoil,Tf,Tv,Tvl,'experimental')
ramp2_filt.plotShengLB(airfoil)
if save_shengLB
    saveas(gcf,'../fig/CN_ShengLB_r026','png')
end
CN2_ShengLB = ramp2_filt.CN_LB;

Tf = 0; Tv = 3; Tvl = 2;

ramp3_filt.BLSheng(airfoil,Tf,Tv,Tvl,'experimental')
ramp3_filt.plotShengLB(airfoil)
if save_shengLB
    saveas(gcf,'../fig/CN_ShengLB_r092','png')
end
CN3_ShengLB = ramp3_filt.CN_LB;

ramp2_filt.plotSeparation(airfoil,'convectime',0)
if save_shengLB
    saveas(gcf,'../fig/f_ShengLB_r026','png')
end

fprintf('Sheng-LB: r=.020, r=.026, r=0.092 \n')
fprintf('%.1f%%, %.1f%%, %.1f%% \n',max(CN1_ShengLB)/max(ramp1_filt.CN)*100-100,max(CN2_ShengLB)/max(ramp2_filt.CN)*100-100,max(CN3_ShengLB)/max(ramp3_filt.CN)*100-100) 
%% Test Expfit-LB

save_expfit = false;

Tf = 0; Tv = 0.5; Tvl = 3;

ramp1_filt.BLexpfit(airfoil,Tf,Tv,Tvl,'experimental')
ramp1_filt.plotLBExpfit(airfoil,CN1_ShengLB)
if save_expfit
    saveas(gcf,'../fig/CN_ExpfitLB_r020','png')
end

CN_ExpfitLB = ramp1_filt.CN_LB;

Tf = 0; Tv = 0.5; Tvl = 2;

ramp2_filt.BLexpfit(airfoil,Tf,Tv,Tvl,'experimental')
ramp2_filt.plotLBExpfit(airfoil,CN2_ShengLB)
if save_expfit
    saveas(gcf,'../fig/CN_ExpfitLB_r026','png')
end

Tf = 0; Tv = 3; Tvl = 2;

ramp3_filt.BLexpfit(airfoil,Tf,Tv,Tvl,'experimental')
ramp3_filt.plotLBExpfit(airfoil,CN3_ShengLB)
if save_expfit
    saveas(gcf,'../fig/CN_ExpfitLB_r092','png')
end

%% Comparison of models
n = length(CN1_LB);

figure
plot(ramp1_filt.S,ramp1_filt.CN,'LineWidth',2,'DisplayName','exp')
hold on
plot(ramp1_filt.S(1:n),CN1_LB,'LineWidth',2,'DisplayName','LB')
plot(ramp1_filt.S(1:n),CN1_ShengLB,'LineWidth',2,'DisplayName','Sheng-LB')
plot(ramp1_filt.S(1:n),CN_ExpfitLB,'LineWidth',2,'DisplayName','LBExpfit')
grid on
xlabel('t_c')
ylabel('C_N')
ax = gca;
ax.FontSize=20;
legend('FontSize',20,'Location','SouthEast')