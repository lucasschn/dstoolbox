close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')

%% Import ramp-up motion

leishman = load('originalLB');
sheng = load('modLBsheng');

%% Obtain CN for both models

[CN_LB,alpha] = getCN(leishman);
CN_SH = getCN(sheng);

figure
plot(alpha(1:length(CN_LB)),CN_LB,'DisplayName','Leishman')
hold on 
plot(alpha(1:length(CN_LB)),CN_SH,'DisplayName','Sheng')
plot(alpha,leishman.ms013.CN,'DisplayName','exp')
grid on
xlabel('alpha (°)')
ylabel('C_N')
legend('Location','NorthWest','FontSize',20)

function [CN,alpha] = getCN(dataset)
CN = dataset.ms013.CN_LB;
alpha = dataset.ms013.alpha;
end