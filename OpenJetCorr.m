%% Wind Tunnel Wall Corrections - AGARD-AG-336
% (Correction des effets de paroi en soufflerie)
% page 103

close all; clear; clc;
run('data/2008_simcos/matlab/labbook.m')
I = param;
I.indx = 13;

load('data/2008_simcos/matfiles/static_pressure.mat')
%% Correction for open jet wind tunnels - static
delta = 0.16;
S_wing = 0.3;
S_windt = 0.7875;

del_A = rad2deg((1 + 0.3)*((mCl*delta*S_wing)./S_windt));
mA_corr = (mA + del_A); 
mCl_corr = mCl + mCl.*cosd(del_A);

save('static_corr','mA_corr','mCl_corr')

alpha_temp = 1:1:30;
CLa = 2*pi*deg2rad(alpha_temp);
%% Figuring
FontSizeAx = 20;
FontSizeLb = 24;
afFigurePosition = [1 1 25 20];
axespos = [0.21 0.27 0.75 0.65];
ylabelpos = [-0.16 0.4];
xlabelpos = [0.5 -0.14];
col = cool(6);

figure;
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
plot(mA,mCl,'color',col(1,:),'linewidth', 3)
hold on, plot(mA,mCl_corr,'color',col(3,:),'linewidth', 3)
hold on, plot(alpha_temp,CLa,'color',col(6,:),'linewidth', 1.5, 'Linestyle', '--')
legend ({'Static','Static Corrected', '$$2 \pi \alpha$$'}, 'location', 'NorthWest','Interpreter','latex')
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx,'TickLabelInterpreter','latex','DefaultFigureWindowStyle','docked')
xlabel('$$\alpha$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_L$$','Interpreter', 'LaTeX','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos);

%% Correction for open jet wind tunnels - dynamic
load(['data/2008_simcos/matfiles/pressure_',LB(I.indx).dpt,'.mat'])

del_A_dyn = rad2deg((1 + 0.3)*((Cl*delta*S_wing)./S_windt));
Cl_corr = Cl + Cl.*cosd(del_A_dyn);
Alpha_corr = (Alpha + del_A_dyn); 

save('dynamic_corr','Alpha_corr','Cl_corr')

figure;
set(gcf, 'Units', 'centimeters','PaperPositionMode', 'auto','Position', afFigurePosition);
set(gcf,'DefaultAxesFontSize',FontSizeAx,'DefaultAxesFontName','TimesNewRoman','DefaultAxesGridLineStyle','-.','DefaultAxesLineWidth',2,'DefaultAxesFontWeight','Normal')
plot(Alpha,Cl,'color',col(1,:),'linewidth', 3)
hold on, plot(Alpha,Cl_corr,'color',col(3,:),'linewidth', 3)
hold on, plot(alpha_temp,CLa,'color',col(6,:),'linewidth', 1.5, 'Linestyle', '--')
legend ({'Dynamic','Dynamic Corrected', '$$2 \pi \alpha$$'}, 'location', 'NorthWest','Interpreter','latex')
set(gca, 'FontName', 'TimesNewRoman', 'FontSize', FontSizeAx,'TickLabelInterpreter','latex','DefaultFigureWindowStyle','docked')
xlabel('$$\alpha$$','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 0,'Units', 'Normalize','Position',xlabelpos,'Interpreter', 'LaTeX');
ylabel('$$C_L$$','Interpreter', 'LaTeX','FontName', 'TimesNewRoman','FontUnit', 'points','FontSize', FontSizeLb,'FontWeight', 'normal','Rotation', 90,'Units', 'Normalize','Position',ylabelpos);