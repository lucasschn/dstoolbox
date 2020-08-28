close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
run('/Users/lucas/src/codes_smarth/labbook.m')
% fy is normal to the airfoil, pointing downwards, fx is parallel to the
% airfoil's chord

n = 26;
c = 0.15; % chord length
s = 0.6; % span
rho = 1000;
U = 0.5;

% /!\ ms115 is an old static curve with alpha<25, higher AoA are missing.
% Use calibrate_polar.m instead.
load('../data/2019_SH/20191009/loads/ms115mpt001.mat');
fdata = [out(2).t out(2).forces(:,1) out(2).forces(:,2)];
[t_NI, helpi_NI] = unique(fdata(:,1)); % time vector for the force data
calib.fx = mean(fdata(helpi_NI,2));
calib.fy = mean(fdata(helpi_NI,3));

alpha = zeros(1,n);
Cn = zeros(1,n);
Cl = zeros(1,n);
Cd = zeros(1,n);
for k=1:n
    load(sprintf('../data/2019_SH/20191009/loads/ms116mpt%03.0f.mat',k))
    alpha(k) = motangle.angle.position;
    fdata = [out(2).t out(2).forces(:,1) out(2).forces(:,2)];
    [t_NI, helpi_NI] = unique(fdata(:,1)); % time vector for the force data
    fx = fdata(helpi_NI,2);
    fy = fdata(helpi_NI,3);
    L = -cosd(alpha(k)).*(fy-calib.fy) + sind(alpha(k)).*(fx-calib.fx);
    D = -sind(alpha(k)).*(fy-calib.fy) + cosd(alpha(k)).*(fx-calib.fx);
    A = c*s; % [m^2]
    q = 0.5*rho*U^2*A;
    Cl(k) = mean(L)/q;
    Cd(k) = mean(D)/q;
    Cn(k) = -mean(fy-calib.fy)/q;
end

figure
plot(alpha,Cl,'DisplayName','C_l')
hold on 
plot(alpha,Cd,'DisplayName','C_d')
plot(alpha,Cn,'DisplayName','C_n')
plot(alpha,2*pi*deg2rad(alpha),'DisplayName','2\pi\alpha')
grid on 
legend('Location','NorthWest')
xlabel('\alpha (°)')
ylabel('C_f')

save('static_flatplate','alpha','Cn','Cd','Cl')
