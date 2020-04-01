close all
clear all
clc

% experimental parameters
M = 0.12;
V = M*343;
c = 0.55;
r = 0.04;

% numerical parameters
n = 200;

% alpha time-evolution
alphadot = rad2deg(r*2*V/c);
t = linspace(0,35/alphadot,n);
alpha = alphadot*t;

% parameters for Eq.5
s = 2*V/c*t;
Talpha = 4;

% computation of alpha_lag
dalpha = diff(alpha);
alpha_lag = zeros(size(alpha));
for k = 1:length(dalpha)
    alpha_lag(k+1) = alpha_lag(k) + dalpha(k)*(1-exp(-s(k)/Talpha));
end

% plot results
figure 
plot(t,alpha,'DisplayName','\alpha')
hold on 
plot(t,alpha_lag,'DisplayName','\alpha''')
grid on
xlabel('t (s)')
ylabel('\alpha')
legend('Location','NorthWest')