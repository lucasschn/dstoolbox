close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')

t = 1:.1:10;

dt = diff(t);
D = zeros(length(dt));
for k=2:length(dt)
    D(k) = D(k-1)*exp(-dt(k))+ 0.1;
end

figure 
plot(t(1:length(D)),D)
xlabel('t (s)')
ylabel('C_N')
grid on

saveas(gcf,'fig/deficiency_example','png')