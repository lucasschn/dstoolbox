close all
clear all
clc


alpha = 0:.1:20;
steady = SteadyCurve(alpha,2*pi*alpha,15);

S1 = 0.2:.2:1.6;
CNk = zeros(length(S1),length(alpha));

for k=1:length(S1)
    CNk(k,:) = Kirchhoff(steady,alpha,[S1(k),1]);
end

figure
hold on
for k=1:length(S1)
    plot(alpha,CNk(k,:),'DisplayName',sprintf('S1=%0.1f',S1(k)))
end
grid on 
legend show
xlabel('\alpha (°)')
ylabel('C_N')