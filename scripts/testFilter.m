% Simple script to the test the performance of myFilter function in
% filtering a 

close all
clear all
clc 


Ts = 0.01;
fs = 1/Ts;

t = 0:Ts:10; 


f1 = 1; %Hz 
f2 = 35; %Hz 
f3 = 50; %Hz 


x = sin(2*pi*f1*t);
x_noisy = x + sin(2*pi*f2*t+rand(1))+sin(2*pi*f3*t+rand(1));


x_filtered = myFilterTwice(x_noisy,fs,0);
x_movmean = myFilterTwice(x_noisy,fs,1);

figure
plot(t,x,'LineWidth',4)
xlabel('t (s)')
grid on 
ax = gca; 
ax.FontSize = 20; 
saveas(gcf,'../fig/1Hz','png')

figure
plot(t,x_noisy,'LineWidth',2)
xlabel('t (s)')
grid on 
ax = gca; 
ax.FontSize = 20; 
saveas(gcf,'../fig/noisy_1Hz','png')

figure
plot(t,x,'LineWidth',4,'DisplayName','original')
hold on
plot(t,x_filtered,'--','LineWidth',4,'DisplayName','retrieved')
legend show
xlabel('t (s)')
grid on 
ax = gca; 
ax.FontSize = 20; 
saveas(gcf,'../fig/filter_test_without','png')

figure
plot(t,x,'LineWidth',4,'DisplayName','original')
hold on
plot(t,x_movmean,'-','LineWidth',4,'DisplayName','retrieved')
legend show
xlabel('t (s)')
grid on 
ax = gca; 
ax.FontSize = 20; 
saveas(gcf,'../fig/filter_test_with','png')