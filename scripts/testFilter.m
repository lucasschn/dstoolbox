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


x_filtered = myFilterTwice(x_noisy,fs);


figure

subplot(311)
plot(t,x,'LineWidth',4)
title('a)')
ylabel('x')
grid on 
ax = gca; 
ax.FontSize = 20; 

subplot(312)
plot(t,x_noisy,'LineWidth',2)
title('b)')
ylabel('x')
grid on 
ax = gca; 
ax.FontSize = 20; 

subplot(313)
plot(t,x,'LineWidth',4,'DisplayName','original')
hold on
plot(t,x_filtered,'--','LineWidth',4,'DisplayName','retrieved')
legend show
title('c)')
xlabel('t (s)')
ylabel('x')
grid on 
ax = gca; 
ax.FontSize = 20; 


saveas(gcf,'../fig/filter_test_without','png')