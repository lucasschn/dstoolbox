close all
clear all
clc

data = csvread('../data/St_flatplate_chen1996.csv');

alpha = data(:,1);
St = data(:,2);

figure
plot(alpha,St,'d')
grid on 
xlabel('\alpha (°)')
ylabel('St')
