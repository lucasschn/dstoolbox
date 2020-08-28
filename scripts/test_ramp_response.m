close all
clear all
clc

t = 0:0.01:1;
K = 2*pi; 
tp = 4; 

CNprime = K*(t.^2/tp-t.^3/tp.^2+t.^4/tp^3-t.^5/tp^4+t.^6/tp^5-t.^7/tp^6-t.^8/tp^7);

plot(t,CNprime)
grid on 