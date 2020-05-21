close all
clear all
clc

z = tf('z');
G = 1/(1-z^(-1));

U = z*sin(1)/(z^2*2*z*cos(1)+1);

bode(G)
grid on

Y = G*U; 


