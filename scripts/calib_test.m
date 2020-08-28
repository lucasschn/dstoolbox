close all
clear all
clc

load('../data/2019_SH/Postprocessing/matfiles/loads/ms012mpt001_loads.mat','raw','zero','avg')


figure
plot(zero.Cl)