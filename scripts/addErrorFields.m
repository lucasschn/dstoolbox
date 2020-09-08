% This script is made to add some error fields to a res structure
% Author : Lucas Schneeberger
% Date : 05.09.2020

close all
clear all
clc

filename = 'res25uniform';
load(fullfile('..','data','paramsweep',filename))

for k=1:length(res)
    res(k).errPeakLoc = res(k).SmaxCN_LB - res(k).SmaxCN;
end

save(fullfile('..','data','paramsweep',filename),'res')