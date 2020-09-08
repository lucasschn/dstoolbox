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
    if res(k).Tv == 0 || res(k).Tvl  == 0 % means there is no vortex lift
                res(k).SmaxCNv = NaN; 
                warning('Peak location has been replaced by NaN.')
    end  
    if res(k).SmaxCNk > 30 
        res(k).SmaxCNk = NaN;
        warning('Peak location has been replaced by NaN.')
    end
    
end

save(fullfile('..','data','paramsweep',filename),'res')