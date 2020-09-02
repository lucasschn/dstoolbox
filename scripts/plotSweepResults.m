% This script is made to plot the results of paramsweep.m
% Author : Lucas Schneeberger
% Date : 28.08.2020

close all
clear all
clc

load(fullfile('..','data','paramsweep','res'))
y
plotOneRate(res,25,'Tp','CNv','Tf')

function plotOneRate(res,rate,varx,vary,color_var)
res_adot = res(cat(1,res.alphadot)==rate);
if isempty(res_adot)
    warning('The specified pitch rate could not be found in the results structure.')
else
    for k=1:length(res_adot)
        x(k) = eval(sprintf('res_adot(k).%s',varx));
        y(k) = eval(sprintf('res_adot(k).%s',vary));
    end
    figure
    s = scatter(x,y,'filled');
    grid on
    xlabel(getLabelString(varx))
    ylabel(getLabelString(vary))
    title(sprintf('alphadot = %.2f',rate))
end

if nargin > 4
    color = nan(size(res));
    cvar_min = eval(sprintf('min(cat(1,res_adot.%s))',color_var));
    cvar_max = eval(sprintf('max(cat(1,res_adot.%s))',color_var));
    for k = 1:length(res_adot)
        color(k) = (eval(sprintf('res_adot(k).%s',color_var))-cvar_min)/(cvar_max-cvar_min);
    end
    
    s.CData = eval(sprintf('cat(1,res_adot.%s)',color_var));
    
    c = colorbar;
    c.Ticks = linspace(cvar_min,cvar_max,length(unique(cat(1,sprintf('res.%s',color_var)))));
    c.Label.String = sprintf('%s',getLabelString(color_var));
end

end

function plotAllRates(res,x,y)
color = nan(size(res));
adot_min = min(cat(1,res.alphadot));
adot_max = max(cat(1,res.alphadot));
for k = 1:length(res)
    color(k) = (res(k).alphadot-adot_min)/(adot_max-adot_min);
end

figure
s = scatter(cat(1,res.Tp), cat(1,res.maxCN),'filled');
grid on
xlabel('T_p')
ylabel('max C_N')
s.CData = color;
c = colorbar;
c.Ticks = linspace(0,1,length(unique(cat(1,res.alphadot))));
c.Label.String = 'pitch rate (°/s)';
c.TickLabels = unique(color*(adot_max-adot_min)+adot_min);
end

function label = getLabelString(var)
switch var
    case 'Tp'
        label = 'T_p';
    case 'Tf'
        label = 'T_f';
    case 'Tv'
        label = 'T_v';
    case 'Tvl'
        label = 'T_{vl}';
    case 'maxCN_LB'
        label = ' max C_N^{LB}';
    case 'maxCNf'
        label = ' max C_N^{f}';
    case 'maxCNk'
        label = ' max C_N^{k}';
    case 'SmaxCN_LB'
        label = 't_{C,maxLB}';
    case 'maxCNv'
        label = 'max C_N^v';
    otherwise
        label = var;
end
end