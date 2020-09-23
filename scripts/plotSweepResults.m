% This script is made to produce scatter plots and histogram plots in order
% to analyse the results of paramsweep.m.
% For some error variables to be correctly set, you need to previously have run the
% script addErrorFields on the results mat-file.
%
% Author : Lucas Schneeberger
% Date : 28.08.2020

close all
clear all
clc

load(fullfile('..','data','paramsweep','res25uniform4'))

plotHistogram(res,'Tv','errSecondPeakLoc',1)

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
    grid minor
    setDataTip(s,res_adot)
    
    if length(res_adot)>5e3
        s.SizeData = 10;
    end
    xlabel(getLabelString(varx))
    ylabel(getLabelString(vary))
    ax = gca;
    ax.FontSize = 20;
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

function plotAllRates(res,varx,vary)
color = nan(size(res));
adot_min = min(cat(1,res.alphadot));
adot_max = max(cat(1,res.alphadot));
for k = 1:length(res)
    color(k) = (res(k).alphadot-adot_min)/(adot_max-adot_min);
end

x = eval(sprintf('cat(1,res.%s)',varx));
y = eval(sprintf('cat(1,res.%s)',vary));

figure
s = scatter(x,y,'filled');
grid on
xlabel(getLabelString(varx))
ylabel(getLabelString(vary))
ax = gca;
    ax.FontSize = 20;
s.CData = color;
c = colorbar;
c.Ticks = linspace(0,1,length(unique(cat(1,res.alphadot))));
c.Label.String = 'pitch rate (°/s)';
c.TickLabels = unique(color*(adot_max-adot_min)+adot_min);
end

function plotHistogram(res,varx,vary,threshold)
if contains(varx,'+')
    varcell = split(varx,'+');
    for k=1:length(varcell)
        x(k,:) = eval(sprintf('cat(1,res.%s)',varcell{k}));
    end
    x = sum(x,1);
else
    x = eval(sprintf('cat(1,res.%s)',varx));
end

y = eval(sprintf('cat(1,res.%s)',vary));

bin_vect = 0:0.1:max(x);

figure
h = histogram(x,bin_vect,'DisplayName','total');
hold on
% take the absolute value of vary small than the threshold
% this is useful for peak location error and does not impact overall error
% (always positive)
histogram(x(abs(y) < threshold),bin_vect,'DisplayName',sprintf('|%s| < %g',getLabelString(vary),threshold))
grid on
ax = gca; 
ax.FontSize=20;
xlabel(getLabelString(varx))
ylabel('nb of samples')
axis([0, max(x), 0, Inf])
legend('Location','NorthWest')
ax1 = gca;

figure
rel_freq = histcounts(x(abs(y)<threshold),bin_vect)./histcounts(x,bin_vect);
b = bar(0.05:0.1:max(x)-0.05,rel_freq,1,'DisplayName',sprintf('|%s| < %g',getLabelString(vary),threshold));
b.EdgeColor = 'black';
b.FaceAlpha = h.FaceAlpha;
grid on
legend('Location','NorthEast')
xlim([bin_vect(1) bin_vect(end)])
xlabel(getLabelString(varx))
ylabel('relative nb of samples (%)')
ax2 = gca;
ax2.XTick = ax1.XTick;
ax2.FontSize = 20; 
yticklabels(num2cell(ax2.YTick*100))
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
        label = 'timing of C_N^{LB} peak';
    case 'SmaxCNf'
        label = 'timing of C_N^f peak';
    case 'maxCNf'
        label = 'height of C_N^f peak';
    case 'SmaxCNv'
        label = 'timing of C_N^v peak';
    case 'maxCNv'
        label = 'height of C_N^v peak';
    case 'err'
        label = 'mean square error';
    case 'errPeakLoc'
        label = 'error on the primary peak timing (convectime time units)';
    case 'errPeakHeight'
        label = 'error on the primary peak height';
    case 'errCNk_PeakLoc'
        label = 'error on Kirchhoff peak timing (convectime time units)';
    case 'errCNk_PeakHeight'
        label = 'error on Kirchhoff peak height';
    case 'errCNf_PeakLoc'
        label = 'error on Kirchhoff w/ added mass peak timing (convectime time units)';
    case 'errCNf_PeakHeight'
        label = 'error on Kirchhoff w/ added mass peak height';
    case 'errCNv_PeakLoc'
        label = 'error on vortex lift peak timing (convectime time units)';
    case 'errCNv_PeakHeight'
        label = 'error on vortex lift peak height';
    case 'errFirstPeakLoc'
        label = 'error on first total lift peak timing';
    case 'errFirstPeakHeight'
        label = 'error on first total lift peak height';
    case 'errSecondPeakLoc'
        label = 'error on second total lift peak timing (if any)';
    case 'errSecondPeakHeight'
        label = 'error on second total lift peak height';
    otherwise
        label = var;
end
end

function setDataTip(scatterplot,res)
scatterplot.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Sample',1:length(res));
scatterplot.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Tp',cat(1,res.Tp));
scatterplot.DataTipTemplate.DataTipRows(5) = dataTipTextRow('Tf',cat(1,res.Tf));
scatterplot.DataTipTemplate.DataTipRows(6) = dataTipTextRow('Tv',cat(1,res.Tv));
scatterplot.DataTipTemplate.DataTipRows(7) = dataTipTextRow('Tvl',cat(1,res.Tvl));
end
