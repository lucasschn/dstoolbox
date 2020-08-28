clear all
close all
clc
 
C = linspecer(12,'qualitative'); 

doplot = [1,1,1,1,1,1,1,1,1,1,1,1];
%%%%%%%%%[1,2,3,4,5,6,7,8,9,10,11,12
figure
ax = gca; 
for k=1:11
    % plot(ax,1:10,rand(1,10))
    if doplot(k)
        set(ax,'NextPlot','add', 'ColorOrder',C);
        plot(ax,0:.1:10,sin(k/10.*(0:.1:10)),'LineWidth',2)        
    end    
end
legend(ax,'show')