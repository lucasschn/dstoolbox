close all
clear all
clc

%% Load 1deg curves
for k=1:5
    assignin('base',sprintf('s%d',k),eval(sprintf('load(''../data/2020_SH/20200510/calibration/static_flatplate%d.mat'')',k)));
end

%% Load 0.5deg curve
highres = load('../upstrokepolar');

highres.polarforces.Cn = highres.polarforces.Cl.*cosd(highres.polarforces.alpha)+highres.polarforces.Cd.*sind(highres.polarforces.alpha);

%% Compute avg CN
if isequal(s1.alpha,s2.alpha,s3.alpha,s4.alpha,s5.alpha)
    savg.alpha = s1.alpha;
    savg.std = zeros(size(s1.alpha));
    savg.CN = zeros(size(s1.alpha));
    for k=1:5
        CN = eval(sprintf('s%d.CN',k));
        savg.CN = (savg.CN*(k-1)+CN)/k;
    end
    for k=1:length(s1.alpha)
        savg.std(k) = std([s1.CN(k) s2.CN(k) s3.CN(k) s4.CN(k) s5.CN(k)]);
    end
else
    error('Static curves have different alpha ranges.')
end

%% Plot figure

xconf = [savg.alpha, flip(savg.alpha)];
yconf = [savg.CN+savg.std, flip(savg.CN-savg.std)];

figure
p = fill(xconf,yconf,[1 .8 .8]);
p.EdgeColor = 'none';
p.DisplayName = 'std';
hold on 
plot(savg.alpha,savg.CN,'DisplayName','avg','LineWidth',2)
if 0
    for k=1:5
        eval(sprintf('plot(s%d.alpha,s%d.CN,''DisplayName'',''s%d'')',k,k,k))
    end
end
if 1
    plot(highres.polarforces.alpha,highres.polarforces.Cn,'DisplayName','highres')
end
grid on
ax = gca; 
ax.FontSize = 20; 
xlabel('\alpha (°)')
ylabel('C_N (-)')
legend('Location','SouthEast')



%% Save avg curve and figure
%save('../static_flatplate','-struct','savg')
%saveas(gcf,'../fig/static_flatplate','png')
