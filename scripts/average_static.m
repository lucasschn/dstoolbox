close all
clear all
clc

for k=1:5
    assignin('base',sprintf('s%d',k),eval(sprintf('load(''../data/2020_SH/20200510/calibration/static_flatplate%d.mat'')',k)));
end

figure
for k=1:5
    eval(sprintf('plot(s%d.alpha,s%d.CN,''DisplayName'',''s%d'')',k,k,k))
    hold on
end
grid on
xlabel('\alpha °')
ylabel('C_N')
legend('Location','SouthEast')

if isequal(s1.alpha,s2.alpha,s3.alpha,s4.alpha,s5.alpha)
    savg.alpha = s1.alpha;
    savg.CN = zeros(size(s1.alpha));
    for k=1:5
        CN = eval(sprintf('s%d.CN',k));
        savg.CN = (savg.CN*(k-1)+CN)/k;
    end
else
    error('Static curves have different alpha ranges.')
end

plot(savg.alpha,savg.CN,'DisplayName','avg','LineWidth',2)

save('static_flatplate','-struct','savg')