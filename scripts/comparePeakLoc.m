close all
clear all
clc

load(fullfile('..','data','paramsweep','res25uniform3'))

CNf_PeakLoc = cat(1,res.SmaxCNf);
CNv_PeakLoc = cat(1,res.SmaxCNv);
CNtot_PeakLoc = cat(1,res.SmaxCN_LB);
dist = CNv_PeakLoc - CNf_PeakLoc;
[dist_sorted,sorted_by_dist] = sort(dist);

figure
scatter(,CNf_PeakLoc(sorted_by_dist(1:50:end)),'DisplayName','Kirchhoff w/ impulsive peak loc')
hold on
scatter(dist_sorted(1:50:end),CNv_PeakLoc(sorted_by_dist(1:50:end)),'DisplayName','Vortex peak loc')
scatter(dist_sorted(1:50:end),CNtot_PeakLoc(sorted_by_dist(1:50:end)),'DisplayName','Total lift peak loc')
grid on 
legend('Location','NorthEast')
