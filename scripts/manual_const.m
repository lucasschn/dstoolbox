close all
clear all
clc

r = [6.5E-03	1.3E-02	2.0E-02	2.6E-02	3.3E-02	3.9E-02	4.6E-02	5.2E-02	5.9E-02	6.5E-02	7.2E-02	7.9E-02	8.5E-02	9.2E-02	9.8E-02	1.0E-01	1.1E-01];

Tp = [1.8	2.9	2.9	2	1.7	1.5	0	0	1.5	1	1	0.7	1.3	0.7	1.3	0.7	1.7];
Tf = [2	3	2.7	2	1.5	1.5	1.6	1.7	1.4	1.2	1.2	0.7	1.4	1	1.3	0.7	1.5];
Tv = [2.5	2.5	2.5	2.5	2.5	2.5	4	2.5	2.5	2.5	2.5	2.5	2.5	2.5	2.5	2.5	2.5];
Tvl = [2	1	0.9	0.9	2	2	2	2	2	2	2	2	2	2	2	2	2];

figure
plot(r,Tp,'d','DisplayName','T_p','MarkerFaceColor',[1 0 0])
hold on
plot(r,Tf,'d','DisplayName','T_f','MarkerFaceColor',[0 1 0])
plot(r,Tp+Tf,'d','DisplayName','T_p+T_f','MarkerFaceColor',[1 1 0])
grid minor 
xlabel('r')
ylabel('T')
legend('Location','NorthEast')