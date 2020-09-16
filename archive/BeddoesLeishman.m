close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')
addpath('../plot_dir/')
addpath('../src/model/')
addpath('../src/common/')
addpath('../src/lib/')
%% Static data

data = load('../naca0012');
data = data.naca0012;

naca0012 = Airfoil('naca0012',0.1);
M = 0.3;
a = 340.3;
nu = 1.48e-5;
Re = M*a*naca0012.c/nu;

% Look at OpenJetCorr if using simcos data
naca0012.steady = SteadyCurve(data.alpha_st,data.CN_st);
naca0012.steady.plotCN()
saveas(gcf,'../fig/static_naca0012','png')

%% Dynamic data

% pitching motion
mean_rad = deg2rad(data.mean);
amp_rad = deg2rad(data.amp);
pitching = PitchingMotion('alpha',data.alpha_xp,'CN',data.CN_xp,'M',M,'k',0.1);
pitching.setName();
pitching.setSinus(naca0012,mean_rad,amp_rad,2*pi/0.032,-pi/2);

% Resample CNsteady to fit experimental unsteady alphas
pitching.setCNsteady(naca0012.steady)

% model parameters
% naca0012.steady.fitKirchhoff();
naca0012.steady.computeSlope();
naca0012.steady.setAlpha0();
naca0012.steady.fitKirchhoff();
naca0012.steady.plotKirchhoff();
Tp = 1.7;
Tf = 3;
Tv = 6;
Tvl = 6;

for kp=1:length(Tp)
    for kf=1:length(Tf)
        for kv=1:length(Tv)
            params = sprintf('S1=%0.1f, S2=%0.1f, Tp=%0.1f, Tf=%0.1f, Tv=%0.1f, Tvl=%0.1f',naca0012.steady.S1,naca0012.steady.S2,Tp(kp),Tf(kf),Tv(kv),Tvl);
            disp(params);
            pitching.BeddoesLeishman(naca0012,Tp(kp),Tf(kf),Tv(kv),Tvl,'experimental');
%             pitching.plotAlphas()
%             pitching.plotSeparation(naca0012,'normal',0) % last argument is saving figure or not  
            fig2 = figure;
            if length(Tv)==1
                if length(Tf) ==1 && mod(length(Tp),3)==0
                    subplot(3,length(Tp)/3,kp);
                elseif length(Tf) ==1 && mod(length(Tp),2)==0
                    subplot(2,length(Tp)/2,kp);
                else
                    subplot(length(Tp),length(Tf),length(Tf)*(kp-1)+kf)
                end
            elseif length(Tf)==1
                subplot(length(Tp),length(Tv),length(Tv)*(kp-1)+kv)
            elseif length(Tp)==1
                subplot(length(Tf),length(Tv),length(Tv)*(kf-1)+kv)
            end
%             plot(naca0012.steady.alpha,pitching.CNk,'DisplayName','C_{N,steady}')
            hold on
            plot_dir(pitching.alpha(1:length(pitching.CNI)),pitching.CNI,'DisplayName','C_N^I');
            plot_dir(pitching.alpha(1:length(pitching.CNC)),pitching. CNC,'DisplayName','C_N^C');
            plot_dir(pitching.alpha(1:length(pitching.CNp)),pitching.CNp,'DisplayName','C_N^p');
            %             plot(PitchingMotion.alphaxp(1:length(CNprime)),CNprime,'DisplayName','C_N''')
            %plot_dir(pitching.alpha(1:length(pitching.CNf)),pitching.CNf,'DisplayName','C_N^f');
            %             plot(PitchingMotion.alphaxp(1:length(CNv)),CNv,'DisplayName','C_N^v')
            plot_dir(pitching.alpha(1:length(pitching.CN_LB)),pitching.CN_LB,'DisplayName','C_{N,LB}','LineWidth',2,'Color','r');
            plot_dir(pitching.alpha,pitching.CN,'DisplayName','C_{N,xp}','LineWidth',2,'Color','b');
            
            if kp==1 && kf==1 && kv==1
                legend('Location','NorthWest','FontSize',20)
            end
            title(params)
            grid on
            xlabel('\alpha (�)')
            ylabel('C_N (-)')
        end
    end
end
