close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% Static data

data = load('naca0012');
data = data.naca0012;

naca0012 = Airfoil('naca0012',0.457);
M = 0.3;
a = 340.3;
nu = 1.48e-5;
Re = M*a*naca0012.c/nu;

% Look at OpenJetCorr if using simcos data
naca0012.steady = SteadyCurve(data.alphaSteady,data.CNsteady);

% experimental separation point
fexp = (2*sqrt(naca0012.steady.CN./(2*pi*naca0012.steady.alpha_rad))-1).^2;

%% Dynamic data

% pitching motion
mean_rad = deg2rad(data.mean);
amp_rad = deg2rad(data.amp);
pitching = PitchingMotion('alpha',data.alphaxp,'CN',data.CNxp,'V',M*a,'k',0.1);
pitching.setSinus(naca0012,mean_rad,amp_rad);

% Resample CNsteady
pitching.setCNsteady(pitching,naca0012.steady)

% model parameters
% values for seppoint in deg
naca0012.steady.fitKirchhoff();
naca0012.steady.plotKirchhoff();
% values for seppoint in rad
Tp = 6;
Tf = 3;
Tv = 3.5;

for kp=1:length(Tp)
    for kf=1:length(Tf)
        for kv=1:length(Tv)
            params = sprintf('S1=%0.1f, S2=%0.1f, Tp=%0.1f, Tf=%0.1f, Tv=%0.1f',naca0012.steady.S1,naca0012.steady.S2,Tp(kp),Tf(kf),Tv(kv));
            disp(params);
            [CN_LB,CNk,CNI,CNC,CNp,CNprime,CNf,CNv,f,fp,fpp,alphaf,alphaE] = ComputeUnsteadyLift(pitching,naca0012,M,Tp(kp),Tf(kf),Tv(kv));
            
%             fig4 = figure;
%             semilogy(steady.alpha,fexp,'DisplayName','experimental')
%             hold on
%             semilogy(steady.alpha,f,'DisplayName','model')
%             legend show
%             xlabel('\alpha')
%             ylabel('f')
%             grid on
%             
%             fig5 = figure;
%             plot(steady.alpha,steady.CN,'x')
%             hold on
%             plot(steady.alpha,CNk)
%             plot(steady.alpha1*ones(2),ylim,'r--')
%             legend('static data','Kirchhoff')
%             xlabel('\alpha')
%             ylabel('C_N')
%             grid on
%             title('Lift attached flow')
%             
%             fig6 = figure;
%             plot(PitchingMotion.alphaxp,'DisplayName','\alpha_{xp}')
%             hold on
%             plot(alphaf,'DisplayName','\alpha_f')
%             plot(alphaE,'DisplayName','\alpha_E')
%             grid on
%             legend('Location','Best')
%             ylabel('\alpha (°)')
%             
            fig7 = figure;
            plot(naca0012.steady.alpha,f,'DisplayName','f')
            hold on
            plot(pitching.alpha(1:length(fp)),fp,'DisplayName','f''')
            plot(pitching.alpha(1:length(fpp)),fpp,'DisplayName','f''''')
            legend show
            %legend('f','f ''','f ''''')
            title(params)
            grid on
            xlabel('\alpha (°)')
            ylabel('separation point (x/c)')
%             
            fig8 = figure;
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
            plot(naca0012.steady.alpha,CNk,'DisplayName','C_{N,steady}')
            hold on
            %             plot(PitchingMotion.alphaxp(1:length(CNI)),CNI,'DisplayName','C_N^I')
            %             plot(PitchingMotion.alphaxp(1:length(CNC)),CNC,'DisplayName','C_N^C')
            %             plot(PitchingMotion.alphaxp(1:length(CNp)),CNp,'DisplayName','C_N^p')
            %             plot(PitchingMotion.alphaxp(1:length(CNprime)),CNprime,'DisplayName','C_N''')
            %             plot(PitchingMotion.alphaxp(1:length(CNf)),CNf,'DisplayName','C_N^f')
            %             plot(PitchingMotion.alphaxp(1:length(CNv)),CNv,'DisplayName','C_N^v')
            plot(pitching.alpha(1:length(CN_LB)),CN_LB,'DisplayName','C_{N,LB}','LineWidth',2)
            plot(pitching.alpha,pitching.CN,'DisplayName','C_{N,xp}','LineWidth',2)
            
            if kp==1 && kf==1 && kv==1
                legend('Location','Best')
            end
            title(params)
            grid on
            xlabel('\alpha (°)')
            ylabel('C_N (-)')
        end
    end
end

%% Save figures
% saveas(fig4,'fig/f_fit.png')
% % saveas(fig,'fig/alpha.png')
% saveas(fig5,'fig/CN_fit.png')
saveas(fig7,'fig/f_curves.png')
% saveas(fig8,'fig/CN_stall.png')
