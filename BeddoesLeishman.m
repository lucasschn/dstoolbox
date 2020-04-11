close all
clear all
clc
set(0,'DefaultFigureWindowStyle','docked')

%% Static data

data = load('naca0012');
data = data.naca0012;

naca0012 = Airfoil('naca0012',0.1);
M = 0.3;
a = 340.3;
nu = 1.48e-5;
Re = M*a*naca0012.c/nu;

% Look at OpenJetCorr if using simcos data
naca0012.steady = SteadyCurve(data.alphaSteady,data.CNsteady);

%% Dynamic data

% pitching motion
mean_rad = deg2rad(data.mean);
amp_rad = deg2rad(data.amp);
pitching = PitchingMotion('alpha',data.alphaxp,'CN',data.CNxp,'V',M*a,'k',0.1);
pitching.setName();
pitching.setSinus(naca0012,mean_rad,amp_rad);

% Resample CNsteady to fit experimental unsteady alphas
pitching.setCNsteady(naca0012.steady)

% model parameters
naca0012.steady.fitKirchhoff();
naca0012.steady.plotKirchhoff();
Tp = 3.7;
Tf = 3;
Tv = 4;

for kp=1:length(Tp)
    for kf=1:length(Tf)
        for kv=1:length(Tv)
            params = sprintf('S1=%0.1f, S2=%0.1f, Tp=%0.1f, Tf=%0.1f, Tv=%0.1f',naca0012.steady.S1,naca0012.steady.S2,Tp(kp),Tf(kf),Tv(kv));
            disp(params);
            computeUnsteadyLift(pitching,naca0012,Tp(kp),Tf(kf),Tv(kv));
            pitching.plotAlphas()
            pitching.plotSeparation(naca0012,'normal',0) % last argument is saving figure or not  
            
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
            plot(naca0012.steady.alpha,pitching.CNk,'DisplayName','C_{N,steady}')
            hold on
            %             plot(PitchingMotion.alphaxp(1:length(CNI)),CNI,'DisplayName','C_N^I')
            %             plot(PitchingMotion.alphaxp(1:length(CNC)),CNC,'DisplayName','C_N^C')
            %             plot(PitchingMotion.alphaxp(1:length(CNp)),CNp,'DisplayName','C_N^p')
            %             plot(PitchingMotion.alphaxp(1:length(CNprime)),CNprime,'DisplayName','C_N''')
            %             plot(PitchingMotion.alphaxp(1:length(CNf)),CNf,'DisplayName','C_N^f')
            %             plot(PitchingMotion.alphaxp(1:length(CNv)),CNv,'DisplayName','C_N^v')
            plot(pitching.alpha(1:length(pitching.CN_LB)),pitching.CN_LB,'DisplayName','C_{N,LB}','LineWidth',2)
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
