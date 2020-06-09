function setLinFit(airfoil,varargin)

        % computes Talpha based on a vector of reduced pitch rates r and
        % corresponding dynamic stall angles alpha_ds. Talpha is the slope
        % of the curve fitting alpha_ds as a function of r.
        
        r0 = airfoil.r0;
        alpha_ss = airfoil.steady.alpha_ss; 
        
        % initialization of the vectors
        r = -ones(size(varargin));
        alpha_ds = -ones(size(varargin));
        
        for k=1:length(varargin)
            % call ramp the current ramp
            ramp = varargin{k};
            % Define reduced pitch rate if necessary ...
            if isempty(ramp.r)
                % compute it with alphadot
                ramp.setPitchRate(airfoil);
            end
            % ... and assign it
            r(k) = ramp.r;
            % Define experimental stall if necessary ...
            if isempty(ramp.i_CConset)
                ramp.findExpOnset();
            end
            % ... and assign it
            alpha_ds(k) = ramp.alpha_CConset;
        end
        
        %% Compute the polynomials for the linear fitting
        pr = polyfit(r(r>=r0),alpha_ds(r>=r0),1);
        pl = [(polyval(pr,r0)-alpha_ss)/r0 alpha_ss];
        Talpha = pi/180*pr(1);
        alpha_ds0 = pr(2);
        save(sprintf('../linfit_%s',airfoil.name),'Talpha','alpha_ds0','r0')
        
        %% Compute alpha_ds_lag and alpha_crit
        
        alpha_lag_ds = -ones(size(varargin));
        alpha_crit = -ones(size(varargin));
        % looking for the value of alpha_lag(alpha) at the point alpha_ds
        for k=1:length(varargin)
            ramp = varargin{k};
            ramp.computeAlphaLag(airfoil,Talpha);
            if isempty(ramp.alpha)
                alpha_lag_ds(k) = interp1(ramp.analpha,ramp.analpha_lag,alpha_ds(k));
            elseif isempty(ramp.i_continuous_grow)
                alpha_lag_ds(k) = interp1(ramp.alpha,ramp.alpha_lag,alpha_ds(k));
            else % if alpha_continuous_grow is defined
                alpha_lag_ds(k) = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),alpha_ds(k));
            end
        % define the critical alpha
        alpha_crit(k) = computeAlphaCrit(airfoil,alpha_ds0,ramp.r);
        end

        
        
        %% Plot the alpha_ds_r figure
        
        % plots the results for Sheng's model
        figure
        % blue dots
        plot(r,alpha_ds,'.','DisplayName','\alpha_{ds} (exp)','MarkerSize',20)
        hold on
        grid on
        xlabel('reduced pitch rate r (-)','FontSize',20);
        ylabel('\alpha_{ds} (°)','FontSize',20);
        ax = gca;
        axis([0 .06 0 30]);
        ax.FontSize = 20;
        % red line
        plot(sort([0 r r0]),[alpha_ss polyval(pl,r(r<r0)) polyval(pr,r0) polyval(pr,r(r>=r0))],'Color','r','DisplayName','Linear fitting','LineWidth',2)
        title(sprintf('%s ($T_{\\alpha} = %.2f$)',airfoil.name,Talpha),'interpreter','latex','FontSize',20)
        % yellow dots
        plot(r,alpha_lag_ds,'.','DisplayName','\alpha_{ds} (lagged)','MarkerSize',20)
        % purple line
        plot(sort([0 r r0]),sort([alpha_ss alpha_crit alpha_ds0]),'--','DisplayName','\alpha_{crit}','LineWidth',2);
        legend('Location','NorthWest','FontSize',20)        
end

