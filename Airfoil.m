classdef Airfoil < handle
    properties
        name
        c
        alpha_ds0
        pr
        pl
        Talpha
        r0 = 0.03 % reduced pitch rate for bilinear transition  
        steady % corresponding steady curve
    end
    methods
        % Unique constructor with airfoil's name and chord length. Airfoil
        % user-defined name property can be distinct from the instance
        % name. This allows characters in the name that are not allowed for
        % variables. The chord length has to be positive, otherwise an
        % error is thrown.
        function obj = Airfoil(name,c)
            obj.name = name;
            if c > 0
                obj.c = c;
            else
                error('Please enter a valid chord length.')
            end
        end
        function b = b(obj)
            b = obj.c/2;
        end
        function setTalpha()
            obj.Talpha = load('Talpha_flatplate.mat');
        end
        function computeTalpha(obj,r,alpha_ds)
            % computes Talpha based on a vector of reduced pitch rates r and
            % corresponding dynamic stall angles alpha_ds. Talpha is the slope
            % of the curve fitting alpha_ds as a function of r.
            obj.pr = polyfit(r(r>=obj.r0),alpha_ds(r>=obj.r0),1);
            D1 = obj.pr(1); % deg
            obj.alpha_ds0 = obj.pr(2);
            obj.Talpha = pi/180*D1; % rad
            obj.pl = polyfit(r(r<obj.r0),alpha_ds(r<obj.r0),1);
        end
        function fig = Sheng(obj,varargin)
            %% extract r and alpha_ds from arguments
            % argument is a set of RampUpMotions
            r = -ones(size(varargin));
            alpha_ds = -ones(size(varargin));
            for k=1:nargin-1 % first argument is self
                ramp = varargin{k};
                if isempty(ramp.r)
                    % compute it with alphadot
                    ramp.setPitchRate(obj);
                end
                r(k) = ramp.r;
                % Define experimental stall if necessary
                if isempty(ramp.i_CConset)
                    ramp.findExpOnset();
                end
                alpha_ds(k) = ramp.alpha_CConset;         
            end
            
            % compute alpha_lag using Talpha and finds alpha_lagonset
            obj.computeTalpha(r,alpha_ds)
            alpha_lag_ds = -ones(size(varargin));
            alpha_crit = -ones(size(varargin));
            for k=1:nargin-1
                ramp = varargin{k};
                ramp.findModelOnset(obj); % alpha_lagonset = alpha_lag_ds only if Talpha is correct
                % looking for the value of alpha_lag(alpha) at the point alpha_ds
                if isempty(ramp.alpha)
                    alpha_lag_ds(k) = interp1(ramp.analpha,ramp.analpha_lag,alpha_ds(k));
                elseif isempty(ramp.i_continuous_grow)
                    alpha_lag_ds(k) = interp1(ramp.alpha,ramp.alpha_lag,alpha_ds(k));
                else % if alpha_continuous_grow is defined
                    alpha_lag_ds(k) = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),alpha_ds(k));
                end
                if ramp.r >= obj.r0
                    alpha_crit(k) = obj.alpha_ds0;
                else
                    alpha_crit(k) = obj.steady.alpha_ss + (obj.alpha_ds0 - obj.steady.alpha_ss)*ramp.r/obj.r0;
                end
            end
            
            fig = obj.plotSheng(r,alpha_ds,alpha_lag_ds,alpha_crit);
            
        end
        function fig = plotSheng(obj,r,alpha_ds,alpha_lag_ds,alpha_crit)
            fig = obj.plotDS(r,alpha_ds);
            figure(fig)
            plot(r,[polyval(obj.pl,r(r<obj.r0)) polyval(obj.pr,r(r>=obj.r0))],'Color','r','DisplayName','Linear fitting','LineWidth',2)
            title(sprintf('%s ($T_{\\alpha} = %.2f$)',obj.name,obj.Talpha),'interpreter','latex','FontSize',20)
            plot(r,alpha_lag_ds,'.','DisplayName','\alpha_{ds} (lagged)','MarkerSize',20)
            plot(r,alpha_crit,'--','DisplayName','\alpha_{ds,0}','LineWidth',2);
            legend('Location','NorthWest','FontSize',20)
        end
    end
    methods (Static)
       function fig = plotDS(r,alpha_ds)
            fig = figure;
            plot(r,alpha_ds,'.','DisplayName','\alpha_{ds} (exp)','MarkerSize',20)
            hold on
            grid on
            xlabel('reduced pitch rate r (-)','FontSize',20);
            ylabel('\alpha_{ds} (°)','FontSize',20);
            ax = gca;
            ax.FontSize = 20;
       end
    end
end
