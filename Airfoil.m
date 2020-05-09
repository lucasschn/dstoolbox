classdef Airfoil < handle
    properties
        name
        c
        D1
        r
        alpha_ds
        Talpha
        r0 = 0.01 % reduced pitch rate at which linear fitting begins
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
        function fig = computeTalpha(obj,varargin)
            % computes Talpha based on a vector of reduced pitch rates r and
            % corresponding dynamic stall angles alpha_ds. Talpha is the slope
            % of the curve fitting alpha_ds as a function of r.
            alpha_ss = obj.steady.alpha_static_stall;
            alpha_ds_r = @(x,r) x(1)-(x(1)-alpha_ss)*exp(-x(2)*r);
            
            xopt = lsqcurvefit(alpha_ds_r,[obj.alpha_ds(end) 1],obj.r,obj.alpha_ds,[0 0 0 ],[Inf Inf Inf]);
            obj.D1 = 100*xopt(1) + 100*(xopt(1)-alpha_ss)*exp(-xopt(2)*obj.r); % deg
            %             obj.alpha_ds0 = p(2);
            %obj.Talpha = pi/180*obj.D1; % rad
            obj.Talpha = 2*0.5/obj.c*[obj.findTalpha(varargin{1}),obj.findTalpha(varargin{2}), obj.findTalpha(varargin{3}), obj.findTalpha(varargin{4}), obj.findTalpha(varargin{5}), obj.findTalpha(varargin{6}), obj.findTalpha(varargin{7}), obj.findTalpha(varargin{8})]; 
            fig = obj.plotDS();
            figure(fig)
            subplot(211)
            hold on
            plot(obj.r,alpha_ds_r(xopt,obj.r),'LineWidth',2,'DisplayName','exponential fit')
            subplot(212)
            plot(obj.r,obj.Talpha,'LineWidth',2,'DisplayName','T_\alpha')
            xlabel('reduced pitch rate r (-)','FontSize',20);
            ylabel('T_\alpha','FontSize',20)
            grid on
        end
        function [tau_sol,t0] = findTalpha(obj,ramp)
            t0 = interp1(ramp.analpha,ramp.t,0);
            t_ds = ramp.t(ramp.i_CConset)-t0;
            alpha_ss = obj.steady.alpha_static_stall;
            K = ramp.alphadot;
            syms tau
            
            sol = solve(alpha_ss == K*(t_ds - tau*(1-exp(-t_ds/tau))),'Real',true,'IgnoreAnalyticConstraints',true);
            
            tau_sol = double(sol); % in dimensional time here
        end
        function fig = Sheng(obj,varargin)
            %% extract r and alpha_ds from arguments
            % argument is a set of RampUpMotions
            obj.r = -ones(size(varargin));
            obj.alpha_ds = -ones(size(varargin));
            for k=1:nargin-1 % first argument is self
                ramp = varargin{k};
                if isempty(ramp.r)
                    % compute it with alphadot
                    ramp.setPitchRate(obj);
                end
                if ramp.r>=obj.r0
                    obj.r(k) = ramp.r;
                    % Define experimental stall if necessary
                    if isempty(ramp.i_CConset)
                        ramp.findExpOnset();
                    end
                    obj.alpha_ds(k) = ramp.alpha_CConset;
                end
                
            end
            
            % compute alpha_lag using Talpha and finds alpha_lagonset
            fig = obj.computeTalpha(varargin{:});
            alpha_lag_ds = -ones(size(varargin));
            for k=1:nargin-1
                ramp = varargin{k};
                ramp.computeAlphaLag(obj);
                %ramp.findModelOnset(obj); % alpha_lagonset = alpha_lag_ds only if Talpha is correct
                % looking for the value of alpha_lag(alpha) at the point alpha_ds
                if ramp.r >= obj.r0
                    if isempty(ramp.alpha)
                        alpha_lag_ds(k) = interp1(ramp.analpha,ramp.analpha_lag,obj.alpha_ds(k));
                    elseif isempty(ramp.i_continuous_grow)
                        alpha_lag_ds(k) = interp1(ramp.alpha,ramp.alpha_lag,obj.alpha_ds(k));
                    else % if alpha_continuous_grow is defined
                        alpha_lag_ds(k) = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),obj.alpha_ds(k));
                    end
                end
            end
            
            fig = obj.plotSheng(fig,alpha_lag_ds);
            
        end
        function fig = plotSheng(obj,fig,alpha_lag_ds)
            figure(fig)
            subplot(211)
            hold on
%             plot(obj.r,obj.D1.*obj.r+obj.alpha_ds0,'DisplayName','Linear fitting','LineWidth',2)
%             title(sprintf('%s ($T_{\\alpha} = %.2f$)',obj.name,obj.Talpha),'interpreter','latex','FontSize',20)
            plot(obj.r,alpha_lag_ds,'.','DisplayName','\alpha_{ds} (lagged)','MarkerSize',20)
            plot(obj.r,ones(size(obj.r)).*obj.steady.alpha_static_stall,'--','DisplayName','\alpha_{ss}','LineWidth',2);
            legend('FontSize',20,'Location','East')
        end
        function fig = plotDS(obj)
            fig = figure;
            subplot(211)
            plot(obj.r,obj.alpha_ds,'.','DisplayName','\alpha_{ds} (exp)','MarkerSize',20)
            grid on
            ylabel('\alpha_{ds} (°)','FontSize',20);
            ax = gca;
            ax.FontSize = 20;
        end
    end
    end
