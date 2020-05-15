classdef Airfoil < handle
    properties
        steady % corresponding steady curve
        name
        c
        figExpfit
        Talpha
        t_ds
        t_ss
        r
        alpha_ds
        A
        B
        % Sheng 2008
        r0 = 0.04 % reduced pitch rate for bilinear transition  
        alpha_ds0
        pr
        pl
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
        function setTalpha(obj)
            obj.Talpha = load('Talpha_flatplate.mat');
        end
        function computeTalphaSheng(obj,r,alpha_ds)
            % computes Talpha based on a vector of reduced pitch rates r and
            % corresponding dynamic stall angles alpha_ds. Talpha is the slope
            % of the curve fitting alpha_ds as a function of r.
            obj.pr = polyfit(r(r>=obj.r0),alpha_ds(r>=obj.r0),1);
            obj.alpha_ds0 = obj.pr(2);
            obj.Talpha = pi/180*obj.pr(1); % deg; % rad
            obj.pl = [(polyval(obj.pr,obj.r0)-obj.steady.alpha_ss)/obj.r0 obj.steady.alpha_ss];
            save('linfit_flatplate.mat','-struct','obj','Talpha','alpha_ds0')
        end
        function setExpfit(obj,A,B)
            % method to set the fit from given data
            obj.A = A;
            obj.B = B;
        end
        function fitExpfit(obj)
            % In Sheng 2006, the curve fitting alpha_ds as a function of r is
            % linear. Here it is exponential.
            alpha_ss = obj.steady.alpha_ss;
            alpha_ds_r = @(x,r) x(1)-(x(1)-alpha_ss)*exp(-x(2)*r);
            % compute the exponential fit
            xopt = lsqcurvefit(alpha_ds_r,[obj.alpha_ds(end) 1],obj.r,obj.alpha_ds,[0 0],[Inf Inf]);
            obj.A = xopt(1);
            obj.B = xopt(2);
            % plot expfit
            figure(obj.figExpfit)
            subplot(311)
            plot(obj.r,obj.alpha_ds,'.','DisplayName','\alpha_{ds} (exp)','MarkerSize',20)
            grid on
            ylabel('\alpha_{ds} (�)','FontSize',20);
            ax = gca;
            ax.FontSize = 20;
            hold on
            plot(obj.r,alpha_ds_r(xopt,obj.r),'LineWidth',2,'DisplayName','exponential fit')
            subplot(313)
            plot(obj.r,obj.t_ds-obj.t_ss,'.','MarkerSize',20)
        end
        function computeTalphaExpfit(obj,varargin)
            % computes Talpha based on a vector of reduced pitch rates r and
            % corresponding dynamic stall angles alpha_ds. varargin is a
            % collection of ramps.
            
            % determination of Talpha so that alpha_lag(t_ds) = alpha_ss
            obj.Talpha = -ones(size(varargin));
            for k = 1:length(varargin)
                ramp = varargin{k};
                obj.Talpha(k) = obj.findTalpha(ramp);
                obj.r(k) = ramp.r;
            end
            
            Talpha = @(x,r) 5*(1-exp(-r/x(1)).*cos(x(2)*r)) % cos without phase so that Talpha(0)=0
            
            xopt = lsqcurvefit(Talpha,[.04 80],obj.r,obj.Talpha);
            
            % plot Talpha and its fit on the lower graph
            figure(obj.figExpfit)
            subplot(312)
            plot(obj.r,obj.Talpha,'.','MarkerSize',20,'DisplayName','T_\alpha')
            hold on
            plot(obj.r,Talpha(xopt,obj.r),'DisplayName','fit for T_\alpha')
            xlabel('reduced pitch rate r (-)','FontSize',20);
            ylabel('T_\alpha','FontSize',20)
            grid on
        end
        function [Talpha,t0] = findTalpha(obj,ramp)
            t0 = interp1(ramp.analpha,ramp.t,0);
            t_ds = ramp.t(ramp.i_CConset)-t0;
            K = ramp.alphadot;
            syms tau % dimensional time constant
            sol = solve(obj.steady.alpha_ss == K*(t_ds - tau*(1-exp(-t_ds/tau))),'Real',true,'IgnoreAnalyticConstraints',true);
            Talpha = 2*ramp.V/obj.c *  double(sol); % in adimensional time here
        end

        function alpha_lag_ds = computeAlphaLagDS(obj,ramp) 
            % looking for the value of alpha_lag(alpha) at the point alpha_ds

            % initialization of the vectors
            obj.r = -ones(size(varargin));
            obj.alpha_ds = -ones(size(varargin));
            alpha_lag_ds = -ones(size(varargin));

            for k=1:length(varargin)
                % call ramp the current ramp
                ramp = varargin{k};
                % Define reduced pitch rate if necessary ...
                if isempty(ramp.r)
                    % compute it with alphadot
                    ramp.setPitchRate(obj);
                end
                % ... and assign it
                obj.r(k) = ramp.r;
                % Define experimental stall if necessary ...
                if isempty(ramp.i_CConset)
                    ramp.findExpOnset();
                end
                % ... and assign it                
                obj.alpha_ds(k) = ramp.alpha_CConset;
            end

            % fill the vector with the value for each ramp (for both models)
            if isempty(ramp.alpha)
                alpha_lag_ds(k) = interp1(ramp.analpha,ramp.analpha_lag,alpha_ds(k));
            elseif isempty(ramp.i_continuous_grow)
                alpha_lag_ds(k) = interp1(ramp.alpha,ramp.alpha_lag,alpha_ds(k));
            else % if alpha_continuous_grow is defined
                alpha_lag_ds(k) = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),alpha_ds(k));
            end
        end

        function Sheng(obj,airfoil,varargin)
            %% extract r and alpha_ds from arguments
            % argument is a set of RampUpMotions
            
            % initialization of the vectors
            alpha_crit = -ones(size(varargin)); 
            % compute alpha_lag using Talpha and finds alpha_lagonset
            obj.computeTalphaSheng(r,alpha_ds)
            % definition of alpha_crit depending on r (only for Sheng)
            for k=1:length(varargin)
                ramp = varargin{k};
                ramp.findModelOnset(obj); % alpha_lagonset = alpha_lag_ds only if Talpha is correct                
                if ramp.r >= obj.r0
                    alpha_crit(k) = obj.alpha_ds0;
                else % if r<r0
                    alpha_crit(k) = obj.steady.alpha_ss + (obj.alpha_ds0 - obj.steady.alpha_ss)*ramp.r/obj.r0;
                end             
            end
            obj.plotSheng(r,alpha_ds,alpha_lag_ds,alpha_crit);          
        end
        
        function Expfit(obj,airfoil,varargin)

            % initialization of the vectors
            obj.t_ss = -ones(size(varargin));
            obj.t_ds = -ones(size(varargin));
            % compute alpha_lag using Talpha and finds alpha_lagonset
            obj.fitExpfit();
            obj.computeTalphaExpfit(varargin{:});
            alpha_lag_ds = -ones(size(varargin));
            for k=1:length(varargin)
                ramp = varargin{k};
                alpha_lag_ds = obj.computeAlphaLagDS(ramp)
                % assign the two last (they have to defined already) (only for Sheng)
                obj.t_ss(k) = ramp.t(find(ramp.alpha>airfoil.steady.alpha_ss,1));
                obj.t_ds(k) = ramp.t(ramp.i_CConset);
                ramp.computeAlphaLag(obj,interp1(obj.r,obj.Talpha,ramp.r));
                %ramp.findModelOnset(obj); % alpha_lagonset = alpha_lag_ds only if Talpha is correct
                % looking for the value of alpha_lag(alpha) at the point alpha_ds
                if isempty(ramp.alpha)
                    alpha_lag_ds(k) = interp1(ramp.analpha,ramp.analpha_lag,obj.alpha_ds(k));
                elseif isempty(ramp.i_continuous_grow)
                    alpha_lag_ds(k) = interp1(ramp.alpha,ramp.alpha_lag,obj.alpha_ds(k));
                else % if alpha_continuous_grow is defined
                    alpha_lag_ds(k) = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),obj.alpha_ds(k));
                end
            end
            obj.plotExpfit(alpha_lag_ds);             
        end

        function plotSheng(obj,r,alpha_ds,alpha_lag_ds,alpha_crit)
            % plots the results for Sheng's model
            obj.figSheng = figure;
            % blue dots
            plot(r,alpha_ds,'.','DisplayName','\alpha_{ds} (exp)','MarkerSize',20)
            hold on
            grid on
            xlabel('reduced pitch rate r (-)','FontSize',20);
            ylabel('\alpha_{ds} (�)','FontSize',20);
            ax = gca;
            axis([0 .06 0 30]);
            ax.FontSize = 20;
            % red line
            plot(sort([0 r obj.r0]),[obj.steady.alpha_ss polyval(obj.pl,r(r<obj.r0)) polyval(obj.pr,obj.r0) polyval(obj.pr,r(r>=obj.r0))],'Color','r','DisplayName','Linear fitting','LineWidth',2)
            title(sprintf('%s ($T_{\\alpha} = %.2f$)',obj.name,obj.Talpha),'interpreter','latex','FontSize',20)
            % yellow dots
            plot(r,alpha_lag_ds,'.','DisplayName','\alpha_{ds} (lagged)','MarkerSize',20)
            % purple line
            plot(sort([0 r obj.r0]),sort([obj.steady.alpha_ss alpha_crit obj.alpha_ds0]),'--','DisplayName','\alpha_{ds,0}','LineWidth',2);
            legend('Location','NorthWest','FontSize',20)
        end
        function plotExpfit(obj,alpha_lag_ds)
            figure(obj.figExpfit)
            subplot(311)
            hold on
            %             plot(obj.r,obj.D1.*obj.r+obj.alpha_ds0,'DisplayName','Linear fitting','LineWidth',2)
            %             title(sprintf('%s ($T_{\\alpha} = %.2f$)',obj.name,obj.Talpha),'interpreter','latex','FontSize',20)
            plot(obj.r,alpha_lag_ds,'.','DisplayName','\alpha_{ds} (lagged)','MarkerSize',20)
            plot(obj.r,ones(size(obj.r)).*obj.steady.alpha_ss,'--','DisplayName','\alpha_{ss}','LineWidth',2);
            legend('FontSize',20,'Location','East')
        end     
    end
end
