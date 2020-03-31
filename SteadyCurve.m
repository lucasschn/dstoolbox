classdef SteadyCurve < handle
    properties
        alpha
        alpha_rad
        CN
        CN0 = 0;
        S1
        S2
        alpha1
        f
        CNalpha
    end
    methods
        % constructor
        function obj = SteadyCurve(alphaSteady,CNSteady,alpha_static_stall)
            obj.alpha = alphaSteady;
            obj.alpha_rad = deg2rad(obj.alpha);
            obj.CN = CNSteady;
            obj.CNalpha = diff(obj.CN)./diff(obj.alpha_rad); % in 1/rad
            if nargin < 3
                obj.alpha1 = StallAngle(obj);
            else
                obj.alpha1 = alpha_static_stall;
            end
        end
        function alpha_static_stall = StallAngle(obj)
            % Static stall angle
            dalpha = diff(obj.alpha);
            ialpha1 = 1;
            while (ialpha1==1 || dalpha(ialpha1-1)<0.01)
                ialpha1 = find(obj.CNalpha(ialpha1:end)<0,1)+ ialpha1;
            end
            alpha_static_stall=obj.alpha(ialpha1);
        end
        function CNslope = Slope(obj)
            % CN slope for attached flow
            CNslopes = obj.CNalpha(obj.alpha<10);
            alphaslopes = obj.alpha(1:length(CNslopes)+1);
            CNslope = sum(diff(alphaslopes).*CNslopes)/sum(diff(alphaslopes)); % mean weighted by the distance between two successive alphas
            
        end
        function plot(obj)
            figure
            plot(obj.alpha,obj.CN)
            grid on
            xlabel('\alpha (°)')
            ylabel('C_N')
        end
        function fitKirchhoff(obj)
            stall_slope_minus = obj.CNalpha(find(obj.alpha<obj.alpha1,1,'last'));
            stall_slope_plus = obj.CNalpha(find(obj.alpha>obj.alpha1,1,'first'));
            S10 = 0.3*deg2rad(obj.alpha1)/(2*sqrt(0.7))*((1+sqrt(0.7))/2-2*stall_slope_minus/(Slope(obj)*(1+sqrt(0.7)))).^(-1);
            S20 = 0.66*deg2rad(obj.alpha1)/(2*sqrt(0.7))*((1+sqrt(0.7))/2+2*stall_slope_plus/(Slope(obj)*(1+sqrt(0.7)))).^(-1);
            opts = optimset('Display','iter-detailed');
            
            Kfunc = @(x,alpha) Kirchhoff(obj,obj.alpha,x);
            
            [fitparams,res,~,exitflag] = lsqcurvefit(Kfunc,[S10 S20],obj.alpha,obj.CN,[0 0],[2 2],opts);
            
            switch(exitflag)
                case 1
                    disp('lsqcurvefit converged to a solution.')
                    obj.S1 = fitparams(1);
                    obj.S2 = fitparams(2);
                    sprintf('Norm of the residual is %0.2e',res)
                case 2
                    disp('Change in X too small.')
                case 3
                    disp('Change in RESNORM too small.')
                    obj.S1 = fitparams(1);
                    obj.S2 = fitparams(2);
                    sprintf('Norm of the residual is %0.2e',res)
                case 4
                    disp('Computed search direction too small.')
                    warning('S1 and S2 have not been assigned, as lsqcurvefit has not conevrged to a solution')
                case 0
                    disp('Too many function evaluations or iterations.')
                    warning('S1 and S2 have not been assigned, as lsqcurvefit has not conevrged to a solution')
                case -1
                    disp('Stopped by output/plot function.')
                    warning('S1 and S2 have not been assigned, as lsqcurvefit has not conevrged to a solution')
                case -2
                    disp('Bounds are inconsistent.')
                    warning('S1 and S2 have not been assigned, as lsqcurvefit has not conevrged to a solution')
            end
        end
        function plotKirchhoff(obj)
            figure
            plot(obj.alpha,obj.CN,'DisplayName','data')
            hold on
            plot(obj.alpha,Kirchhoff(obj,obj.alpha),'DisplayName','Kirchhoff model')
            grid on
            legend('Location','Best')
            xlabel('\alpha (°)')
            ylabel('C_N')
        end
        function setCN0(obj,CN0)
            if isnan(CN0)
                error('CN= cannot be NaN.')
            else
                obj.CN0 = CN0;
            end
        end
    end
end