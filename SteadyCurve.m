classdef SteadyCurve < handle
    properties
        alpha
        alpha_rad
        CN
        CN0 = 0;
        alpha_static_stall % denoted alpha_1 in Beddoes-Leishman
        f
        CNalpha % 1/deg
        slope % pre-stall slope in 1/deg
        slope_rad % pre-stall slope in 1/rad
        fexp
        % Kirchhoff 
        S1
        S2
    end
    methods
        % constructor
        function obj = SteadyCurve(alphaSteady,CNSteady,alpha_static_stall)
            obj.alpha = alphaSteady;
            obj.alpha_rad = deg2rad(obj.alpha);
            obj.CN = CNSteady;
            obj.CNalpha = diff(obj.CN)./diff(obj.alpha); % in 1/deg
            if nargin < 3
                obj.alpha_static_stall = StallAngle(obj);
            else
                obj.alpha_static_stall = alpha_static_stall;
            end
        end
        function alpha_static_stall = StallAngle(obj)
            % Static stall angle
            dalpha = diff(obj.alpha);
            ialphass = 1;
            while (ialphass==1 || dalpha(ialphass-1)<0.01)
                ialphass = find(obj.CNalpha(ialphass:end)<0,1)+ ialphass;
            end
            alpha_static_stall=obj.alpha(ialphass);
        end
        function computeSlope(obj)
            % CN slope for attached flow. Should be around 2pi/beta
            % converted to degrees. 
            CNslopes = obj.CNalpha(obj.alpha<10); % 1/deg
            if isempty(CNslopes)
                obj.slope_rad = 2*pi; %1/rad
                obj.slope = obj.slope_rad*pi/180; %1/deg
            else
                alphaslopes = obj.alpha(1:length(CNslopes)+1); % deg
                obj.slope = sum(diff(alphaslopes).*CNslopes)/sum(diff(alphaslopes)); % mean weighted by the distance between two successive alphas
                obj.slope_rad = obj.slope*180/pi;
            end
        end
        function plot(obj)
            figure
            plot(obj.alpha,obj.CN)
            grid on
            xlabel('\alpha (°)')
            ylabel('C_N')
        end
        function fitKirchhoff(obj)
            obj.computeSlope();
            stall_slope_minus = obj.CNalpha(find(obj.alpha<obj.alpha_static_stall,1,'last'));
            stall_slope_plus = obj.CNalpha(find(obj.alpha>obj.alpha_static_stall,1,'first'));
            S10 = 0.3*deg2rad(obj.alpha_static_stall)/(2*sqrt(0.7))*(((1+sqrt(0.7))/2)^2-stall_slope_minus/obj.slope_rad).^(-1);
            S20 = 0.66*deg2rad(obj.alpha_static_stall)/(2*sqrt(0.7))*(((1+sqrt(0.7))/2)^2-stall_slope_plus/obj.slope_rad).^(-1);
            opts = optimset('Display','off');
            
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
                    warning('S1 and S2 have not been assigned, as lsqcurvefit has not converged to a solution')
                case 0
                    disp('Too many function evaluations or iterations.')
                    warning('S1 and S2 have not been assigned, as lsqcurvefit has not converged to a solution')
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
            if nargin == 2
                if isnan(CN0)
                    error('CN= cannot be NaN.')
                else
                    obj.CN0 = CN0;
                end
            elseif nargin == 1
                obj.CN0 = interp1(obj.alpha,obj.CN,0);
            end
                
        end
        function computeSeparation(obj)
            % computes the experimental separation point using Kirchhoff
            % model
            obj.fexp = max([zeros(size(obj.CN)),min([ones(size(obj.CN)), (2*sqrt((obj.CN-obj.CN0)./(obj.slope*obj.alpha))-1).^2],[],2)],[],2);
            
        end
    end
end