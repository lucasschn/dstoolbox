classdef RampUpMotion < handle
    properties
        name
        alpha
        alpha_continuous_grow
        i_continuous_grow
        alpha_rad
        analpha
        alpha_lag
        analpha_lag
        alpha_CConset
        i_CConset
        i_CNonset
        alpha_CLonset
        alpha_lagonset % corresponds to alpha'_ds
        alphadot %°/s
        CL
        CD
        CN
        CNslope1
        CNslope2
        CC
        Ts
        t
        V
        r % reduced pitch rate
        rt % instantaneous red. pitch rate (shoudl be const)
    end
    methods
        % constructor (not possible to overload in Matlab)
        function obj = RampUpMotion(varargin)
            
            obj.name = varargin{1};
            
            if nargin >= 2 %(name,alpha)
                
                if length(varargin{2}) > 1
                    obj.alpha = varargin{2};
                    obj.alpha_rad = deg2rad(obj.alpha);
                    
                elseif nargin == 3 %(name,alpha,Ts)
                    obj.Ts = varargin{3};
                    obj.Ts = varargin{3};
                    
                elseif nargin == 4 %(name,alpha,CN,Ts)
                    obj.Ts = varargin{2};
                    obj.CN = varargin{3};
                    obj.Ts = varargin{4};
                    obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
                    
                elseif nargin == 5 %(name,alpha,CN,V,Ts)
                    obj.Ts = varargin{2};
                    obj.CN = varargin{3};
                    obj.V = varargin{4};
                    obj.Ts = varargin{5};
                    obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
                else
                    error('The second argument must be the measured incidence angles vector.')
                end
            end
        end
%         function name = objName(obj)
%             name = inputname(1); % 1 stands for arg nb 1
%         end
        function isolateRamp(obj)
            dalpha = diff(obj.alpha);
            i_end = find((abs(dalpha)<1e-2) .* (obj.alpha(1:end-1)>10),ceil(1/obj.Ts));
            i_end = i_end(end);
            t_end = obj.t(i_end);
            fprintf('Data will be cutoff at %.2fs \n',t_end)
            obj.alpha = obj.alpha(1:i_end);
            obj.t = obj.t(1:i_end);
            obj.CN = obj.CN(1:i_end);
            obj.CC = obj.CC(1:i_end);
            obj.CL = obj.CL(1:i_end);
            obj.CD = obj.CD(1:i_end);
            % then find valid values for continuously increasing alphas
            i_grow = findGrowingIndices(obj.alpha);
            for l=1:length(i_grow)
                ls(l) = length(i_grow{l});
            end
            [~,imax] = max(ls);
            obj.i_continuous_grow = i_grow{imax};
            obj.alpha_continuous_grow = obj.alpha(i_grow{imax});
        end
        function sett(obj,t)
            if length(t)==length(obj.alpha)
                obj.t = t;
                obj.Ts = mean(diff(t));
            else
                error('t and alpha must be of same length.')
            end
        end
        function setCL(obj,CL)
            if length(CL)==length(obj.alpha)
                obj.CL = CL;
            else
                error('CL and alpha must be of same length.')
            end
        end
        function setCD(obj,CD)
            if length(CD)==length(obj.alpha)
                obj.CD = CD;
            else
                error('CD and alpha must be of same length.')
            end
        end
        function setCN(obj,CN)
            if length(CN)==length(obj.alpha)
                obj.CN = CN;
            else
                error('CN and alpha must be of same length.')
            end
        end
        function setalphadot(obj,alphadot)
            dalphadt = diff(obj.alpha)./diff(obj.t);
            if nargin > 1
                obj.alphadot = alphadot;
            else
                obj.alphadot = max(dalphadt);
            end
            if ~isempty(obj.alpha)
                analpha0 = lsqcurvefit(@(x,xdata) obj.alphadot*xdata+x,0,obj.t(dalphadt>=5),obj.alpha(dalphadt>=5));
            else
                analpha0 = 0;
            end
            if isempty(obj.t)
                obj.t = 1:0.1:300;
            end
            obj.analpha = obj.alphadot*obj.t + analpha0;
        end
        function computeAirfoilFrame(obj)
            if ~isempty(obj.CD) && ~isempty(obj.CL)
                if ~isempty(obj.CC) || ~isempty(obj.CN)
                    warning('The data for CN and CC will be erased.')
                end
                obj.CN = obj.CL.*cos(obj.alpha_rad) + obj.CD.*sin(obj.alpha_rad);
                obj.CC = obj.CD.*cos(obj.alpha_rad) - obj.CL.*sin(obj.alpha_rad);
            else
                error('CL and CD must be defined to compute CN and CC.')
            end
        end
        function computeFlowFrame(obj)
            if ~isempty(obj.CC) && ~isempty(obj.CN)
                if ~isempty(obj.CD) || ~isempty(obj.CL)
                    warning('The data for CL and CD will be erased.')
                end
                obj.CL = obj.CN.*cos(obj.alpha_rad) - obj.CC.*sin(obj.alpha_rad);
                obj.CD = obj.CC.*cos(obj.alpha_rad) + obj.CN.*sin(obj.alpha_rad);
            else
                error('CL and CD must be defined to compute CN and CC.')
            end
        end
        function setPitchRate(obj,airfoil)
            obj.r = deg2rad(obj.alphadot)*airfoil.c/(2*obj.V);
            dalphadt = diff(obj.alpha)./diff(obj.t);
            obj.rt = deg2rad(dalphadt)*airfoil.c/(2*obj.V);
        end
        function findExpOnset(obj)
            % finds experimental dynamic stall onset for a specific ramp-up
            % experiment
            obj.CNslope1 = lsqcurvefit(@(x,xdata) x*xdata,0,obj.alpha(obj.alpha<15),obj.CN(obj.alpha<15));
            err = abs(obj.CN - obj.CNslope1*obj.alpha);
            i_CNonset = find(err>1e-2);
            [CNmax,iCNmax] = max(obj.CN);
            %obj.CNslope2 = lsqcurvefit(@(x,xdata) x*xdata+CNmax-x*xdata(iCNmax),0,obj.alpha(1:iCNmax),obj.CN(1:iCNmax));
            %             dCLdalpha = diff(obj.CL(obj.i_continuous_grow))./diff(obj.alpha_continuous_grow);
%             figure
%             plot(obj.alpha(1:length(dCLdalpha)),dCLdalpha)
%             hold on
%             xlabel('\alpha')
%             ylabel('dC_N/d\alpha')
%             grid on
            % onset based on CC curve
            [~,obj.i_CConset] = min(obj.CC);
            obj.alpha_CConset = obj.alpha(obj.i_CConset); % Sheng uses an inverted definition of CC
        end
        function computeAlphaLag(obj,airfoil)
            % computes the vector alpha_lag for an experimental time evolution of alpha
            s = obj.t*obj.V/(airfoil.c/2);
            if ~isempty(obj.alpha) % compute alpha_lag from alpha using Eq.5
                dalpha = diff(obj.alpha);
                obj.alpha_lag = zeros(size(obj.alpha));
                for k = 1:length(dalpha)
                    obj.alpha_lag(k+1) = obj.alpha_lag(k) + dalpha(k)*(1-exp(-s(k)/airfoil.Talpha));
                end
                if ~isempty(obj.alpha_lag) % compute analpha_lag from alpha_lag 
                    dalpha_lagdt = diff(obj.alpha_lag)./diff(obj.t);
                    alphadot_lag = max(dalpha_lagdt);
                    analpha_lag0 = lsqcurvefit(@(x,xdata) alphadot_lag*xdata+x,0,obj.t(dalpha_lagdt>=5),obj.alpha_lag(dalpha_lagdt>=5));
                    obj.analpha_lag =  alphadot_lag*obj.t + analpha_lag0;
                end
            elseif ~isempty(obj.analpha) % compute analpha_lag from analpha
                dalpha = diff(obj.analpha);
                obj.analpha_lag = zeros(size(obj.analpha));
                for k = 1:length(dalpha)
                    obj.analpha_lag(k+1) = obj.analpha_lag(k) + dalpha(k)*(1-exp(-s(k)/airfoil.Talpha));
                end
            end
        end
        function findModelOnset(obj,airfoil)
            % finds the Sheng-predicted dynamic stall angle for a specific
            % time-evolution of alpha
            obj.computeAlphaLag(airfoil);
            i_lagonset = find(obj.analpha_lag>airfoil.alpha_ds0,1);
            if ~isempty(obj.alpha)
                obj.alpha_lagonset = obj.alpha(i_lagonset);
            else 
                obj.alpha_lagonset = obj.analpha(i_lagonset);
            end
        end
        function plotCL(obj)
            figure
            plot(obj.alpha,obj.CL)
            grid on
            xlabel('\alpha (°)')
            ylabel('C_L (-)')
            title(obj.name)
        end
        function plotCD(obj)
            figure
            plot(obj.alpha,obj.CD)
            grid on
            xlabel('\alpha (°)')
            ylabel('C_D (-)')
            title(obj.name)
        end
        function plotCN(obj)
            figure
            plot(obj.alpha,obj.CN,'DisplayName','exp')
            hold on 
            plot(obj.alpha,obj.alpha*obj.CNslope1,'--','DisplayName','1st part fit')
            %plot(obj.alpha,obj.alpha*obj.CNslope2 + obj.CNslope2,'--','DisplayName','2nd part fit')
            grid on
            legend('Location','SouthEast')
            xlabel('\alpha (°)')
            ylabel('C_N (-)')
            title(obj.name)
        end
        function fig = plotCC(obj)
            fig = figure;
            plot(obj.alpha,obj.CC)
            hold on
            plot(obj.alpha_CConset,min(obj.CC),'rx')
            grid on
            xlabel('\alpha (°)')
            ylabel('C_C (-)')
            title(obj.name)
        end
        function plotAlpha(obj)
            figure
            plot(obj.t,obj.alpha,'DisplayName','exp')
            hold on
            if ~isempty(obj.analpha)
                
                plot(obj.t,obj.analpha,'--','DisplayName','ideal')
            end
            if ~isempty(obj.alpha_CConset)
                plot(obj.t(obj.i_CConset),obj.alpha_CConset,'rx','DisplayName','\alpha_{ds,CC}')
            end
            legend show
            grid on
            ylabel('\alpha (°)')
            xlabel('t (s)')
        end
        function plotAlphaLag(obj)
            figure
            if ~isempty(obj.alpha)
                plot(obj.t,obj.alpha,'DisplayName','\alpha')
                hold on
                plot(obj.t,obj.alpha_lag,'DisplayName','\alpha''')
            end
            if ~isempty(obj.analpha)
                plot(obj.t,obj.analpha,'--','DisplayName','ideal \alpha')
                hold on
                plot(obj.t,obj.analpha_lag,'--','DisplayName','ideal \alpha''')
            end
            if ~isempty(obj.alpha_CConset)
                plot(obj.t(obj.i_CConset),obj.alpha_CConset,'rx','DisplayName','\alpha_{ds,CC}')
                hold on
                plot(obj.t(obj.i_CConset),obj.alpha_lag(obj.i_CConset),'bx','DisplayName','\alpha''_{ds,CC}')
            end
            grid on
            xlabel('t (s)')
            ylabel('\alpha')
            legend('Location','SouthEast')
        end
        function plotPitchRate(obj)
            figure
            plot(obj.t(1:length(obj.rt)),obj.rt,'DisplayName','r(t)')
            hold on
            plot(obj.t,obj.r*ones(size(obj.t)),'--','DisplayName','r')
            plot(obj.t(obj.i_CConset),obj.rt(obj.i_CConset),'rx')
            legend show
            xlabel('t (s)')
            ylabel('r (-)')
            grid on
        end
    end
end