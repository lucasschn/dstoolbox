classdef RampUpMotion < AirfoilMotion
    properties
        alpha_continuous_grow
        i_continuous_grow
        % experimental pitch angle evolutions
        alpha_lag
        analpha_lag
        % experimental dynamic stall angles
        alpha_CConset
        alpha_CLonset
        % their indices
        i_CConset
        i_CNonset
        % model outputs
        alpha_lagonset % corresponds to alpha'_ds
        alphadot %°/s
        alpha_onset % modelled one
        % Fitting parameters
        CNslope1
        CNslope2
        % experimental parameters
        r % reduced pitch rate
        f_pts

    end
    methods
        % convenient constructor with name/value pair of any attribute of RampUpMotion
        function obj = RampUpMotion(varargin)
            obj@AirfoilMotion(varargin)
            p = inputParser;
            % Add name / default value pairs
            prop = properties('RampUpMotion'); % makes a cell array of all properties of the specified ClassName
            for k=1:length(prop)
                p.addParameter(prop{k},[]);
            end
            p.parse(varargin{:}); % {:} is added to take the content of the cells
            for k=1:length(prop)
                eval(sprintf('obj.%s = p.Results.%s;',prop{k},prop{k}));
            end
            if ~isempty(obj.alpha)
                obj.alpha_rad = deg2rad(obj.alpha);
            end
            if isempty(obj.Ts) && ~isempty(obj.t)
                obj.Ts = mean(diff(obj.t));
            end
        end
        function isolateRamp(obj)
            dalpha = diff(obj.alpha);
            i_end = find((abs(dalpha/obj.Ts)<1e-2) .* (obj.alpha(1:end-1)>10),ceil(1/obj.Ts));
            i_end = i_end(end);
            t_end = obj.t(i_end);
            fprintf('Data will be cutoff at %.2fs \n',t_end)
            if 1
            obj.alpha = obj.alpha(1:i_end);
            obj.t = obj.t(1:i_end);
            obj.CN = obj.CN(1:i_end);
            obj.CC = obj.CC(1:i_end);
            obj.CL = obj.CL(1:i_end);
            obj.CD = obj.CD(1:i_end);
            end
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
        function setAlphaDot(obj,alphadot)
            % sets alphadot in degrees
            dalphadt = diff(obj.alpha)./diff(obj.t);
            if nargin > 1
                obj.alphadot = alphadot; % deg/s
            else
                obj.alphadot = max(dalphadt); % deg/s
            end
            if ~isempty(obj.alpha)
                analpha0 = lsqcurvefit(@(x,xdata) obj.alphadot*xdata+x,0,obj.t(dalphadt>=5),obj.alpha(dalphadt>=5));
            else
                analpha0 = 0;
            end
            if isempty(obj.t)
                obj.t = linspace(0,round(35/obj.alphadot,1),500);
            end
            obj.analpha = obj.alphadot*obj.t + analpha0;
        end
        function setPitchRate(obj,airfoil)
            obj.r = deg2rad(obj.alphadot)*airfoil.c/(2*obj.V);
            dalphadt = diff(obj.alpha)./diff(obj.t);
            obj.rt = deg2rad(dalphadt)*airfoil.c/(2*obj.V);
        end
        function findExpOnset(obj)
            % finds experimental dynamic stall onset for a specific ramp-up
            % experiment
%             obj.CNslope1 = lsqcurvefit(@(x,xdata) x*xdata,0,obj.alpha(obj.alpha<15),obj.CN(obj.alpha<15));
%             err = abs(obj.CN - obj.CNslope1*obj.alpha);
%             i_CNonset = find(err>1e-2);
%             [CNmax,iCNmax] = max(obj.CN);
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
                obj.analpha_lag(1) = obj.analpha(1); % = analpha0, as obj.t(1)=0
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
            if isempty(i_lagonset)
                warning('The airfoil "%s" does not show stall in the experiment %s.',airfoil.name,obj.name)
            else
                if ~isempty(obj.alpha)
                obj.alpha_lagonset = obj.alpha(i_lagonset);
                obj.alpha_onset = interp1(obj.alpha_lag(obj.i_continuous_grow),obj.alpha_continuous_grow,obj.alpha_lagonset);
                elseif ~isempty(obj.analpha)
                obj.alpha_lagonset = obj.analpha(i_lagonset);
                [c,ia] = unique(obj.analpha_lag);
                obj.alpha_onset = interp1(c,obj.analpha(ia),obj.alpha_lagonset);
                else
                error('Impossible to define stall angle. %s has no angle of attack defined.',obj.name)
                end
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
        function plotCNLag(obj)
            figure
            if ~isempty(obj.alpha_lag)
                plot(obj.alpha,obj.CN,'DisplayName','C_N(\alpha)')
                hold on
                plot(obj.alpha_lag,obj.CN,'--','DisplayName','C_N(\alpha'')')
            elseif ~isempty(obj.analpha_lag)
                plot(obj.analpha,obj.CN,'DisplayName','C_N(\alpha)')
                hold on
                plot(obj.analpha_lag,obj.CN,'--','DisplayName','C_N(\alpha'')')
                warning('%s : CN curve was displayed for analytical alpha as the experimental alpha was not defined.',obj.name)
            end
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
            title(sprintf('%s ($\\dot{\\alpha} = %.2f ^{\\circ}$/s)',obj.name,obj.alphadot),'interpreter','latex')
        end
        function plotAlpha(obj)
            plotAlpha@AirfoilMotion(obj)
            if ~isempty(obj.alpha_CConset)
                plot(obj.t(obj.i_CConset),obj.alpha_CConset,'rx','DisplayName','\alpha_{ds,CC}')
            end  
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
            if ~isempty(obj.i_CConset)
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