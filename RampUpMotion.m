classdef RampUpMotion < AirfoilMotion
    properties
        alpha_continuous_grow
        i_continuous_grow
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
            p = inputParser;
            % Add name / default value pairs
            mco = ?RampUpMotion;
            prop = mco.PropertyList; % makes a cell array of all properties of the specified ClassName
            for k=1:length(prop)
                if ~prop(k).Constant
                    p.addParameter(prop(k).Name,[]);
                end
            end
            p.parse(varargin{:}); % {:} is added to take the content of the cells
            for k=1:length(prop)
                if ~prop(k).Constant
                    obj.set(prop(k).Name,p.Results.(prop(k).Name))
                end
            end
            obj.fillProps() 
        end
        function isolateRamp(obj)
            dalpha = diff(obj.alpha);
            i_end = find((abs(dalpha/obj.Ts)<1e-2) .* (obj.alpha(1:end-1)>10),ceil(5/obj.Ts)); %during 1s
            i_end = i_end(end);
            t_end = obj.t(i_end);
            fprintf('Data will be cutoff at %.2fs \n',t_end)
            obj.alpha = obj.alpha(1:i_end);
            obj.alpha_rad = obj.alpha_rad(1:i_end);
            if ~isempty(obj.analpha)
                obj.analpha = obj.analpha(1:i_end);
            end
            if ~isempty(obj.analpha_rad)
                obj.analpha_rad = obj.analpha_rad(1:i_end);
            else 
                obj.analpha_rad = deg2rad(obj.analpha_rad);
            end
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
            % experiment. Is then used in Sheng algorithm to define the
            % vector alpha_ds.
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
        function findModelOnset(obj,airfoil)
            % finds the Sheng-predicted dynamic stall angle for a specific
            % time-evolution of alpha. Predicts alpha_onset from
            % alpha_lagonset from previously defined alpha_ds0.
            i_lagonset = find(obj.alpha_lag>airfoil.alpha_ds0,1);
            if isempty(i_lagonset)
                warning('The airfoil "%s" does not show stall in the experiment %s.',airfoil.name,obj.name)
            else
                if ~isempty(obj.alpha)
                obj.alpha_lagonset = obj.alpha(i_lagonset);
                obj.alpha_onset = interp1(obj.alpha_lag(obj.i_continuous_grow),obj.alpha_continuous_grow,obj.alpha_lagonset);
                elseif ~isempty(obj.analpha) % if alpha is empty
                obj.alpha_lagonset = obj.analpha(i_lagonset);
                [c,ia] = unique(obj.analpha_lag);
                obj.alpha_onset = interp1(c,obj.analpha(ia),obj.alpha_lagonset);
                else
                error('Impossible to define stall angle. %s has no angle of attack defined.',obj.name)
                end
            end
        end
        function plotCL(obj,xaxis)
            figure
            if (~exist('xaxis','var') || strcmp(xaxis,'alpha'))
                plot(obj.alpha,obj.CL)
            elseif strcmp(xaxis,'convectime')
                plot(obj.S,obj.CL)
            end
            grid on
            xlabel('\alpha (°)')
            ylabel('C_L (-)')
            title(sprintf('%s ($\\dot{\\alpha} = %.2f ^{\\circ}$/s)',obj.name,obj.alphadot),'interpreter','latex')
        end
        function plotCD(obj)
            figure
            plot(obj.alpha,obj.CD)
            grid on
            xlabel('\alpha (°)')
            ylabel('C_D (-)')
            title(sprintf('%s ($\\dot{\\alpha} = %.2f ^{\\circ}$/s)',obj.name,obj.alphadot),'interpreter','latex')
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
            title(sprintf('%s ($\\dot{\\alpha} = %.2f ^{\\circ}$/s)',obj.name,obj.alphadot),'interpreter','latex')
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