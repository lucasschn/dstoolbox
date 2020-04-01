classdef Airfoil < handle
    properties
        name
        c
        alpha_ds0
        D1
        Talpha
        r0 = 0.01 % reduced pitch rate at which linear fitting begins
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
        function Sheng(varargin)
            obj = varargin{1};
            % argument is a series of RampUpMotions
            r = [];
            alpha_ds = [];
            for k=2:nargin % first argument is self
                ramp = varargin{k};
                if isempty(ramp.r)
                    % with alphdot
                    ramp.setPitchRate(obj);
                end
                %r(k-1) = ramp.rt(ramp.i_CConset);
                if ramp.r>=obj.r0
                    r(end+1) = ramp.r;
                    alpha_ds(end+1) = ramp.alpha_CConset;
                end
            end
            p = polyfit(r,alpha_ds,1);
            obj.D1 = p(1);
            obj.alpha_ds0 = p(2);
            obj.Talpha = pi/180*obj.D1; % Talpha seems too low
%             obj.Talpha = 3.97;
            figure
            plot(r,alpha_ds,'o','DisplayName','\alpha_{ds} (exp)')
            hold on
            grid on
            xlabel('reduced pitch rate r (-)');
            ylabel('\alpha_{ds} (°)');
            title(obj.name)
            plot(r,polyval(p,r),'DisplayName','Linear fitting')
            alpha_lag_ds = zeros(size(alpha_ds));
            for k=1:length(alpha_ds)
                % retrieve which ramp corresponds to the current alpha_ds
                for kk=2:nargin
                    tmp = varargin{kk};
                    if tmp.r == r(k)
                        ramp = varargin{kk};
                    end
                end
                % compute alpha_lag and finds alpha_lagonset
                ramp.findModelOnset(obj);
                % looking for the value of alpha_lag(alpha) at the point alpha_ds
                if ramp.r >= obj.r0
                    if isempty(ramp.alpha)
                        alpha_lag_ds(k) = interp1(ramp.analpha,ramp.analpha_lag,alpha_ds(k));
                    elseif isempty(ramp.i_continous_grow)
                        alpha_lag_ds(k) = interp1(ramp.alpha,ramp.alpha_lag,alpha_ds(k));
                    else % if alpha_continuous_grow is defined
                        alpha_lag_ds(k) = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),alpha_ds(k));
                    end
                end
            end
            plot(r,alpha_lag_ds,'o','DisplayName','\alpha_{ds} (lagged)')
            plot(r,ones(size(r)).*obj.alpha_ds0,'--','DisplayName','\alpha_{ds,0}');
            legend('Location','SouthEast')
        end
    end
end
