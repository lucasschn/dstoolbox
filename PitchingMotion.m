classdef PitchingMotion < AirfoilMotion
    properties
        CNsteady
        mean_rad
        amp_rad
        freq
        omega
        phi % phase at t=0 in radians
        k % reduced freq
        f_pts
        % LE separation
        CNprime
    end
    methods
        % convenient constructor with name/value pair of any attribute of
        % PitchingMotion
        function obj = PitchingMotion(varargin)
            p = inputParser;
            mco = ?PitchingMotion;
            prop = mco.PropertyList; % makes a cell array of all properties of the specified ClassName
            for k=1:length(prop)
                if ~prop(k).Constant
                    p.addParameter(prop(k).Name,[]);
                end
            end
            p.parse(varargin{:}); % {:} is added to take the content of the cells
            % Add name / default value pairs
            for k=1:length(prop)
                if ~prop(k).Constant
                    obj.set(prop(k).Name,p.Results.(prop(k).Name))
                end
            end
            obj.fillProps
        end

        function setSinus(obj,airfoil,varargin)
            % setSinus reconstructs the experimental angle of attack using a sinusoid of the form alpha(t) = mean_rad +
            % amp_rad*sin(omega*t+phi). If no optional argument is given,
            % the sinusoid is identified using a FFT. Otherwise, the
            % mean in radians (mean_rad), the amplitude of the main
            % harmonics in radians (amp_rad), the angular frequency (omega)
            % and the phase-lag (phi) must be given as argument.
            if isempty(varargin)
                [obj.mean_rad,obj.amp_rad,obj.f_pts,obj.phi]=findSinus(obj.alpha_rad);
                [obj.t,obj.Ts] = findTime(obj,airfoil.c);
                obj.omega = 2*pi*obj.f_pts/obj.Ts; % rev/si si/time   
            else
                obj.mean_rad = varargin{1};
                obj.amp_rad = varargin{2};
                obj.omega = varargin{3};
                if nargin==6 %(obj + airfoil + 4*varargin)
                    obj.phi = varargin{4};
                    [~,~,obj.f_pts] = findSinus(obj.alpha_rad);
                else
                    [~,~,obj.f_pts,obj.phi] = findSinus(obj.alpha_rad);
                end
                [obj.t,obj.Ts] = findTime(obj,airfoil.c);
            end
            obj.analpha_rad = reshape(obj.mean_rad + obj.amp_rad*sin(obj.omega*obj.t+obj.phi),size(obj.alpha_rad));
            obj.analpha = rad2deg(obj.analpha_rad);
        end
        function setPitchRate(obj,airfoil)
            % compute instantaneous reduced pitch rate
            if ~isempty(obj.omega)
                obj.omega = 2*pi*obj.freq;
            end
            obj.k = obj.omega*airfoil.c/(2*obj.V);
            obj.rt = obj.k*obj.amp_rad*cos(obj.omega*obj.t);
        end
        function setCN(obj,CN)
            if length(CN)==length(obj.alpha)
                obj.CN = CN;
            else
                error('CN and alpha must be of same length.')
            end
        end
        function setCNsteady(obj,varargin)
            if isa(varargin{1},'SteadyCurve')
                steady = varargin{1};
                obj.CNsteady = interp1(steady.alpha,steady.CN,obj.alpha);
            else
                if length(varargin{1})==length(obj.alpha)
                    obj.CNsteady = varargin{1};
                else
                    error('CN and alpha must be of same length. Provide a SteadyCurve object for resampling.')
                end
            end
        end
        function computeAnalyticalImpulsiveLift(obj,Talpha)
            % analytical alphas in degrees
            danalpha_rad = diff(obj.analpha_rad);
            ddanalpha_rad = diff(danalpha_rad); % not used, for debugging purposes
            danalphadt = obj.amp_rad*cos(obj.omega*obj.t+obj.phi)*obj.omega; % is equal to danalpha/Ts, but with length equal to that of alpha
            ddanalphadt2 = -(obj.amp_rad)*sin(obj.omega*obj.t+obj.phi)*obj.omega.^2; % is equal to diff(danalpha)*Ts^2, but with length equal to that of alpha
            D = zeros(size(ddanalphadt2));
            for n=2:length(ddanalphadt2)
                D(n) = D(n-1)*exp(-obj.Ts/(Talpha))+(ddanalphadt2(n)*obj.Ts)*exp(-obj.Ts/(2*Talpha));
            end
            obj.CNI = 4*Talpha/obj.M*(danalphadt-D);
        end
        function computeAnalyticalCirculatoryLift(obj,airfoil)
            % all angles in radians
            dalpha_rad = diff(obj.analpha_rad); % not used, for debugging purpose
            deltaalpha = obj.amp_rad*cos(obj.omega*obj.t+obj.phi)*obj.omega*obj.Ts; % deg, is equal to dalpha
            rule = 'mid-point';
            [X,Y] = obj.computeDuhamel(deltaalpha,rule);
            obj.CNC = airfoil.steady.slope*(reshape(obj.analpha(1:length(X)),size(X))-X-Y); % alpha is in degrees, slope is in 1/deg
        end
        function computeLEseparation(obj,airfoil,Tp,alphamode)
            % Computes the delayed normal coefficient, CNprime, depending
            % on the time constant Tp for a given airfoil undergoing the
            % instanciated pitching motion.
            obj.Tp = Tp;
            Dp = zeros(size(obj.CNp));
            for n=2:length(obj.CNp)
                Dp(n) =  Dp(n-1)*exp(-obj.DeltaS/obj.Tp) + (obj.CNp(n)-obj.CNp(n-1))*exp(-obj.DeltaS/(2*obj.Tp));
            end
            switch alphamode
                case 'analytical'
                    obj.CNprime = airfoil.steady.slope*obj.analpha(1:length(Dp)) - Dp; % we pretend the flow is attached over the whole alpha-range
                case 'experimental'
                    obj.CNprime = airfoil.steady.slope*obj.alpha(1:length(Dp)) - Dp; % we pretend the flow is attached over the whole alpha-range
            end
        end
        function computeSepLag(obj,airfoil)
            obj.alphaf = obj.CNprime/airfoil.steady.slope;
            obj.fp = seppoint(airfoil.steady,obj.alphaf);
        end
        function plotAttachedLift(obj,airfoil)
            figure;
            plot(airfoil.steady.alpha,airfoil.steady.CN,'x')
            hold on
            plot(airfoil.steady.alpha,obj.CNk)
            plot(airfoil.steady.alpha_static_stall*ones(2),ylim,'r--')
            legend('static data','Kirchhoff')
            xlabel('\alpha')
            ylabel('C_N')
            grid on
            title('Lift attached flow')
        end
        function plotStallOnset(obj,airfoil)
            CN_static_stall = interp1(airfoil.steady.alpha,airfoil.steady.CN,airfoil.steady.alpha_static_stall);
            figure
            plot(obj.alpha(1:length(obj.CNprime)),obj.CNprime,'DisplayName','CN''')
            hold on
            plot(obj.alpha,CN_static_stall*ones(size(obj.alpha)),'r--','DisplayName','C_N critical')
            grid on
            xlabel('\alpha (°)')
            ylabel('C_N''')
        end
        function plotLB(obj,xaxis)
            figure
            switch xaxis
                case 'alpha'
                    plot(obj.alpha(1:length(obj.CN)),obj.CN,'DisplayName','exp')
                    hold on 
                    plot(obj.alpha(1:length(obj.CN_LB)),obj.CN_LB,'DisplayName','LB')
                    xlabel('\alpha (°)')
                case 'convectime'
                    plot(obj.S(obj.t < 1/obj.freq),obj.CN(obj.t < 1/obj.freq),'DisplayName','exp')
                    hold on
                    plot(obj.S(obj.t < 1/obj.freq),obj.CN_LB(obj.t < 1/obj.freq),'DisplayName','LB')
                    xlabel('t_c (-)')
            end
            grid on 
            legend('Location','NorthEast','FontSize',20)
            ylabel('C_N')
            ax = gca;
            ax.FontSize = 20;
        end
    end
end