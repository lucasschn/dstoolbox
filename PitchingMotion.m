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
        % Beddoes-Leishman
        DeltaS
        Tp
        Tf
        Tv
        % Attached flow behaviour
        CNI
        CNC
        CCC
        CNp
        alphaE
        alphaE_rad
        % LE separation
        CNprime
        % TE separation
        alphaf
        alphaf_rad
        CNk
        f
        fp
        fpp
        CNf
        CCf
        % Dynamic Stall
        CNv
        CN_LB
    end
    properties (Constant = true)
        % Beddoes constants
        A1 = 0.3;
        A2 = 0.7;
        b1 = 0.14;
        b2 = 0.53;
    end
    methods
        % convenient constructor with name/value pair of any attribute of
        % PitchingMotion
        function obj = PitchingMotion(varargin)
            obj@AirfoilMotion(varargin)
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
                    eval(sprintf('obj.%s = p.Results.%s;',prop(k).Name,prop(k).Name))
                end
            end
            if ~isempty(obj.alpha)
                obj.alpha_rad = deg2rad(obj.alpha);
            end
            if ~isempty(obj.Ts) && isempty(obj.t)
                obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
            end
            if isempty(obj.Ts) && ~isempty(obj.t)
                obj.Ts = mean(diff(obj.t));
            end
            if ~isempty(obj.V) && isempty(obj.M)
                a = 340.3;
                obj.M = obj.V/a;
            end
        end
        function b = b(obj)
            b = obj.c/2;
        end
        function beta = beta(obj)
            if obj.M < 1
                beta = sqrt(1-obj.M^2);
            else
                beta = sqrt(obj.M^2-1);
            end
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
                obj.phi = varargin{4};
                [~,~,obj.f_pts] = findSinus(obj.alpha_rad);
                
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
        function BeddoesLeishman(obj,airfoil,Tp,Tf,Tv)
            obj.computeAttachedFlow(airfoil,'experimental');
            obj.computeLEseparation(airfoil,Tp);
            obj.computeTEseparation(airfoil,Tf);
            obj.computeDS(Tv);
        end
        function computeAttachedFlow(obj,airfoil,alphamode)
            obj.S = 2*obj.V*obj.t/airfoil.c;
            obj.computeImpulsiveLift(airfoil,alphamode);
            obj.computeCirculatoryLift(airfoil,alphamode);
            % effective angle of attack
            obj.alphaE = obj.CNC/airfoil.steady.slope; % slope is in deg
            obj.alphaE_rad = deg2rad(obj.alphaE);
            % chord force
            obj.CCC = obj.CNC.*tan(obj.alphaE_rad);
            % Potential normal coefficient
            if length(obj.CNI)<length(obj.CNC)
                obj.CNp = obj.CNI + obj.CNC(1:length(obj.CNI));
            else
                obj.CNp = obj.CNI(1:length(obj.CNC))+obj.CNC;
            end
        end
        function computeImpulsiveLift(obj,airfoil,alphamode)
            % Impulsive (non-circulatory) normal coeff
            Kalpha = 0.75/(1-obj.M+pi*obj.beta*obj.M^2*(obj.A1*obj.b1+obj.A2*obj.b2)); % time constant
            a = 340.3;
            Tl = airfoil.c/a;
            switch(alphamode)
                case 'experimental'
                    % experimental alphas, angles in degrees
                    dalpha = diff(obj.alpha_rad); % should be degrees, but only radians work!
                    ddalpha = diff(dalpha); % same unit as above
                    D = zeros(size(ddalpha));
                    for n=2:length(ddalpha)
                        D(n) = D(n-1)*exp(-obj.Ts/(Kalpha*Tl))+(ddalpha(n)/obj.Ts)*exp(-obj.Ts/(2*Kalpha*Tl));
                    end          
                    obj.CNI = 4*Kalpha*Tl/obj.M*(dalpha(1:end-1)/obj.Ts-D);
                case 'analytical'
                    % analytical alphas in degrees
                    danalpha = diff(obj.analpha);
                    ddanalpha = diff(danalpha); % not used, for debugging purposes
                    danalphadt = rad2deg(obj.amp_rad)*cos(obj.omega*obj.t+obj.phi)*obj.omega; % is equal to danalpha/Ts, but with length equal to that of alpha
                    ddanalphadt2 = -rad2deg(obj.amp_rad)*sin(obj.omega*obj.t+obj.phi)*obj.omega.^2; % is equal to diff(danalpha)*Ts^2, but with length equal to that of alpha
                    D = zeros(size(ddanalphadt2));
                    for n=2:length(ddanalphadt2)
                        D(n) = D(n-1)*exp(-obj.Ts/(Kalpha*Tl))+(ddanalphadt2(n)*obj.Ts)*exp(-obj.Ts/(2*Kalpha*Tl));
                    end
                    obj.CNI = 4*Kalpha*Tl/obj.M*(danalphadt-D);
%                 case 'statespace'
%                     Kalpha
%                     A = -1/(Kalpha*Talpha);
%                     B = [1 0];
%                     C = [c13 0];
%                     D = []
%                     obj.CNI = 4/M*
            end
        end
        function computeCirculatoryLift(obj,airfoil,alphamode)
            % Circulatory normal coeff
            obj.DeltaS = mean(diff(obj.S));
            switch alphamode
                case 'experimental'
                    deltaalpha = diff(obj.alpha_rad); % should be in degrees, but only radians work
                case 'analytical'
                    % all angles in degrees
                    dalpha = diff(obj.analpha); % not used, for debugging purpose
                    deltaalpha = rad2deg(obj.amp_rad)*cos(obj.omega*obj.t+obj.phi)*obj.omega*obj.Ts; % deg, is equal to dalpha
                otherwise
                    error('The type of alpha to be taken for computations has to be specified.')
            end
            X = zeros(size(deltaalpha));
            Y = zeros(size(deltaalpha));
            for n=2:length(deltaalpha)
                % Choose the best-suited algorithm for approximation of the
                % Duhamel integral. Rectangle is not recommended is the
                % step size in convective time is small, but can be 10
                % times faster to execute. Useful for numerous datapoints. (See section
                % 8.14.1 of Principles of Helicopter Aerodynamics by J.
                % Gordon Leishman)
                rule = 'mid-point'; 
                switch rule
                    case 'mid-point'  % approximation of Duhamel's integral using mid-point rule.
                        % Error of order DeltaS^3 (< 1% if both b1*DeltaS and b2*DeltaS are <0.25)
                        X(n) = X(n-1)*exp(-obj.b1*obj.beta^2*obj.DeltaS) + obj.A1*deltaalpha(n)*exp(-obj.b1*obj.beta^2*obj.DeltaS/2);
                        Y(n) = Y(n-1)*exp(-obj.b2*obj.beta^2*obj.DeltaS) + obj.A2*deltaalpha(n)*exp(-obj.b2*obj.beta^2*obj.DeltaS/2);
                    case 'rectangle' % approximation of the integral using rectangle rule.
                        % Faster, but error of order DeltaS (< 5% if both b1*DeltaS and b2*DeltaS are <0.05)
                        X(n) = X(n-1)*exp(-obj.b1*obj.beta^2*obj.DeltaS) + obj.A1*deltaalpha(n);
                        Y(n) = Y(n-1)*exp(-obj.b2*obj.beta^2*obj.DeltaS) + obj.A2*deltaalpha(n);
                end
            end
            obj.CNC = airfoil.steady.slope*(reshape(obj.alpha(1:length(X)),size(X))-X-Y); % alpha is in degrees, slope is in 1/deg  
        end
        function computeLEseparation(obj,airfoil,Tp)
            % Computes the delayed normal coefficient, CNprime, depending
            % on the time constant Tp for a given airfoil undergoing the
            % instanciated pitching motion.
            obj.Tp = Tp;
            Dp = zeros(size(obj.CNp));
            for n=2:length(obj.CNp)
                Dp(n) =  Dp(n-1)*exp(-obj.DeltaS/obj.Tp) + (obj.CNp(n)-obj.CNp(n-1))*exp(-obj.DeltaS/(2*obj.Tp));
            end
            obj.CNprime = airfoil.steady.slope*obj.analpha(1:length(Dp)) - Dp; % we pretend the flow is attached over the whole alpha-range
        end
        function computeTEseparation(obj,airfoil,Tf)
            obj.Tf = Tf;
            obj.DeltaS = mean(diff(obj.S));
            obj.f = seppoint(airfoil.steady,airfoil.steady.alpha);
            
            % Kirchhoff law
            obj.CNk = Kirchhoff(airfoil.steady,airfoil.steady.alpha);
            
            obj.alphaf = obj.CNprime/airfoil.steady.slope;
            obj.alphaf_rad = deg2rad(obj.alphaf);
            
            obj.fp = seppoint(airfoil.steady,obj.alphaf); % effective separation point
            
            Df=zeros(size(obj.fp));
            for n=2:length(obj.fp)
                Df(n) = Df(n-1)*exp(-obj.DeltaS/Tf) + (obj.fp(n)-obj.fp(n-1))*exp(-obj.DeltaS/(2*Tf));
            end
            
            obj.fpp = obj.fp - Df;
            
            if length(obj.CNI)<length(obj.CNC)
                obj.CNf = ((1+sqrt(obj.fpp))/2).^2.*obj.CNC(1:length(obj.CNI))+obj.CNI;
            else
                obj.CNf = ((1+sqrt(obj.fpp))/2).^2.*obj.CNC+obj.CNI(1:length(obj.CNC));
            end
            
            eta = 0.95;
            obj.CCf = eta*airfoil.steady.slope*obj.alphaE(1:length(obj.fpp)).^2.*sqrt(obj.fpp);
        end
        function computeDS(obj,Tv)
            % Computes the final Beddoes-Leishman predicted CN for the instanciated pitching motion, after
            % having computed the attached-flow behavior and both TE and LE
            % separation with the homolog methods.
            obj.Tv = Tv;
            
            % vortex normal coeff
            KN = (1+sqrt(obj.fpp)).^2/4;
            Cv = obj.CNC(1:length(KN)).*(1-KN);
            
            obj.CNv=zeros(size(Cv));
            for n=2:length(Cv)
                obj.CNv(n) = obj.CNv(n-1)*exp(-obj.DeltaS/Tv) + (Cv(n)-Cv(n-1))*exp(-obj.DeltaS/(2*Tv));
            end
            
            % normal force coefficient
            obj.CN_LB = obj.CNf + obj.CNv;
        end
        function plotAlphas(obj)
            figure
            plot(obj.t,obj.alpha,'LineWidth',2,'DisplayName','\alpha_{xp}')
            hold on
            plot(obj.t(1:length(obj.alphaf)),obj.alphaf,'DisplayName','\alpha_f')
            plot(obj.t(1:length(obj.alphaE)),obj.alphaE,'DisplayName','\alpha_E')
            grid on
            legend('Location','Best')
            xlabel('time (s)')
            ylabel('\alpha (°)')
        end
        function plotSeparation(obj,airfoil,mode,save)
            figure
            switch mode
                case 'normal'
                    plot(airfoil.steady.alpha,obj.f,'DisplayName','f')
                    hold on
                    plot(obj.alpha(1:length(obj.fp)),obj.fp,'DisplayName','f''')
                    plot(obj.alpha(1:length(obj.fpp)),obj.fpp,'DisplayName','f''''')
                case 'log'
                    semilogy(airfoil.steady.alpha,airfoil.steady.fexp,'DisplayName','experimental')
                    hold on
                    semilogy(airfoil.steady.alpha,obj.f,'DisplayName','model')
            end
            legend show
            title(sprintf('%s (Tp = %.2f, Tf = %.2f, Tv = %.2f)',obj.name,obj.Tp,obj.Tf,obj.Tv))
            grid on
            xlabel('\alpha (°)')
            ylabel('separation point (x/c)')
            if save
                saveas(gcf,'fig/f_curves.png')
            end
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
        
    end
end