classdef PitchingMotion < AirfoilMotion
    properties
        CNsteady
        mean_rad
        amp_rad
        freq
        omega
        phi % phase at t=0 in radians     
        k % reduced freq
        % Beddoes-Leishman
        DeltaS
        Tp
        Tf
        Tv
        % Attached flow behaviour
        CNI
        CNC
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
        % Dynamic Stall
        CNv
        CNf
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
        function setSinus(obj,airfoil,varargin)
            % setSinus uses mean_rad0, amp_rad0 and the optional argument freq0 as first guesses to
            % fin a sinusoid of the form alpha(t) = mean_rad +
            % amp_rad*sin(omega*t+phi) to the experimental angle of attack
            % belonging to the instanciated PitchingMotion object.
            mean_rad0 = varargin{1};
            amp_rad0 = varargin{2};
            if nargin < 5
                if isempty(obj.freq)
                    if isempty(obj.k) && isempty(obj.V)
                        error('k and V are unassigned, impossible to compute the sinus')
                    else
                        freq0 = obj.k*obj.V/(pi*airfoil.c);
                    end
                end
            else
                freq0 = varargin{3};
            end
            
            alpha_norm0 = (obj.alpha_rad(1)-mean_rad0)/amp_rad0;
            if alpha_norm0 > 1
                alpha_norm0 = 1;
            elseif alpha_norm0 < -1
                alpha_norm0 = -1;
            end
            if diff(obj.alpha_rad(1:2)) > 0 % if phi is in 1st or 4th quadrant
                phi0 = asin(alpha_norm0);
            else % if phi is in 2nd or 3rd quadrant
                phi0 = pi - asin(alpha_norm0);
            end
            omega0 = 2*pi*freq0;
            if nargin < 6
                pks_plus = findpeaks(obj.alpha);
                pks_minus = findpeaks(-obj.alpha);
                nbperiod = (length(pks_plus) + length(pks_minus))/2;
                T0 = 1/freq0; % duration of a period
                nb_samples_per_period = length(obj.alpha_rad)/nbperiod;
                Ts0 = T0/nb_samples_per_period; % sampling period
            else 
                fs0 = varargin{4};
                Ts0 = 1/fs0;
            end
            % we have to fit the time vector as a function of alpha
            x0 = [Ts0, omega0, mean_rad0, amp_rad0, phi0];
            alphaopt0 = x0(3) + x0(4)*sin(x0(2)*x0(1)*[0:(length(obj.alpha_rad)-1)]+x0(5));
            S = 0;
            for ks=1:length(alphaopt0)
                S = S + (alphaopt0(ks)-obj.alpha_rad(ks)).^2;
            end
            if S > 1 % This does not work properly, I don't know why
                opts = optimset('Display','iter');
                LB = 0.9*[0, 0, mean_rad0, amp_rad0, phi0];
                UB = [1e-9, 5e4, 1.2*mean_rad0, 1.2*amp_rad0, 5e4, 1.2*phi0];
                sinparams = lsqcurvefit(@(x,xdata) x(3)+x(4)*sin(x(2)*x(1)*xdata+x(5)),x0,reshape(0:(length(obj.alpha_rad)-1),size(obj.alpha_rad)),obj.alpha_rad,LB,UB,opts);
                if ~isempty(obj.amp_rad)
                    warning('amp_rad value of PitchingMotion instance is not empty. Erasing with new value.')
                end
                if ~isempty(obj.mean_rad)
                    warning('mean_rad value of PitchingMotion instance is not empty. Erasing with new value.')
                end
                obj.Ts = sinparams(1);
                obj.omega = sinparams(2);
                obj.mean_rad = sinparams(3);
                obj.amp_rad = sinparams(4);
                obj.phi = sinparams(5);
            else
                obj.Ts = x0(1);
                obj.omega = x0(2);
                obj.mean_rad = x0(3);
                obj.amp_rad = x0(4);
                obj.phi = x0(5);
            end
            obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha_rad)-1);
            obj.analpha_rad = obj.mean_rad + obj.amp_rad*sin(obj.omega*obj.t+obj.phi);
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
            obj.computeAttachedFlow(airfoil,'analytical');
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
            obj.alphaE_rad = deg2rad(obj.alphaE_rad);
            
            % Potential normal coefficient
            if length(obj.CNI)<length(obj.CNC)
                obj.CNp = obj.CNI + obj.CNC(1:length(obj.CNI));
            else
                obj.CNp = obj.CNI(1:length(obj.CNC))+obj.CNC;
            end
        end
        function computeImpulsiveLift(obj,airfoil,alphamode)
            % Impulsive (non-circulatory) normal coeff
            betasquared = 1-obj.M^2;
            Kalpha = 0.75/(1-obj.M+pi*betasquared*obj.M^2*(obj.A1*obj.b1+obj.A2*obj.b2));% time constant
            a = 340.3;
            Tl = airfoil.c/a;
            switch(alphamode)
                case 'experimental'
                    % angles in radian
                    deltaalpha = diff(obj.alpha);
                    ddalpha = diff(deltaalpha);
                    D = zeros(size(ddalpha));
                    for n=2:length(ddalpha)
                        D(n) = D(n-1)*exp(-obj.Ts/(Kalpha*Tl))+(ddalpha(n)/obj.Ts)*exp(-obj.Ts/(2*Kalpha*Tl));
                    end
                    % experimental alpha
                    obj.CNI = 4*Kalpha*Tl/obj.M*(deltaalpha(1:end-1)/obj.Ts-D);
                case 'analytical'
                    % analytical alphas
                    danalphadt = obj.amp_rad*cos(obj.omega*obj.t-obj.phi)*obj.omega;
                    ddanalphadt2 = obj.amp_rad*sin(obj.omega*obj.t-obj.phi)*obj.omega.^2;
                    D = zeros(size(ddanalphadt2));
                    for n=2:length(ddanalphadt2)
                        D(n) = D(n-1)*exp(-obj.Ts/(Kalpha*Tl))+(ddanalphadt2(n)*obj.Ts)*exp(-obj.Ts/(2*Kalpha*Tl));
                    end
                    % analytical alpha
                    obj.CNI = 4*Kalpha*Tl/obj.M*(danalphadt-D);
            end 
        end
        function computeCirculatoryLift(obj,airfoil,alphamode)
            % Circulatory normal coeff
            obj.DeltaS = mean(diff(obj.S));
            betasquared = 1-obj.M^2;
            switch alphamode
                case 'experimental'
                deltaalpha = diff(obj.alpha);
                case 'analytical'
                deltaalpha = diff(obj.analpha);  
                otherwise
                    error('The type of alpha to be taken for computations has to be specified.')
            end
            X = zeros(size(deltaalpha));
            Y = zeros(size(deltaalpha));
            for n=2:length(deltaalpha)
                X(n) = X(n-1)*exp(-obj.b1*betasquared*obj.DeltaS) + obj.A1*deltaalpha(n)*exp(-obj.b1*betasquared*obj.DeltaS/2);
                Y(n) = Y(n-1)*exp(-obj.b2*betasquared*obj.DeltaS) + obj.A2*deltaalpha(n)*exp(-obj.b2*betasquared*obj.DeltaS/2);
            end
            obj.CNC = airfoil.steady.slope*(reshape(obj.alpha(1:length(X)),size(X))-X-Y);
            
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
            
            obj.alphaf = obj.CNprime/airfoil.steady.slope; %twice the pulsation
            obj.alphaf_rad = deg2rad(obj.alphaf);
            
            obj.fp = seppoint(airfoil.steady,obj.alphaf); % effective separation point
            
            Df=zeros(size(obj.fp));
            for n=2:length(obj.fp)
                Df(n) = Df(n-1)*exp(-obj.DeltaS/Tf) + (obj.fp(n)-obj.fp(n-1))*exp(-obj.DeltaS/(2*Tf));
            end
            
            obj.fpp = obj.fp - Df;
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
            
            if length(obj.CNI)<length(obj.CNC)
                obj.CNf = ((1+sqrt(obj.fpp))/2).^2.*obj.CNC(1:length(obj.CNI))+obj.CNI;
            else
                obj.CNf = ((1+sqrt(obj.fpp))/2).^2.*obj.CNC+obj.CNI(1:length(obj.CNC));
            end
            
            % normal force coefficient
            obj.CN_LB = obj.CNf + obj.CNv;
        end
        function plotAlphas(obj)
            figure
            plot(obj.alpha,'DisplayName','\alpha_{xp}')
            hold on
            plot(obj.alphaf,'DisplayName','\alpha_f')
            plot(obj.alphaE,'DisplayName','\alpha_E')
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
            plot(airfoil.steady.alpha,steady.CN,'x')
            hold on
            plot(airfoil.steady.alpha,obj.CNk)
            plot(airfoil.steady.alpha1*ones(2),ylim,'r--')
            legend('static data','Kirchhoff')
            xlabel('\alpha')
            ylabel('C_N')
            grid on
            title('Lift attached flow')
        end
    end
end
