classdef AirfoilMotion < matlab.mixin.SetGet
    properties
        name
        model
        % experimental pitch angle evolutions
        alpha
        alpha_rad
        analpha
        analpha_rad
        alphadot_rad
        % exprimental load curves
        CL
        CD
        CN
        CC
        % experimental parts
        CNsteady % static value of CN(alpha(n)) according to the dynamic alpha(n)
        Ts
        t
        S % convective time
        rt % instantaneous red. pitch rate
        % experimental flow parameters
        V
        M % Mach number
        %% Beddoes-Leishman
        % Beddoes constants
        A1 = 0.3;
        A2 = 0.7;
        A3 = 0; % for third order according to Sheng 2008 and 2011
        b1 = 0.14;
        b2 = 0.53;
        b3 = 0;
        % Time constants
        DeltaS
        Tp
        Tf
        Tv
        Tvl
        % Attached flow behaviour
        CNI
        CNC
        CCC
        CNp
        alphaE
        alphaE_rad
        % TE separation
        alphaf
        alphaf_rad
        CNprime
        CNk
        f
        fp
        fpp
        CNf
        CCf
        % Dynamic Stall
        tau_v
        CNv        
        CN_LB
        %% Goman-Khrabrov
        tau1
        tau2
        x
        alpha_shift
        CN_GK
        %% Sheng
        alpha_lag
        analpha_lag
    end
    properties (Constant = true)
        % speed of sound
        a = 340.3;
    end
    methods
        function obj = AirfoilMotion(varargin)
            % convenient constructor with name/value pair of any attribute of
            % AirfoilMotion
            p = inputParser;
            % Add name / default value pairs
            mco = ?AirfoilMotion;
            prop = mco.PropertyList; % makes a cell array of all properties of the specified ClassName; % makes a cell array of all properties of the specified ClassName
            for k=1:length(prop)
                if ~prop(k).Constant && ~prop(k).HasDefault
                    p.addParameter(prop(k).Name,[]);
                end
            end
            p.parse(varargin{:}); % {:} is added to take the content of the cells
            % Add name / default value pairs
            for k=1:length(prop)
                if ~prop(k).Constant && ~prop(k).HasDefault
                    obj.set(prop(k).Name,p.Results.(prop(k).Name))
                end
            end
        end
        function fillProps(obj)
            % fills the properties that can be deduced from one another
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
                obj.M = obj.V/obj.a;
            end
            if isempty(obj.V) && ~isempty(obj.M)
                obj.V = obj.M*obj.a;
            end
        end
        function setName(obj,name)
            if nargin > 1
                obj.name = name;
            else
                obj.name = inputname(1);
            end
        end
        function setCN(obj,CN)
            if length(CN)==length(obj.alpha)
                obj.CN = CN;
            else
                error('CN and alpha must be of same length.')
            end
        end
        function setCNsteady(obj,varargin)
            % Sets the static CN value for each time step. Pass a StedayCurve as an argument to automatically compute
            % the static CN for each angle of attack alpha(t) of the
            % current AirfoilMotion. Provide a vector describing CN(t) to
            % manually set the static CN values for each alpha(t).
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
        function beta = beta(obj)
            if obj.M < 1
                beta = sqrt(1-obj.M^2);
            else
                beta = sqrt(obj.M^2-1);
            end
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
        function BeddoesLeishman(obj,airfoil,Tp,Tf,Tv,Tvl,alphamode)
            obj.model = 'LB';
            airfoil.steady.fitKirchhoff()
            obj.setCNsteady(airfoil.steady)
            obj.computeAttachedFlow(airfoil,alphamode);
            obj.computeLEseparation(Tp,alphamode);
            obj.computeTEseparation(airfoil,Tf,'BL');
            obj.computeDS(airfoil,Tv,Tvl,'BL');
        end
        function BLBangga(obj,airfoil,Tp,Tf,Tv,Tvl,alphamode)
            obj.setCNsteady(airfoil.steady)
            obj.computeAttachedFlow(airfoil,alphamode);
            obj.computeLEseparation(airfoil,Tp,alphamode);
            obj.computeTEseparation(airfoil,Tf,'BLBangga');
            obj.computeDS(airfoil,Tv,Tvl,'BLBangga');
        end
        function BLSheng(obj,airfoil,Tf,Tv,Tvl,alphamode)
            obj.model = 'Sheng-LB';
            airfoil.steady.fitKirchhoff()
            obj.computeAttachedFlow(airfoil,alphamode);
            obj.computeTEseparation(airfoil,Tf,'BLSheng');
            obj.computeDS(airfoil,Tv,Tvl,'BLSheng');
        end
        function BLexpfit(obj,airfoil,Tf,Tv,Tvl,alphamode)
            obj.model = 'Expfit';
            airfoil.steady.fitKirchhoff()
            obj.computeAttachedFlow(airfoil,alphamode);
            obj.computeTEseparation(airfoil,Tf,'BLSheng with expfit');
            obj.computeDS(airfoil,Tv,Tvl,'BLSheng with expfit');
        end
        function GomanKhrabrov(obj,steady,tau1,tau2)
            % Resample CNsteady to fit experimental unsteady alphas
            obj.tau1 = tau1;
            obj.tau2 = tau2;
            steady.setCN0();
            obj.setCNsteady(steady)
            steady.computeSlope();
            steady.computeSeparation();
            obj.alphadot_rad = diff(obj.alpha_rad)./diff(obj.t);
            obj.alpha_shift = obj.alpha_rad(1:length(obj.alphadot_rad))-tau2*obj.alphadot_rad;
            u = obj.alpha_shift;
            u(u<min(steady.alpha_rad)) = min(steady.alpha_rad);
            u(u>max(steady.alpha_rad)) = max(steady.alpha_rad);
            x0 = interp1(steady.alpha_rad,steady.fexp,u);
            sys = ss(-1/tau1,1/tau1,1,0);
            obj.x = lsim(sys,x0,obj.t(1:length(x0)));
            obj.CN_GK = steady.slope/4*(1+sqrt(obj.x)).^2+steady.CN0;
        end
        function computeAttachedFlow(obj,airfoil,alphamode)
            obj.S = 2*obj.V*obj.t/airfoil.c;
            obj.DeltaS = mean(diff(obj.S));
            obj.computeImpulsiveLift(airfoil,alphamode);
            obj.computeCirculatoryLift(airfoil,alphamode);
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
            Tl = airfoil.c/obj.a;
            TlKalpha = Kalpha*Tl;
            switch(alphamode)
                case 'experimental'
                    obj.computeExperimentalImpulsiveLift(TlKalpha)
                case 'analytical'
                    obj.computeAnalyticalImpulsiveLift(TlKalpha)
            end
        end
        function computeExperimentalImpulsiveLift(obj,TlKalpha)
            % experimental alphas, angles in rad
            dalpha = diff(obj.alpha_rad); % should be in degrees but only radians work!
            ddalpha = diff(dalpha); % same unit as above
            D = zeros(size(ddalpha));
            for n=2:length(ddalpha)
                D(n) = D(n-1)*exp(-obj.Ts/TlKalpha)+(ddalpha(n)/obj.Ts)*exp(-obj.Ts/(2*TlKalpha));
            end
            dalphadt = (dalpha(2:end) + dalpha(1:end-1))/(2*obj.Ts); % way smoother than Euler method (dalphadt = dalpha(1:end-1)/obj.Ts)!
            obj.CNI = 4*TlKalpha/obj.M*(dalphadt-D);
        end
        function computeCirculatoryLift(obj,airfoil,alphamode)
            % Circulatory normal coeff
            switch alphamode
                case 'experimental'
                    obj.computeExperimentalCirculatoryLift(airfoil);
                case 'analytical'
                    obj.computeAnalyticalCirculatoryLift(airfoil);
                otherwise
                    error('The type of alpha to be taken for computations has to be specified.')
            end
        end
        function computeExperimentalCirculatoryLift(obj,airfoil)
            deltaalpha = diff(obj.alpha); % deg
            rule = 'mid-point';            
            order = 2;
            switch order
                case 2
                    [X,Y] = obj.computeDuhamel(deltaalpha,rule);
                    % effective angle of attack
                    airfoil.steady.computeSlope(5)
                    obj.alphaE = (reshape(obj.alpha(1:length(X)),size(X))-X-Y);
                    obj.alphaE_rad = deg2rad(obj.alphaE);
                    obj.CNC = airfoil.steady.slope*obj.alphaE; % alpha is in degrees, slope is in 1/deg
                case 3
                    [X,Y,Z] = obj.computeDuhamelThirdOrder(deltaalpha);
                    % effective angle of attack
                    obj.alphaE = (reshape(obj.alpha(1:length(X)),size(X))-X-Y-Z);
                    obj.alphaE_rad = deg2rad(obj.alphaE);
                    obj.CNC = airfoil.steady.slope*obj.alphaE;
            end
        end
        function [X,Y] = computeDuhamel(obj,deltaalpha,rule)
            % Choose the best-suited algorithm for approximation of the
            % Duhamel integral. Rectangle is not recommended is the
            % step size in convective time is small, but can be 10
            % times faster to execute. Useful for numerous datapoints. (See section
            % 8.14.1 of Principles of Helicopter Aerodynamics by J.
            % Gordon Leishman)
            X = zeros(size(deltaalpha));
            Y = zeros(size(deltaalpha));
            switch rule
                case 'mid-point'  % approximation of Duhamel's integral using mid-point rule.
                    % Error of order DeltaS^3 (< 1% if both b1*DeltaS and b2*DeltaS are <0.25)
                    for n=2:length(deltaalpha)
                        X(n) = X(n-1)*exp(-obj.b1*obj.beta^2*obj.DeltaS) + obj.A1*deltaalpha(n)*exp(-obj.b1*obj.beta^2*obj.DeltaS/2);
                        Y(n) = Y(n-1)*exp(-obj.b2*obj.beta^2*obj.DeltaS) + obj.A2*deltaalpha(n)*exp(-obj.b2*obj.beta^2*obj.DeltaS/2);
                    end
                case 'rectangle' % approximation of the integral using rectangle rule.
                    % Faster, but error of order DeltaS (< 5% if both b1*DeltaS and b2*DeltaS are <0.05)
                    for n=2:length(deltaalpha)
                        X(n) = X(n-1)*exp(-obj.b1*obj.beta^2*obj.DeltaS) + obj.A1*deltaalpha(n);
                        Y(n) = Y(n-1)*exp(-obj.b2*obj.beta^2*obj.DeltaS) + obj.A2*deltaalpha(n);
                    end
            end
            
        end
        function [X,Y,Z] = computeDuhamelThirdOrder(obj,deltaalpha)
            % Implements a third order circulatory force model, following
            % Sheng 2008 & 2011. Changes Beddoes constants!
            obj.A1 = 0.165;
            obj.A2 = 0.335; 
            obj.A3 = 0.5;
            obj.b1 = 0.05; 
            obj.b2 = 0.222;         
            obj.b3 = 0.8/obj.M;
            X = zeros(size(deltaalpha));
            Y = zeros(size(deltaalpha));
            Z = zeros(size(deltaalpha));
            for n=2:length(deltaalpha)
                X(n) = X(n-1)*exp(-obj.b1*obj.beta^2*obj.DeltaS) + obj.A1*deltaalpha(n)*exp(-obj.b1*obj.beta^2*obj.DeltaS/2);
                Y(n) = Y(n-1)*exp(-obj.b2*obj.beta^2*obj.DeltaS) + obj.A2*deltaalpha(n)*exp(-obj.b2*obj.beta^2*obj.DeltaS/2);
                Z(n) = Z(n-1)*exp(-obj.b3*obj.beta^2*obj.DeltaS) + obj.A3*deltaalpha(n)*exp(-obj.b3*obj.beta^2*obj.DeltaS/2);
            end
            obj.A1 = 0.3;
            obj.A2 = 0.7;
            obj.b1 = 0.14;
            obj.b2 = 0.53;
        end
        function computeLEseparation(obj,Tp,alphamode)
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
                    obj.CNprime = obj.CNp - Dp;
                case 'experimental'
                    obj.CNprime = obj.CNp - Dp;
            end
        end
        function computeTEseparation(obj,airfoil,Tf,model)
            obj.Tf = Tf;
            obj.f = seppoint(airfoil.steady,airfoil.steady.alpha);
            
            % Kirchhoff law
            obj.CNk = kirchhoff(airfoil.steady,airfoil.steady.alpha);
            
            % Here a model for fp is selected
            obj.computeSepLag(airfoil,model)
            
            Df=zeros(size(obj.fp));
            for n=2:length(obj.fp)
                Df(n) = Df(n-1)*exp(-obj.DeltaS/Tf) + (obj.fp(n)-obj.fp(n-1))*exp(-obj.DeltaS/(2*Tf));
            end
            obj.fpp = obj.fp - Df;
            
            n = min([length(obj.CNI),length(obj.CNC)]);
            % TODO: Replace Kirchhoff model by interpolation of experimental
            % static curve
            airfoil.steady.computeSlope(5)
            slope5 = airfoil.steady.slope;
            airfoil.steady.computeSlope(10) % this has to be the same slope has used to computed fexp from CN in the steady curve. Otherwise, fp(end) will not match f(end)
%             slope10 = airfoil.steady.slope;
%             airfoil.steady.computeSlope(13)
%             slope13 = airfoil.steady.slope;
            % these two make discontinuities
            % slope_adapted = [slope5*ones(size(obj.alpha(obj.alpha<=5))); slope10*ones(size(obj.alpha(obj.alpha(1:n) > 5)))];
%             CNcrit = interp1(airfoil.steady.alpha,airfoil.steady.CN,airfoil.steady.alpha_ss);
%             slope_adapted = slope5.*(obj.CNprime<CNcrit)+slope10.*(obj.CNprime>CNcrit);
            slope_adapted = slope5;
            obj.CNf = slope_adapted.*(obj.alphaE(1:n).*((1+sqrt(obj.fpp(1:n)))/2).^2-airfoil.steady.alpha0)+obj.CNI(1:n);
            % CNss goes to CNalpha*(1+sqrt(f_\infty)/2
            eta = 0.95; 
            obj.CCf = eta*airfoil.steady.slope*obj.alphaE(1:n).^2.*sqrt(obj.fpp(1:n));
        end
        function computeSepLag(obj,airfoil,model)
            switch model
                case 'BL'
                    obj.alphaf = obj.CNprime/airfoil.steady.slope; % effective separation point
                    obj.alphaf_rad = deg2rad(obj.alphaf);
                    obj.fp = seppoint(airfoil.steady,obj.alphaf);
                case 'BLBangga'
                    obj.alphaf = obj.CNprime/airfoil.steady.slope; % effective separation point
                    obj.alphaf_rad = deg2rad(obj.alphaf);
                    n = min([length(obj.alphaf),length(obj.CNsteady)]);
                    % Bangga 2020, Eq. 33 (I think Bangga's equation is
                    % wrong, it should not be alphaf but alpha. Otherwise we would get f>1)
                    % take only alpha>0, we are not interested in knowing the values of the separation function on the negative side)                                    
                    visc_ratio = obj.CNsteady(obj.alpha(1:n)>0)./obj.CNprime(obj.alpha(1:n)>0); % here the zero-lift AoA for inviscid flow is needed.                    
                    if any(visc_ratio<0)
                        error('Viscous ratio cannot be negative.')
                    elseif any(visc_ratio>1)
                        error('Viscous ratio cannot be bigger than 1.')
                    else
                        obj.fp = (2*sqrt(visc_ratio)-1).^2;
                    end
                case 'BLSheng'
                    load(sprintf('/Users/lucas/Documents/EPFL/PDM/linfit_%s.mat',airfoil.name),'Talpha','alpha_ds0');
                    obj.computeAlphaLag(airfoil,Talpha)                                                         
                    da1 = alpha_ds0 - airfoil.steady.alpha_ss;                   
                    da2 = (alpha_ds0 - airfoil.steady.alpha_ss)*obj.rt/airfoil.r0;
                    delta_alpha1 = da1.*(obj.rt>=airfoil.r0)+da2.*(obj.rt<airfoil.r0); 
                    n = min([length(obj.alpha_lag),length(delta_alpha1)]);
                    obj.fp = seppoint(airfoil.steady,obj.alpha_lag(1:n)-delta_alpha1(1:n));
                case 'BLSheng with expfit'
                    [~,Talpha] = Expfit(airfoil,obj);
                    obj.computeAlphaLag(airfoil,Talpha)
                    obj.fp = seppoint(airfoil.steady,obj.alpha_lag); % effective separation point
            end
        end        
        function computeDS(obj,airfoil,Tv,Tvl,model)
            % Computes the final Beddoes-Leishman predicted CN for the instanciated pitching motion, after
            % having computed the attached-flow behavior and both TE and LE
            % separation with the homolog methods.
            obj.Tv = Tv; 
            obj.Tvl = Tvl;
            % vortex normal coeff
            KN = 1/4*(1+sqrt(obj.fpp)).^2;
            n = min([length(obj.CNC),length(KN)]);
            Cv = obj.CNC(1:n).*(1-KN(1:n));
            % Here a model for tau_v is selected
            obj.computeVortexTime(airfoil,model)
            
            obj.CNv=zeros(size(Cv));
            for k=2:length(Cv)
                if obj.tau_v(k) < Tvl
                    obj.CNv(k) = obj.CNv(k-1)*exp(-obj.DeltaS/Tv) + (Cv(k)-Cv(k-1))*exp(-obj.DeltaS/(2*Tv));
                else
                    obj.CNv(k) = obj.CNv(k-1)*exp(-obj.DeltaS/Tv);
                end
            end
            n = min([n,length(obj.CNf)]);
            % normal force coefficient
            obj.CN_LB = obj.CNf(1:n) + obj.CNv(1:n);
        end
        function computeVortexTime(obj,airfoil,model)
           switch model 
               case {'BL','BLBangga'}
                   CNcrit = interp1(airfoil.steady.alpha,airfoil.steady.CN,airfoil.steady.alpha_ss);
                   isstalled = obj.CNprime > CNcrit;
               case 'BLSheng'
                   load(sprintf('/Users/lucas/Documents/EPFL/PDM/linfit_%s.mat',airfoil.name),'alpha_ds0','r0');
                   if exist('r0','var')
                       airfoil.r0 = r0;
                   end
                   alpha_crit = computeAlphaCrit(airfoil,alpha_ds0,obj.r);
                   isstalled = obj.alpha_lag > alpha_crit;
               case 'BLSheng with expfit'
                   alpha_crit = airfoil.steady.alpha_ss;
                   isstalled = obj.alpha_lag > alpha_crit;
           end
           dt = mean(diff(obj.t));
           dalpha = diff(obj.alpha);
           obj.tau_v = zeros(size(obj.t));
           % Following Bangga 2020, Eq. 35 (first tau is always zero, last two points have no CNprime defined)
           for k = 2:length(obj.S)-2
               if isstalled(k)
                   obj.tau_v(k) = obj.tau_v(k-1)+0.45*dt/airfoil.c*obj.V;
               elseif dalpha(k) >= 0
                   obj.tau_v(k) = 0;
               else
                   obj.tau_v(k) = obj.tau_v(k-1);
               end
           end
        end
        function computeAlphaLag(obj,airfoil,Talpha)
            % computes the delayed angle of attack alpha_lag w.r.t. the
            % experimental angle of attack alpha. The time constant for the
            % delay is the corresponding Talpha(r), r being the reduced
            % pitch rate of the ramp-up motion.
            obj.S = obj.t*obj.V/(airfoil.c/2);
            if isempty(obj.DeltaS)
                obj.DeltaS = mean(diff(obj.S));
            end
            if ~isempty(obj.alpha) % compute alpha_lag from alpha using Eq.5
                dalpha = diff(obj.alpha);
                obj.alpha_lag = zeros(size(obj.alpha));
                Dalpha = zeros(size(obj.alpha));
                for k = 1:length(dalpha)
                    Dalpha(k+1) = Dalpha(k)*exp(-obj.DeltaS/Talpha) + dalpha(k)*exp(-obj.DeltaS/(2*Talpha));
                end
                obj.alpha_lag = obj.alpha - Dalpha;
                % compute analpha_lag from alpha_lag
                dalpha_lagdt = diff(obj.alpha_lag)./diff(obj.t);
                alphadot_lag = max(dalpha_lagdt);
                opts = optimset('Diagnostics','off','Display','off');
                analpha_lag0 = lsqcurvefit(@(x,xdata) alphadot_lag*xdata+x,0,obj.t(dalpha_lagdt>=5),obj.alpha_lag(dalpha_lagdt>=5),[],[],opts);
                obj.analpha_lag =  alphadot_lag*obj.t + analpha_lag0;
            elseif ~isempty(obj.analpha) % compute analpha_lag from analpha
                warning('analpha_lag was computed from analpha')
                dalpha = diff(obj.analpha);
                obj.analpha_lag = zeros(size(obj.analpha));
                Dalpha = zeros(size(obj.alpha));
                for k = 1:length(dalpha)
                    Dalpha(k+1) = Dalpha(k)*exp(-obj.DeltaS/Talpha) + dalpha(k)*exp(-obj.DeltaS/(2*Talpha));
                end
                obj.alpha_lag = obj.alpha - Dalpha;
            end
        end
        function plotAlpha(obj)
            figure
            if ~isempty(obj.alpha)
                plot(obj.S,obj.alpha,'LineWidth',5,'DisplayName','exp')
                hold on
            end
            if ~isempty(obj.analpha)
                plot(obj.S,obj.analpha,'--','LineWidth',5,'DisplayName','ideal')
                hold on
            end
            axis([0 Inf -5 35])
            legend('Location','SouthEast','FontSize',20)
            ax = gca;
            ax.FontSize = 20; 
            grid on
            xlabel('t_c')
            ylabel('\alpha (°)')
            title(obj.name)
        end
        function plotPolars(obj)
            figure
            suptitle(sprintf('Polar plots %s',obj.name))
            subplot(221)
            plot(obj.alpha,obj.CL)
            xlabel('\alpha (°)')
            ylabel('C_L')
            grid on
            subplot(222)
            plot(obj.alpha,obj.CD)
            xlabel('\alpha (°)')
            ylabel('C_D')
            grid on
            subplot(223)
            plot(obj.alpha,obj.CN)
            xlabel('\alpha (°)')
            ylabel('C_N')
            grid on
            subplot(224)
            if ~isempty(obj.CC)
                plot(obj.alpha,obj.CC)
                xlabel('\alpha (°)')
                ylabel('C_C')
                grid on
            end
        end
        function plotLB(obj,mode)
            figure
            switch mode
                case 'angle'
                    plot(obj.alpha(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
                    hold on
                    plot(obj.alpha(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
                    %plot(obj.alpha,obj.CNsteady)
                    xlabel('\alpha (°)')
                case 'convectime'
                    plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
                    hold on
                    plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
                    %plot(obj.S(1:length(obj.CNsteady)),obj.CNsteady,'LineWidth',2,'DisplayName','static')
                    xlabel('t_c')
            end
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            ax = gca;
            ax.FontSize = 20;
            title(sprintf('$r=%.3f,T_p=%.1f,T_f =%.1f,T_v=%.1f,T_{vl}=%.1f$',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex') 
        end
        function plotSteady(obj,airfoil)
            airfoil.steady.computeSlope(5)
            slope5 = airfoil.steady.slope;
            airfoil.steady.computeSlope(13)
            slope13 = airfoil.steady.slope;
            CN_steady = slope5*((1+sqrt(airfoil.steady.f_inf))/2).^2*30;
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
            hold on
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
            plot(obj.S,CN_steady.*ones(size(obj.S)),'--','LineWidth',2,'DisplayName','steady-state (Eq. 3.2)')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','NorthEast','FontSize',20)
            ax = gca;
            ax.FontSize = 20;
            title(sprintf('$r=%.3f,T_p=%.1f,T_f =%.1f,T_v=%.1f,T_{vl}=%.1f$',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex') 
        end
        function plotShengLB(obj,airfoil)
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'DisplayName','exp','LineWidth',2)
            hold on
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'DisplayName','Sheng-LB','LineWidth',2)
            plot(obj.S(1:length(obj.CNf)),obj.CNf,'DisplayName','C_N^f')
            plot(obj.S(1:length(obj.CNv)),obj.CNv,'DisplayName','C_N^v')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            load(sprintf('../linfit_%s.mat',airfoil.name),'Talpha');
            title(sprintf('$r=%.3f,T\\alpha=%.2f,Tf =%.1f,Tv=%.1f,Tvl=%.1f$',obj.r,Talpha,obj.Tf,obj.Tv,obj.Tvl),'interpreter','latex')
            ax = gca;
            ax.FontSize = 20;
            axis([0 obj.S(end) -0.03 2])
        end
        function plotLBExpfit(obj,airfoil,CN_Sheng)
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'DisplayName','exp','LineWidth',2)
            hold on            
            if nargin > 2
                plot(obj.S(1:length(CN_Sheng)),CN_Sheng,'DisplayName','LB-Sheng','LineWidth',2)
            end            
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'DisplayName','LB-Expfit','LineWidth',2)
            if obj.Tv ~=0 && isempty(CN_Sheng)
                plot(obj.S(1:length(obj.CNf)),obj.CNf,'DisplayName','C_N^f')
                plot(obj.S(1:length(obj.CNv)),obj.CNv,'DisplayName','C_N^v')
            end
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            [~,Talpha] = Expfit(airfoil,obj);
            title(sprintf('$r=%.3f,T\\alpha=%.1f,Tf =%.1f,Tv=%.1f,Tvl=%.1f$',obj.r,Talpha,obj.Tf,obj.Tv,obj.Tvl),'interpreter','latex')
            ax = gca;
            ax.FontSize = 20;
        end
        function plotGK(obj)
            figure
            plot(obj.alpha(1:length(obj.CN_GK)),obj.CN_GK,'LineWidth',2)
            xlabel('\alpha (°)')
            ylabel('C_N')
            grid on
            title(sprintf('\\tau_1 = %.3f \\tau_2 = %.2f',obj.tau1,obj.tau2))
        end
        function plotAlphas(obj)
            figure
            plot(obj.t,obj.alpha,'LineWidth',2,'DisplayName','\alpha_{xp}','LineWidth',2)
            hold on
            plot(obj.t(1:length(obj.alphaf)),obj.alphaf,'DisplayName','\alpha_f','LineWidth',2)
            plot(obj.t(1:length(obj.alphaE)),obj.alphaE,'DisplayName','\alpha_E','LineWidth',2)
            grid on
            ax = gca; 
            ax.FontSize = 20;
            legend('Location','Best','Fontsize',20)
            xlabel('time (s)')
            ylabel('\alpha (°)')
        end
        function plotSeparation(obj,airfoil,mode,save)
            figure
            switch mode
                case 'angle'
                    plot(airfoil.steady.alpha,airfoil.steady.fexp,'DisplayName','fexp','LineWidth',2)
                    hold on
                    plot(airfoil.steady.alpha,obj.f,'DisplayName','f','LineWidth',2)
                    plot(obj.alpha(1:length(obj.fp)),obj.fp,'DisplayName','f''','LineWidth',2)
                    plot(obj.alpha(1:length(obj.fpp)),obj.fpp,'DisplayName','f''''','LineWidth',2)
                    xlabel('\alpha (°)')
                case 'convectime'
                    plot(obj.S,interp1(airfoil.steady.alpha,airfoil.steady.fexp,obj.alpha),'DisplayName','fexp','LineWidth',2)
                    hold on
                    plot(obj.S,interp1(airfoil.steady.alpha,airfoil.steady.f,obj.alpha),'DisplayName','f (Eq. 2.1)','LineWidth',2)
                    if strcmp(obj.model,'LB')
                        disp_name = 'f'' (Eq. 2.10)';
                    else
                        disp_name = 'f'' (Eq. 2.24)';
                    end
                    plot(obj.S(1:length(obj.fp)),obj.fp,'DisplayName',disp_name,'LineWidth',2)
                    plot(obj.S(1:length(obj.fpp)),obj.fpp,'DisplayName','f'''' (Eq. 2.11)','LineWidth',2)
                    xlabel('t_c')
                otherwise 
                    error('Either angle or convectime mode must be specified as an argument')
            end
            legend('FontSize',20)
            title(sprintf('r=%.3f,Tp=%.1f,Tf =%.1f,Tv=%.1f,Tvl=%.1f',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'interpreter','latex')
            ax = gca;
            ax.FontSize = 20; 
            grid on
            ylabel('f')
            if  nargin>3 && save
                saveas(gcf,'../fig/f_curves.png')
            end
        end
        function plotStallOnset(obj,airfoil)
            % plots dynamic stall onset on CNprime curve based on
            % Leishman-Beddoes criterion CNprime > CNss
            CN_static_stall = interp1(airfoil.steady.alpha,airfoil.steady.CN,airfoil.steady.alpha_ss);
            figure
            plot(obj.alpha(1:length(obj.CNprime)),obj.CNprime,'DisplayName','CN''')
            hold on
            plot(obj.alpha,CN_static_stall*ones(size(obj.alpha)),'r--','DisplayName','C_N critical')
            grid on
            xlabel('\alpha (°)')
            ylabel('C_N''')
        end
        function plotSheng(obj,airfoil)
            figure
            plot(obj.rss,obj.alpha_ds,'.','DisplayName','\alpha_{ds} (exp)','MarkerSize',20)
            hold on
            plot(airfoil.r,airfoil.alpha_ds)
            grid on
            ylabel('\alpha_{ds} (°)','FontSize',20);
            ax = gca;
            ax.FontSize = 20;
        end
        function plotAlphaLag(obj,airfoil,model)
            figure
            if ~isempty(obj.alpha)
                plot(obj.t,obj.alpha,'DisplayName','\alpha')
                hold on
                plot(obj.t,obj.alpha_lag,'DisplayName','\alpha''')
            end
            if ~isempty(obj.analpha)
                plot(obj.t,obj.analpha,'--','DisplayName','ideal \alpha')
                hold on
                %plot(obj.t,obj.analpha_lag,'--','DisplayName','ideal \alpha''')
            end
            if ~isempty(obj.i_CConset)
                plot(obj.t(obj.i_CConset),obj.alpha_CConset,'rx','DisplayName','\alpha_{ds,CC}')
                hold on
                %plot(obj.t(obj.i_CConset),obj.alpha_lag(obj.i_CConset),'bx','DisplayName','\alpha''_{ds,CC}')
            end
            
            if nargin > 1
                switch model
                    case 'BLSheng'
                        load('/Users/lucas/Documents/EPFL/PDM/linfit_flatplate')
                        plot(obj.t(obj.i_CConset),computeAlphaCrit(airfoil,alpha_ds0,obj.r),'x','DisplayName','\alpha_{crit}')
                    case 'BLSheng with expfit'
                        plot(obj.t(obj.i_CConset),airfoil.steady.alpha_ss,'x','DisplayName','\alpha_{ss}')
                        [~,Talpha] = Expfit(airfoil,obj);                                                
                end
                t0 = interp1(obj.analpha,obj.t,0);
                tau = Talpha*airfoil.c/(2*obj.V);
                plot(obj.t,obj.alphadot*(obj.t-t0-tau*(1-exp(-(obj.t-t0)/tau))),'--','DisplayName','ideal ramp response')
                axis([0 Inf -1 40])
            end
            grid on
            xlabel('t (s)')
            ylabel('\alpha')
            title(sprintf('%s ($\\dot{\\alpha} = %.2f ^{\\circ}$/s)',obj.name,obj.alphadot),'interpreter','latex')
            legend('Location','SouthEast')
        end
        function plotCNprime(obj)
            n = length(obj.CNprime);
            figure
            plot(obj.S,obj.CNsteady,'DisplayName','static')
            hold on 
            plot(obj.S(1:n),obj.CNprime,'DisplayName','C_N''')
            grid on 
            xlabel('t_c')
            ylabel('C_N')
            legend('Location','SouthEast')
        end
    end
end
