classdef AirfoilMotion < matlab.mixin.SetGet
    properties
        name
        model
        dataset
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
        CNqs % static + pi*c/2V*alphadot(t), according to Theodorsen
        CNth % CNqs + added mass 
        Ts
        t
        S % convective time
        rt % instantaneous red. pitch rate
        % experimental flow parameters
        V
        M % Mach number
        a % speed of sound
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
        % LE separation
        CNprime
        Dp
        % TE separation
        alphaf
        alphaf_rad
        CNk
        f
        fp
        fpp
        fppexp
        CNf
        CCf
        % Dynamic Stall
        tau_v
        Cv
        CNv
        St % Strouhal number
        CNv2
        CN_LB
        % Post-processing
        maxCN_LB
        maxCNk
        maxCNf
        maxCNv
        SmaxCN_LB
        SmaxCNk
        SmaxCNf
        SmaxCNv
        err % mean squared difference between the experimental data and the LB prediction       
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
            
            if isempty(obj.dataset)
                if isa(obj,'RampUpMotion')
                    obj.dataset = 'SH2019';
                    % speed of sound in water
                    obj.a = 1481; %m/s
                else
                    % speed of sound in air
                    obj.a = 343; %m/s
                    if obj.V==50
                        obj.dataset = 'simcos2008';
                    else
                        obj.dataset = 'LB1989';
                    end
                end
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
            % Sets the static CN value for each time step. Pass a SteadyCurve as an argument to automatically compute
            % the static CN for each angle of attack alpha(t) of the
            % current AirfoilMotion. Provide a vector describing CN(t) to
            % manually set the static CN values for each alpha(t).
            % Quasi-steady lift based on Theodorsen theory is also computed
            if isa(varargin{1},'SteadyCurve')
                steady = varargin{1};
                obj.CNsteady = interp1(steady.alpha,steady.CN,obj.alpha); % only the alphadot term
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
            n = length(obj.rt);
            CNsteady_pot = interp1(airfoil.steady.alpha(airfoil.steady.alpha<13),airfoil.steady.CN(airfoil.steady.alpha<13),obj.alpha,'linear','extrap');
            obj.CNqs = CNsteady_pot(1:n) + pi*obj.rt(1:n).*cosd(obj.alpha(1:n));        
            obj.computeAttachedFlow(airfoil,alphamode);
            obj.computeLEseparation(Tp,alphamode);
            obj.computeTEseparation(airfoil,Tf,'BL');
            obj.computeDS(airfoil,Tv,Tvl,'BL');
        end
        function BLBangga(obj,airfoil,Tp,Tf,Tv,Tvl,alphamode)
            obj.setCNsteady(airfoil.steady)
            obj.computeAttachedFlow(airfoil,alphamode);
            obj.computeLEseparation(Tp,alphamode);
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
            % Resample  to fit experimental unsteady alphas
            obj.tau1 = tau1;
            obj.tau2 = tau2;
            steady.setCN0();
            obj.set(steady)
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
            Tl = airfoil. c/obj.a;
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
            dalpha = diff(obj.alpha); % should be in degrees but only radians work!
            ddalpha = diff(dalpha); % same unit as above
            dalphadt = (dalpha(2:end) + dalpha(1:end-1))/(2*obj.Ts); % way smoother than Euler method (dalphadt = dalpha(1:end-1)/obj.Ts)!
            D = zeros(size(ddalpha));
            for n=2:length(ddalpha)
                D(n) = D(n-1)*exp(-obj.Ts/TlKalpha)+(dalphadt(n)-dalphadt(n-1))*exp(-obj.Ts/(2*TlKalpha));
            end
            obj.CNI = pi/180*4*TlKalpha/obj.M*(dalphadt-D); % pi/180 comes from a degree units remaining when checking all units. CNI should be dimensionless. Tl is in seconds and the parenthesis is in deg/s.
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
            obj.Dp = zeros(size(obj.CNp));
            for n=2:length(obj.CNp)
                obj.Dp(n) =  obj.Dp(n-1)*exp(-obj.DeltaS/obj.Tp) + (obj.CNp(n)-obj.CNp(n-1))*exp(-obj.DeltaS/(2*obj.Tp));
            end
            switch alphamode
                case 'analytical'
                    obj.CNprime = obj.CNp - obj.Dp;
                case 'experimental'
                    obj.CNprime = obj.CNp - obj.Dp;
            end
        end
        function computeTEseparation(obj,airfoil,Tf,model)
            obj.Tf = Tf;
            obj.f = seppoint(airfoil.steady,airfoil.steady.alpha);                       
            
            % Here a model for fp is selected
            obj.computeSepLag(airfoil,model)
            
            Df=zeros(size(obj.fp));
            for n=2:length(obj.fp)
                Df(n) = Df(n-1)*exp(-obj.DeltaS/Tf) + (obj.fp(n)-obj.fp(n-1))*exp(-obj.DeltaS/(2*Tf));
            end
            obj.fpp = obj.fp - Df;
            
            n = min([length(obj.CNI),length(obj.CNC)]);
            % slope_adapted = slope5.*(obj.CNprime<CNcrit)+slope10.*(obj.CNprime>CNcrit);
             % Kirchhoff law
            obj.CNk = airfoil.steady.slope.*(obj.alphaE(1:n).*((1+sqrt(obj.fpp(1:n)))/2).^2-airfoil.steady.alpha0);
            obj.CNf = obj.CNk + obj.CNI(1:n);
            % CNss goes to CNalpha*(1+sqrt(f_\infty)/2
            eta = 0.95;
            obj.CCf = eta*airfoil.steady.slope*obj.alphaE(1:n).^2.*sqrt(obj.fpp(1:n));
        end
        function computeSepLag(obj,airfoil,model)
            % change the static curve model here
            interpolate_fexp = false;
            switch model
                case {'BL','BLBangga'}
                    obj.model='BL';
                    obj.alphaf = obj.CNprime/airfoil.steady.slope; % effective separation point
                    obj.alphaf_rad = deg2rad(obj.alphaf);
                    if interpolate_fexp
                        obj.fp = interp1(airfoil.steady.alpha,airfoil.steady.fexp,obj.alphaf);
                    else % use the seppoint function fit evaluated in alphaf
                        obj.fp = seppoint(airfoil.steady,obj.alphaf);
                    end
                case 'BLSheng'
                    obj.model='BLSheng';
                    load(sprintf('/Users/lucas/Documents/EPFL/PDM/linfit_%s.mat',airfoil.name),'Talpha','alpha_ds0');
                    obj.computeAlphaLag(airfoil,Talpha)
                    % delta_alpha1 if r≥r0
                    da1 = alpha_ds0 - airfoil.steady.alpha_ss;
                    % delta_alpha1 if r<r0
                    da2 = (alpha_ds0 - airfoil.steady.alpha_ss)*obj.rt/airfoil.r0; % taking r(t) makes f' tend to finfty when r(t)=0
                    % concatenation of the two parts
                    use_rt = false;
                    if use_rt
                        delta_alpha1 = da1.*(obj.rt>=airfoil.r0)+da2.*(obj.rt<airfoil.r0);
                    else
                        delta_alpha1 = da1.*(obj.r>=airfoil.r0)+da2.*(obj.r<airfoil.r0);
                    end
                    n = min([length(obj.alpha_lag),length(delta_alpha1)]);
                    if interpolate_fexp
                        obj.fp = interp1(airfoil.steady.alpha,airfoil.steady.fexp,obj.alpha_lag(1:n)-delta_alpha1(1:n),'linear','extrap');
                    else
                        obj.fp = seppoint(airfoil.steady,obj.alpha_lag(1:n)-delta_alpha1(1:n));
                    end
                case 'BLSheng with expfit'
                    obj.model='expfit';
                    [~,Talpha] = Expfit(airfoil,obj);
                    obj.computeAlphaLag(airfoil,Talpha)
                    if interpolate_fexp
                        obj.fp = interp1(airfoil.steady.alpha,airfoil.fexp,obj.alpha_lag,'linear','extrap');
                    else
                        obj.fp = seppoint(airfoil.steady,obj.alpha_lag); % effective separation point
                    end
            end
            obj.computeExpSeparation(airfoil);
        end
        function computeExpSeparation(obj,airfoil)
            % computes the experimental separation point using inverted Kirchhof model,
            % ref: Leishman, Principles of Helicopter Aerodynamics 2nd Ed., eq. 7.106 page 405
            unbounded_fppexp = (2*sqrt((obj.CN(1:length(obj.CNv))-obj.CNv)./(airfoil.steady.slope*(obj.alpha(1:length(obj.CNv)) - airfoil.steady.alpha0)))-1).^2;
            obj.fppexp = max([zeros(size(unbounded_fppexp)),min([ones(size(unbounded_fppexp)), unbounded_fppexp],[],2)],[],2);
            obj.fppexp(obj.alpha<5) = 1;
        end
        function computeDS(obj,airfoil,Tv,Tvl,model)
            % Computes the final Beddoes-Leishman predicted CN for the instanciated pitching motion, after
            % having computed the attached-flow behavior and both TE and LE
            % separation with the homolog methods.
            secondary_vortex=false;
            obj.Tv = Tv;
            obj.Tvl = Tvl;
            % vortex normal coeff
            KN = 1/4*(1+sqrt(obj.fpp)).^2;
            n = min([length(obj.CNC),length(KN)]);
            % Cv defines how much lift will eventually go in the LEV. It is
            % proportional to the circulation and decreased by the
            % separation.
            obj.Cv = obj.CNC(1:n).*(1-KN(1:n));
            % Here a model for stall criterion is selected. tau_v is
            % incremented after stall has started
            obj.computeVortexTime(airfoil,model)
            
            obj.CNv=zeros(size(obj.Cv));
            for k=2:length(obj.Cv)
                if obj.tau_v(k) < Tvl % first one is 1!
                    % CNv augments proportionally to delta Cv and decays exponentially at the same time
                    obj.CNv(k) = obj.CNv(k-1)*exp(-obj.DeltaS/Tv) + (obj.Cv(k)-obj.Cv(k-1))*exp(-obj.DeltaS/(2*Tv));
                else
                    % CNv decays exponentially after tau_v has passed Tvl
                    obj.CNv(k) = obj.CNv(k-1)*exp(-obj.DeltaS/Tv);
                end
            end
            if secondary_vortex
                obj.CNv2=zeros(size(obj.Cv));
                obj.St =.19; % Strouhal number
                Tst = 2*(1-obj.fpp)/obj.St; % scalar?
                for k=2:length(obj.Cv)
                    if obj.tau_v(k) > Tst(k) & obj.tau_v(k) < Tst(k) + Tvl
                        % CNv augments proportionally to delta Cv and decays exponentially at the same time
                        obj.CNv2(k) = obj.CNv2(k-1)*exp(-obj.DeltaS/Tv) + (obj.Cv(k)-obj.Cv(k-1))*exp(-obj.DeltaS/(2*Tv));
                    else
                        % CNv decays exponentially after tau_v has passed Tvl
                        obj.CNv2(k) = obj.CNv2(k-1)*exp(-obj.DeltaS/Tv);
                    end
                end
            end
            n = min([n,length(obj.CNf)]);
            % normal force coefficient
            if secondary_vortex
                obj.CN_LB = obj.CNf(1:n) + obj.CNv(1:n) + obj.CNv2(1:n);
            else
                obj.CN_LB = obj.CNf(1:n) + obj.CNv(1:n);
            end
        end
        function estimateTp(obj,steady)
            % Estimation of the deficiency function Dp based on the equation CNprime(tds) = CNss = CN(tds) - Dp(tds)
            Dp_ds = obj.CNp(obj.i_CConset) - steady.CN(steady.alpha == steady.alpha_ss); % Dp(t_ds)
            figure
            plot(obj.S(1:length(obj.Dp)),obj.Dp,'LineWidth',2,'DisplayName','D_p')
            hold on
            plot(obj.S(obj.i_CConset),Dp_ds,'d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','D_p(t_{ds}) = C^p_{N,ds} - C_{N,ss}')
            xlabel('t_c')
            ylabel('D_p')
            grid on
            title(sprintf('r = %.3f, Tp = %.1f',obj.r,obj.Tp))
            legend('Location','NorthEast')
        end
        function removeAddedMass(obj)
            CNCppV = obj.CN(1:length(obj.CNI)) - obj.CNI;
            figure
            plot(obj.S,obj.CN,'LineWidth',2,'DisplayName','C_N^{exp}')
            hold on
            plot(obj.S(1:length(CNCppV)),CNCppV,'LineWidth',2,'DisplayName','C_N^{C+V}')
            plot(obj.S(1:length(obj.CNsteady)),obj.CNsteady,'--','LineWidth',2,'DisplayName','static')
            plot(obj.S(1:length(obj.CNI)),obj.CNf(1:length(obj.CNI))-obj.CNI,'LineWidth',2,'DisplayName','C_N^f')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','SouthEast')
            title(sprintf('r = %.3f',obj.r))
        end
        function Tvl = estimateTvl(obj)
            % gives an estimation for Tvl based on the time difference
            % between the max in CC and the max in CN (Boutet 2020)
            [~,imaxCN] = max(obj.CN);
            Tvl = obj.S(imaxCN) - obj.S(obj.i_CConset);
        end
        function computeVortexTime(obj,airfoil,model)
            % computes tau_v, the adimensional time variable that keeps track of the
            % vortex⁄⁄⁄ passing over the airfoil
            switch model
                case {'BL','BLBangga'}
                    % the airfoil stalls when a certain normal coeff is
                    % exceeded by CNprime
                    CNss = interp1(airfoil.steady.alpha,airfoil.steady.CN,airfoil.steady.alpha_ss);
                    % CNcrit = 1.147; % limit of CNds as r->0
                    CNcrit = CNss;
                    isstalled = obj.CNprime > CNcrit;
                case 'BLSheng'
                    % the airfoil stalls when a certain AoA is exceeded by
                    % alpha'
                    load(sprintf('/Users/lucas/Documents/EPFL/PDM/linfit_%s.mat',airfoil.name),'alpha_ds0','r0');
                    if exist('r0','var')
                        airfoil.r0 = r0;
                    end
                    alpha_crit = computeAlphaCrit(airfoil,alpha_ds0,obj.r);
                    isstalled = obj.alpha_lag > alpha_crit;
                case 'BLSheng with expfit'
                    % the airfoil stalls when the static stall angle is
                    % exceeded by alpha'
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
                Dalpha = zeros(size(obj.analpha));
                for k = 1:length(dalpha)
                    Dalpha(k+1) = Dalpha(k)*exp(-obj.DeltaS/Talpha) + dalpha(k)*exp(-obj.DeltaS/(2*Talpha));
                end
                obj.analpha_lag = obj.analpha - Dalpha;
            end
        end
        function findPeaks(obj)
            [maxCN,imaxCN] = max(obj.CN);
            [obj.maxCN_LB,imaxCN_LB] = max(obj.CN_LB);
            [obj.maxCNk,imaxCNk] = max(obj.CNk);
            [obj.maxCNf,imaxCNf] = max(obj.CNf);
            [obj.maxCNv,imaxCNv] = max(obj.CNv);
            SmaxCN = obj.S(imaxCN);
            obj.SmaxCN_LB = obj.S(imaxCN_LB);
            obj.SmaxCNk = obj.S(imaxCNk);
            obj.SmaxCNf = obj.S(imaxCNf);
            obj.SmaxCNv = obj.S(imaxCNv);
        end
        function save2Excel
            writematrix(obj.r,'../paramsweep.xlsx','Sheet',sprintf('adot = %1.1f',obj.alphadot),'Range','B4')
        end
        function save2struct
            save
        function plotAlpha(obj,mode)
            figure
            if ~isempty(obj.S)
                plot(obj.S,obj.alpha,'LineWidth',5,'DisplayName','exp')
                hold on
            else
                error('Convective time has not been computed.')
            end
            if ~isempty(obj.analpha)  && strcmp(mode,'full')
                plot(obj.S,obj.analpha,'--','LineWidth',5,'DisplayName','ideal')
                hold on
            end
            plot(obj.S(obj.i_CConset),obj.alpha_CConset,'d','MarkerFaceColor','k','DisplayName','\alpha_{ds,CC}')
            axis([0 Inf -5 35])
            if strcmp(mode,'full')
                legend('Location','SouthEast','FontSize',20)
            end
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
        function plotExp(obj)
            plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','dynamic')
            hold on
            plot(obj.S(1:length(obj.CNsteady)),obj.CNsteady,'--','LineWidth',2,'DisplayName','static')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            axis([0 obj.S(end) -0.1 2])
            ax = gca;
            ax.FontSize = 20;
        end
        function plotLB(obj,mode)
            figure
            switch mode
                case 'angle'
                    plot(obj.alpha(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
                    hold on
                    plot(obj.alpha(1:length(obj.CN_LB)),obj.CN_LB,'--','LineWidth',2,'DisplayName','LB')
                    plot(obj.alpha,obj.CNsteady)
                    xlabel('\alpha (°)')
                case 'convectime'
                    plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
                    hold on
                    plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
                    plot(obj.S(1:length(obj.CNsteady)),obj.CNsteady,'--','LineWidth',2,'DisplayName','static')
                    xlabel('t_c')
                otherwise
                    error('mode should be either ''angle'' or ''convectime''')
            end
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            axis([0 obj.S(end) -0.1 2])
            ax = gca;
            ax.FontSize = 20;
            title(sprintf('$r=%.3f,T_p=%.1f,T_f =%.1f,T_v=%.1f,T_{vl}=%.1f$',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex')
        end
        function plotVortex(obj)
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
            hold on
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
            plot(obj.S(1:length(obj.CNv)),obj.CNv,'LineWidth',2,'DisplayName','C_N^v')
            plot(obj.S(1:length(obj.CNC)),obj.CNC,'LineWidth',2,'DisplayName','C_N^C')
            plot(obj.S(1:length(obj.CNf)),obj.CNf,'LineWidth',2,'DisplayName','C_N^f')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            axis([0 obj.S(end) -0.1 2])
            ax = gca;
            ax.FontSize = 20;
            title(sprintf('$r=%.3f,T_p=%.1f,T_f =%.1f,T_v=%.1f,T_{vl}=%.1f$',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex')
            
            i_stop = find(obj.tau_v>obj.Tvl,1);            
            S_stop = obj.S(i_stop);
            
            figure
            subplot(211)
            plot(obj.S,obj.tau_v,'LineWidth',2,'DisplayName','\tau_v')
            hold on             
            xline(S_stop,'r--','DisplayName','\tau_v = T_{vl}','LineWidth',2)
            grid on
            legend('Location','SouthWest')
            
            subplot(212)
            plot(obj.S(1:length(obj.CNv)),obj.CNv,'LineWidth',2,'DisplayName','C_N^v')
            hold on 
            plot(obj.S(1:length(obj.Cv)),obj.Cv,'LineWidth',2,'DisplayName','C_v')
            xline(S_stop,'r--','DisplayName','\tau_v = T_{vl}','LineWidth',2)
            xlabel('t_c')
            grid on
            legend('Location','SouthWest')
        end
        function plotSteady(obj,airfoil)
            %             airfoil.steady.computeSlope(13)
            %             slope13 = airfoil.steady.slope;
            CN_steady = airfoil.steady.slope*((1+sqrt(airfoil.steady.f_inf))/2).^2*30;
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
            hold on
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
            plot(obj.S,CN_steady.*ones(size(obj.S)),'--','LineWidth',2,'DisplayName','steady-state')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            legend('Location','SouthEast','FontSize',20)
            axis([0 obj.S(end) -0.1 2])
            ax = gca;
            ax.FontSize = 20;
            title(sprintf('$r=%.3f,T_p=%.1f,T_f =%.1f,T_v=%.1f,T_{vl}=%.1f$',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex')
        end
        function plotShengLB(obj,airfoil)
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'DisplayName','exp','LineWidth',2)
            hold on
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'DisplayName','Sheng-LB','LineWidth',2)
            plot(obj.S(1:length(obj.CNf)),obj.CNf,'LineWidth',2,'DisplayName','C_N^f')
            plot(obj.S(1:length(obj.CNv)),obj.CNv,'LineWidth',2,'DisplayName','C_N^v')
            xlabel('t_c')
            ylabel('C_N')
            grid on
            ax = gca;
            ax.FontSize = 20;
            legend('Location','SouthEast','FontSize',20)
            load(sprintf('../linfit_%s.mat',airfoil.name),'Talpha');
            title(sprintf('$r=%.3f,T\\alpha=%.2f,Tf =%.2f,Tv=%.2f,Tvl=%.2f$',obj.r,Talpha,obj.Tf,obj.Tv,obj.Tvl),'FontSize',12,'interpreter','latex')
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
            plot(obj.t,obj.alpha,'LineWidth',2,'DisplayName','\alpha_{exp}','LineWidth',2)
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
                    plot(airfoil.steady.alpha,airfoil.steady.fexp,'DisplayName','f_{exp}','LineWidth',2)
                    hold on
                    plot(airfoil.steady.alpha,obj.f,'DisplayName','f','LineWidth',2)
                    plot(obj.alpha(1:length(obj.fp)),obj.fp,'DisplayName','f''','LineWidth',2)
                    plot(obj.alpha(1:length(obj.fpp)),obj.fpp,'DisplayName','f''''','LineWidth',2)
                    xlabel('\alpha (°)')
                    legend('Location','SouthWest','FontSize',20)
                case 'convectime'
                    plot(obj.S,interp1(airfoil.steady.alpha,airfoil.steady.fexp,obj.alpha),'DisplayName','f_{exp}','LineWidth',2)
                    hold on
                    if ~strcmp(obj.model,'BLBangga')
                        plot(obj.S,interp1(airfoil.steady.alpha,airfoil.steady.f,obj.alpha),'DisplayName','f','LineWidth',2)
                    end
                    plot(obj.S(1:length(obj.fp)),obj.fp,'DisplayName','f''','LineWidth',2)
                    plot(obj.S(1:length(obj.fpp)),obj.fpp,'DisplayName','f''''','LineWidth',2)
                    plot(obj.S(1:length(obj.fppexp)),obj.fppexp,'DisplayName','f''''_{exp}','LineWidth',2)
                    xlabel('t_c')
                    legend('Location','NorthEast','FontSize',20)
                otherwise
                    error('Either angle or convectime mode must be specified as an argument')
            end
            ax = gca;
            ax.FontSize = 20;
            axis([0 Inf 0 1])
            grid on
            ylabel('f')
            title(sprintf('r=%.3f,Tp=%.1f,Tf =%.1f,Tv=%.1f,Tvl=%.1f',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex')
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
            if ~any(strcmp(obj.model,{'BLSheng','Expfit'}))
                error('Neither Sheng or Expfit model has been run on that data.')
            end
            figure
            if ~isempty(obj.alpha)
                plot(obj.S,obj.alpha,'LineWidth',2,'DisplayName','\alpha')
                hold on
                plot(obj.S,obj.alpha_lag,'LineWidth',2,'DisplayName','\alpha''')
            end
            if ~isempty(obj.analpha)
                plot(obj.S,obj.analpha,'--','LineWidth',2,'DisplayName','ideal \alpha')
                hold on
                %plot(obj.S,obj.analpha_lag,'--','DisplayName','ideal \alpha''')
            end
            if ~isempty(obj.i_CConset)
                plot(obj.S(obj.i_CConset),obj.alpha_CConset,'rx','DisplayName','\alpha_{ds,CC}')
                hold on
                %plot(obj.S(obj.i_CConset),obj.alpha_lag(obj.i_CConset),'bx','DisplayName','\alpha''_{ds,CC}')
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
                plot(obj.S,obj.alphadot*(obj.t-t0-tau*(1-exp(-(obj.t-t0)/tau))),'--','LineWidth',2,'DisplayName','ideal ramp response')
                axis([0 Inf -1 40])
            end
            grid on
            xlabel('t_c')
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
        %         function plotRemainder(obj)
        %             remain = obj.CN(1:length(obj.CNI)) - obj.CNI - obj.CNsteady;
        %             figure
        %             plot(obj.S(1:length(remain)),remain)
        %             plot(obj.S,
        function plotFatma(obj)
            figure
            plot(obj.S(1:length(obj.CN)),obj.CN,'LineWidth',2,'DisplayName','exp')
            hold on 
            plot(obj.S(1:length(obj.CN_LB)),obj.CN_LB,'LineWidth',2,'DisplayName','LB')
            plot(obj.S(1:length(obj.CNqs)),obj.CNqs,'LineWidth',2,'DisplayName','C_N^{qs}')
            plot(obj.S(1:length(obj.CNv)),obj.CNv,'LineWidth',2,'DisplayName','C_N^v')
            plot(obj.S(1:length(obj.CNf)),obj.CNf,'LineWidth',2,'DisplayName','C_N^f')
            plot(obj.S(1:length(obj.CNsteady)),obj.CNsteady,'LineWidth',2,'DisplayName','static')
            grid on 
            xlabel('t_c')
            ylabel('C_N')
            legend('Location','SouthEast')
            ax = gca;
            ax.FontSize = 20;
            title(sprintf('$r=%.3f,T_p=%.1f,T_f =%.1f,T_v=%.1f,T_{vl}=%.1f$',obj.r,obj.Tp,obj.Tf,obj.Tv,obj.Tvl),'FontSize',14,'interpreter','latex')
        end
    end
end
