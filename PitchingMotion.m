classdef PitchingMotion < handle
    properties
        name
        alpha
        alpha_rad
        analpha
        analpha_rad
        CN
        CC
        CL
        CD
        CNsteady
        mean_rad
        amp_rad
        freq
        omega
        phi % phase at t=0 in radians
        Ts
        t
        S % convective time
        k
        rt
        V
        M
        % Beddoes-Leishman
        % Beddoes constants
        A1 = 0.3;
        A2 = 0.7;
        b1 = 0.14;
        b2 = 0.53;
        % Attached flow behaviour
        CNI
        CNC
        CNp
        alphaE
        alphaE_rad
        
    end
    methods
        % convenient constructor with name/value pair of any attribute of
        % PitchingMotion
        function obj = PitchingMotion(varargin)
            p = inputParser;
            % Add name / default value pairs
            prop = properties('PitchingMotion'); % makes a cell array of all properties of the specified ClassName
            for k=1:length(prop)
                if strcmp(prop{k},'A1')
                    p.addParameter(prop{k},0.3);
                elseif strcmp(prop{k},'A2')
                    p.addParameter(prop{k},0.7);
                elseif strcmp(prop{k},'b1')
                    p.addParameter(prop{k},0.14);
                elseif strcmp(prop{k},'b2')
                    p.addParameter(prop{k},0.53);
                else
                    p.addParameter(prop{k},[]);
                end
            end
            p.parse(varargin{:}); % {:} is added to take the content of the cells
            for k=1:length(prop)
                eval(sprintf('obj.%s = p.Results.%s;',prop{k},prop{k}));
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
        function setSinus(obj,airfoil,mean_rad,amp_rad,freq)
            obj.freq = obj.k*obj.V/(pi*airfoil.c);
            obj.Ts = 1/(obj.freq*length(obj.alpha));
            obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
            if nargin<5
                freq = obj.freq;
            end
            if ~isempty(obj.amp_rad)
                warning('amp_rad value of PitchingMotion instance is not empty. Erasing with new value.')
            end
            if ~isempty(obj.mean_rad)
                warning('mean_rad value of PitchingMotion instance is not empty. Erasing with new value.')
            end
            obj.amp_rad = amp_rad;
            obj.mean_rad = mean_rad;
            obj.omega = 2*pi*freq;
            if diff(obj.alpha(1:2))>0 % if phi is in 1st or 2nd quadrant
                obj.phi = asin((obj.alpha(1)-obj.mean_rad)/obj.amp_rad);
            else % if phi is in 3rd or 4th quadrant
                obj.phi = pi - asin((obj.alpha(1)-obj.mean_rad)/obj.amp_rad);
            end
            obj.Ts = 1/(freq*length(obj.alpha));
            obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
            obj.analpha_rad = mean_rad + amp_rad*sin(obj.omega*obj.t-obj.phi);
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
            if isa(varargin{2},'SteadyCurve')
                steady = varargin{2};
                obj.CNsteady = interp1(steady.alpha,steady.CN,obj.alpha);
            else
                if length(varargin{2})==length(obj.alpha)
                    obj.CNsteady = varargin{2};
                else
                    error('CN and alpha must be of same length. Provide a SteadyCurve object for resampling.')
                end
            end
        end
        function plotAlpha(obj)
            figure
            plot(obj.t,obj.alpha,'DisplayName','\alpha')
            hold on
            plot(obj.t,obj.analpha,'DisplayName','ideal \alpha')
            grid on
            legend show
            xlabel('t (s)')
            ylabel('\alpha (deg)')
            title(obj.name)
        end
        function BeddoesLeishman(obj,airfoil)
            obj.computeAttachedFlow(airfoil);
            obj.computeTEseparation();
            obj.computeLEseparation();
            obj.computeDS();
        end
        function computeAttachedFlow(obj,airfoil)
            obj.S = 2*obj.V*obj.t/airfoil.c;
            obj.computeImpulsiveLift(airfoil,'analytical');
            obj.computeCirculatoryLift(airfoil);
            % effective angle of attack
            obj.alphaE_rad = obj.CNC/airfoil.steady.slope;
            obj.alphaE = rad2deg(obj.alphaE_rad);
            
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
        function computeCirculatoryLift(obj,airfoil)
            % Circulatory normal coeff
            DeltaS = mean(diff(obj.S));
            deltaalpha = diff(obj.alpha);
            X = zeros(size(deltaalpha));
            Y = zeros(size(deltaalpha));
            for n=2:deltaalpha
                X(n) = X(n-1)*exp(-obj.b1*betasquared*DeltaS) + obj.A1*deltaalpha(n)*exp(-obj.b1*betasquared*DeltaS/2);
                Y(n) = Y(n-1)*exp(-obj.b2*betasquared*DeltaS) + obj.A2*deltaalpha(n)*exp(-obj.b2*betasquared*DeltaS/2);
            end
            CNslope = Slope(airfoil.steady);
            obj.CNC = CNslope*(obj.alpha_rad(1:length(X))-X-Y);
            
        end
        function computeTEseparation(obj)
        end
        function computeLEseparation(obj)
        end
        function computeDS(obj)
        end
    end
end
