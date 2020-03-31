classdef PitchingMotion < handle
    properties
        alpha
        alpha_rad
        analpha
        analpha_rad
        CN
        CNsteady
        mean_rad
        amp_rad
        Ts
        t
        k
        rt
        V
    end
    methods
        % constructor (not possible to overload in Matlab)
        function obj = PitchingMotion(varargin)
            if length(varargin{1}) > 1
                obj.alpha = varargin{1};
                obj.alpha_rad = deg2rad(obj.alpha);
                
                if nargin == 2 %(alpha,Ts)
                    obj.Ts = varargin{2};
                    
                elseif nargin == 3 %(alpha,CN,Ts)
                    obj.CN = varargin{2};
                    obj.Ts = varargin{3};
                end
                obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
            else
                error('The first argument must be the measured incidence angles vector.')
            end
        end
        function setSinus(obj,mean_rad,amp_rad,freq)
            obj.amp_rad = amp_rad;
            obj.mean_rad = mean_rad;
            obj.Ts = 1/(freq*length(obj.alpha));
            obj.t = 0:obj.Ts:obj.Ts*(length(obj.alpha)-1);
            obj.analpha_rad = mean_rad + amp_rad*sin(2*pi*freq*obj.t);
        end
        function setPitchRate(obj,freq,airfoil)
            % compute instantaneous reduced pitch rate
            omega = 2*pi*freq;
            obj.k = omega*airfoil.c/(2*obj.V);
            obj.rt = obj.k*obj.amp_rad*cos(omega*obj.t);
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
    end
end
