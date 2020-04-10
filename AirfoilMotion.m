classdef AirfoilMotion < handle
    properties
        name
        % experimental pitch angle evolutions
        alpha
        alpha_rad
        analpha
        analpha_rad
        % exprimental load curves
        CL
        CD
        CN
        CC
        % experimental parts
        Ts
        t
        S % convective time
        rt % instantaneous red. pitch rate
        % experimental flow parameters
        V
        M
    end
    methods
        function obj = AirfoilMotion(varargin)
            p = inputParser;
            % Add name / default value pairs
            prop = properties('AirfoilMotion'); % makes a cell array of all properties of the specified ClassName
            for k=1:length(prop)
                p.addParameter(prop{k},[]);
            end
        end
        function setName(obj)
            obj.name = inputname(1);
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
        function plotAlpha(obj)
            figure
            if ~isempty(obj.alpha)
                plot(obj.t,obj.alpha,'DisplayName','exp')
                hold on
            end
            if ~isempty(obj.analpha) 
                plot(obj.t,obj.analpha,'--','DisplayName','ideal')
                hold on
            end
            legend show
            grid on
            xlabel('t (s)')
            ylabel('\alpha (°)')
            title(obj.name)
        end
    end
end
