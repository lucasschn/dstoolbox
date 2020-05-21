classdef Joukowsky < handle
    properties
        mux = -0.2;
        muy = 0.2;
        z
        x
        y
        zeta
        khi
        eta
        r
        theta
    end
    methods
        function Joukowksky(mux,muy)
            obj.mux = mux;
            obj.muy = muy;
            obj.theta = linspace(0,2*pi,100);
            obj.r = sqrt((1-mux).^2+muy.^2);
            
            obj.khi = obj.r*cos(obj.theta) + mux;
            obj.eta = obj.r*sin(obj.theta) + muy;
            obj.zeta = obj.khi + 1i*obj.eta;
            obj.z = obj.zeta + obj.zeta.^(-1);
            obj.x = real(obj.z);
            obj.y = imag(obj.z);
            % alternative method
            % x = khi.*(khi.^2+eta.^2+1)./(khi.^2+eta.^2);
            % y = eta.*(khi.^2+eta.^2-1)./(khi.^2+eta.^2);
        end
        function plotCircle(obj)
           figure
           plot(obj.khi,obj.eta)
           hold on
           plot(obj.mux,obj.muy,'rx')
           axis([-1.5,1.5,-1.5,1.5])
           xlabel('\chi')
           ylabel('\eta')
           grid on
            
        end
        function plotAirfoil(obj)
            figure
            plot(obj.x,obj.y)
            axis([-3,3,-3,3])
            xlabel('x')
            ylabel('y')
            grid on
        end
    end
end