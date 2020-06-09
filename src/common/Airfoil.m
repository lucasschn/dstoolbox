classdef Airfoil < handle
    properties
        % This attributes must depend on the airfoil only. An airfoil has a
        % name and a chord. Its steady curve is not a property because it
        % depends on Re.
        steady % corresponding steady curve
        name
        c
        % Sheng 2008
        r0 = 0.01; % reduced pitch rate for bilinear transition  
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
        function b = b(obj)
            b = obj.c/2;
        end
    end
end
