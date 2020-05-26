function CNk = kirchhoff(steady,alpha,x)

if nargin > 2
    steady.S1 = x(1);
    steady.S2= x(2);
end

f = seppoint(steady,alpha);

% Kirchhoff law, ref : Leishman, Principles of Helicopter Aerodynamics, 2nd
% Ed., eq. 7.105 page 405
CNk = steady.slope*((1+sqrt(f))/2).^2.*(alpha-steady.alpha0); % slope and alpha must have the same angle unit

end
