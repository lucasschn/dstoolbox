function CNk=Kirchhoff(steady,alpha,x)

if nargin > 2
    steady.S1 = x(1);
    steady.S2= x(2);
end

f = seppoint(steady,alpha);

% Kirchhoff law
CNk = steady.slope*((1+sqrt(f))/2).^2.*alpha + steady.CN0; % slope and alpha must have the same angle unit

end
