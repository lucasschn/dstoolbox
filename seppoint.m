function f = seppoint(steady,alpha)
%SEPPOINT Returns the separation point of the boundary layer. AoAs in
%degrees.
% The function seppoint(S,dalpha) returns scalar or vector values of the boundary layer separation point f of
% in terms of position normalized by chord length x/c.
% The output is therefore between 0 and 1. The input arguments are an alpha
% scalar or vector giving the angle of attack in degrees at which f should be
% computed. The constants S1 and S2 together with the stall angle alpha1 depend on the airfoil and the flow
% and are determined using static data.

f1 = 1 - 0.3*exp((alpha-steady.alpha1)/steady.S1);
f2 = 0.04 +.66*exp((steady.alpha1-alpha)/steady.S2);
f = (alpha<=steady.alpha1).*f1 + (alpha>steady.alpha1).*f2;
end

