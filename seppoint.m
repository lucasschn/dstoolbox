function f = seppoint(steady,alpha)
%SEPPOINT Returns the separation point of the boundary layer. AoAs in
%degrees.
% The function seppoint(steady,dalpha) returns scalar or vector values of the boundary layer separation point f of
% in terms of position normalized by chord length x/c.
% The output is therefore between 0 and 1. The input arguments are an alpha
% scalar or vector giving the angle of attack in degrees at which f should be
% computed. The constants S1 and S2 together with the stall angle alpha_static_stall depend on the airfoil and the flow
% and are determined using static data. The coefficients are given by Sheng
% 2008.

f1 = 1 - .4*exp((alpha-steady.alpha_static_stall)/steady.S1);
f2 = .02 +.58*exp((steady.alpha_static_stall-alpha)/steady.S2);
f = (alpha<=steady.alpha_static_stall).*f1 + (alpha>steady.alpha_static_stall).*f2; % f1 is weighted by 1 if alpha<alpha_ss and f2 by 0 and vice-versa
end

