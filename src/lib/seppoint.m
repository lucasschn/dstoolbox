function f = seppoint(steady,alpha)
%SEPPOINT Returns the separation point of the boundary layer. AoAs in
%degrees.
% The function seppoint(steady,dalpha) returns scalar or vector values of the boundary layer separation point f of
% in terms of position normalized by chord length x/c.
% The output is therefore between 0 and 1. The input arguments are an alpha
% scalar or vector giving the angle of attack in degrees at which f should be
% computed. The constants S1 and S2 together with the stall angle alpha_ss depend on the airfoil and the flow
% and are determined using static data.

f1 = 1 - 0.3*exp((alpha-steady.alpha_ss)/steady.S1);
f2 = 0.04 +.66*exp((steady.alpha_ss-alpha)/steady.S2);
f = (alpha<=steady.alpha_ss).*f1 + (alpha>steady.alpha_ss).*f2; % f1 is weighted by 1 if alpha<alpha_ss and f2 by 0 and vice-versa
end
