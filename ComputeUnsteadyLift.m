function [CN_LB,CNk,CNprime,CNf,CNv,f,fp,fpp,alphaf,alphaE] = ComputeUnsteadyLift(pitching,airfoil,M,Tp,Tf,Tv)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Checks for parameters

%% Attached flow 

pitching.computeAttachedFlow(airfoil)

%% Leading edge separation

% Delayed normal coefficient
DeltaS = mean(diff(pitching.S));
Dp = zeros(size(pitching.CNp));
for n=2:length(pitching.CNp)
    Dp(n) =  Dp(n-1)*exp(-DeltaS/Tp) + (CNp(n)-CNp(n-1))*exp(-DeltaS/(2*Tp));
end 
% CNprime = CNslope*PitchingMotion.pitching.alpha_rad(1:length(Dp)) - Dp; % from CN steady modelled
CNprime = CNslope*analpha(1:length(Dp)) - Dp;

%% Trailing edge separation

f = seppoint(airfoil.steady,airfoil.steady.alpha); % alpha is PitchingMotion.alphaxp

% Kirchhoff law
CNk = Kirchhoff(airfoil.steady,airfoil.steady.alpha);


alphaf_rad = CNprime/CNslope;
alphaf = rad2deg(alphaf_rad);

fp = seppoint(airfoil.steady,alphaf); % effective separation point

Df=zeros(size(fp));
for n=2:length(fp)
    Df(n) = Df(n-1)*exp(-DeltaS/Tf) + (fp(n)-fp(n-1))*exp(-DeltaS/(2*Tf));
end

fpp = fp - Df;

%% Dynamic stall

% vortex normal coeff
KN = (1+sqrt(fpp)).^2/4;
Cv = CNC(1:length(KN)).*(1-KN);

CNv=zeros(size(Cv));
for n=2:length(Cv)
    CNv(n) = CNv(n-1)*exp(-DeltaS/Tv) + (Cv(n)-Cv(n-1))*exp(-DeltaS/(2*Tv));
end

if length(CNI)<length(CNC)
    CNf = ((1+sqrt(fpp))/2).^2.*CNC(1:length(CNI))+CNI;
else
    CNf = ((1+sqrt(fpp))/2).^2.*CNC+CNI(1:length(CNC));
end

% normal force coefficient
CN_LB = CNf + CNv;

end