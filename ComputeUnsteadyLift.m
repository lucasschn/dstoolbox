function computeUnsteadyLift(pitching,airfoil,Tp,Tf,Tv)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Checks for parameters

%% Attached flow 
pitching.computeAttachedFlow(airfoil,'analytical')

%% Leading edge separation
pitching.computeLEseparation(airfoil,Tp)

%% Trailing edge separation
pitching.computeTEseparation(airfoil,Tf)

%% Dynamic stall
pitching.computeDS(Tv)

end