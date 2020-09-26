function out = functionLB(input,LBcoeffs)

airfoil = Airfoil('airfoil',input.c);
airfoil.steady = SteadyCurve(input.static.alpha,input.static.Cn,input.alpha_ss);
ramp = RampUpMotion('alpha',input.dyn.alpha,'t',input.dyn.t,'V',input.U,'alphadot',input.alphadot);
ramp.setCL(input.dyn.Cl)
ramp.setCD(input.dyn.Cd)
ramp.computeAirfoilFrame();
% uncomment this if you don't want the whole signal but just the ramp part
%  ramp.isolateRamp();

ramp.BeddoesLeishman(airfoil,LBcoeffs.Tp,LBcoeffs.Tf,LBcoeffs.Tv,LBcoeffs.Tvl,'experimental')
out = struct(ramp);
end