function out = functionLB(param,static_alpha,static_CN,alpha_ss,raw,zero,LB,c,Tp,Tf,Tv,Tvl)

airfoil = Airfoil('airfoil',param.c);
airfoil.steady = SteadyCurve(static_alpha,static_CN,alpha_ss);
ramp = RampUpMotion('alpha',raw.alpha,'t',raw.t,'V',LB(c).U,'alphadot',LB(c).alphadot);

Cl = raw.Cl-zero.Cl;
Cd = raw.Cd-zero.Cd;
Cl = Cl - mean(Cl(1:50));
Cd = Cd - mean(Cd(1:50));
fs = 1/ramp.Ts;

% filter the load data
Cl_fff = myFilterTwice(Cl,fs);
Cd_fff = myFilterTwice(Cd,fs);
ramp.setCL(Cl_fff);
ramp.setCD(Cd_fff);

ramp.computeAirfoilFrame();
ramp.isolateRamp();

ramp.BeddoesLeishman(airfoil,Tp,Tf,Tv,Tvl,'experimental')
out = struct(ramp);
end