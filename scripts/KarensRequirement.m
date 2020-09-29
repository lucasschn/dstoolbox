% This script is made to test that functionLB() in the three use cases,
% rampup motion, pitching motion and general motion.

% Author : Lucas Schneeberger
% Date : 29.09.2020

clear all
close all
clc

% now testing for simcos data
test = 'simcos';

%% Static properties
switch test
    case {'rampup','general'}
        static = load(fullfile('..','data','static_flatplate.mat'));
        input.static.alpha = static.alpha;
        input.static.Cn = static.CN;
        input.alpha_ss = 13;
    case 'simcos'
        load(fullfile('..','static_corr'))
        input.static.alpha = mA;
        input.static.Cn = mCl_corr;
    otherwise
        error('Your test case have not been recognized. The three options are rampup, simcos, and general.')
end

%% Defining model parameters

% Leishman-Beddoes constants
LBcoeffs.Tp = 3;
LBcoeffs.Tf = 3;
LBcoeffs.Tv = 1;
LBcoeffs.Tvl = 1;

%% Loading dynamic data

switch test
    case 'simcos' % for pitching motion
        run labbook_simcos.m
        c=13;
        data = load(pressuredata(c));
        % Drag missing, compute the pressure drag
        CD = zeros(size(data.Cl));
        for k=1:length(data.Cl)
            CD(k) = sum(data.Cp(k,:)*data.xk);
        end
        load(fullfile('..','dynamic_corr'))
        input.U = param.U0;
        input.c = param.c;
        n = length(data.Cl);
        TS = 1/LB(c).FS;
        input.dyn.t = 0:TS:(n-1)*TS;
        input.dyn.alpha = Alpha;
        input.dyn.Cl = Cl_corr;
        input.dyn.Cd = zeros(size(Cl_corr));

        % sinus parameters
        input.freq = LB(c).fosc;
        input.mean_rad = deg2rad(LB(c).alpha_0);
        input.amp_rad = deg2rad(LB(c).alpha_1);
        out = functionLB(input,LBcoeffs,'pitching');
    case 'rampup' % for rampup motion
        run labbook.m
        c = 71;
        load(loadmat(LB(c).ms,LB(c).mpt),'raw','zero')
        input.U=LB(c).U;
        input.c = param.c;
        input.dyn.t = raw.t;
        input.dyn.alpha = raw.alpha;
        input.dyn.Cl = raw.Cl-mean(raw.Cl(1:50));
        input.dyn.Cd = raw.Cd-mean(raw.Cd(1:50));

        % enter here the pitch rate
        input.alphadot=0.2; % deg/s

        out = functionLB(input,LBcoeffs,'ramp');
    case 'general'
        input.U = 50;
        input.c =0.15;
        n = 100;
        Ts = 0.1;
        t = 0:Ts:(n-1)*Ts;
        input.dyn.t = t ;
        input.dyn.alpha = 10+20*sind(4*t)+10*sind(10*t)+10*sind(20*t);
        input.dyn.Cl = zeros(size(t));
        input.dyn.Cd = zeros(size(t));

        out = functionLB(input,LBcoeffs,'general');
    otherwise
        error('Your test case have not been recognized. The three options are rampup, simcos, and general.')
end

figure, hold all
plot(input.dyn.t,out.CN)
plot(input.dyn.t,out.CN_LB)
