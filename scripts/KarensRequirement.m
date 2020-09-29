

%% Static properties
static = load(fullfile('..','data','static_flatplate.mat'));

input.c = param.c;
input.static.alpha = static.alpha;
input.static.Cn = static.CN;
input.alpha_ss = 13;

%% Defining model parameters

% Leishman-Beddoes constants
LBcoeffs.Tp = 3;
LBcoeffs.Tf = 3;
LBcoeffs.Tv = 1;
LBcoeffs.Tvl = 1;

%% Loading dynamic data

% now testing for simcos data
test = 'general';
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
        n = length(data.Cl);
        TS = 1/LB(c).FS;
        input.dyn.t = 0:TS:(n-1)*TS;
        input.dyn.alpha = Alpha;
        input.dyn.Cl = Cl_corr;
        input.dyn.Cd = CD;

        % sinus parameters
        input.freq = LB(c).fosc;
        input.mean_rad = deg2rad(LB(c).alpha_0);
        input.amp_rad = deg2rad(LB(c).alpha_1);
        out = functionLB(input,LBcoeffs,'ramp');
    case 'rampup' % for rampup motion
        load(loadmat(LB(c).ms,LB(c).mpt),'raw','zero')
        input.U=LB(c).U;

        input.dyn.t = raw.t;
        input.dyn.alpha = raw.alpha;
        input.dyn.Cl = raw.Cl-mean(raw.Cl(1:50));
        input.dyn.Cd = raw.Cd-mean(raw.Cd(1:50));

        % enter here the pitch rate
        input.alphadot=0.2; % deg/s

        out = functionLB(input,LBcoeffs,'pitching');
    case 'general'
        input.U=LB(c).U;
        n = 100;
        Ts = 0.1;
        t = 0:Ts:(n-1)*Ts;
        input.dyn.t = t ;
        input.dyn.alpha = 10+20*sind(4*t)+10*sind(10*t)+10*sind(20*t);
        input.dyn.Cl = zeros(size(t));
        input.dyn.Cd = zeros(size(t));

        out = functionLB(input,LBcoeffs,'pitching');
    otherwise
        error('Your test case have not been recognized. The three options are rampup, simcos, and general.')
end

figure, hold all
plot(input.dyn.t,input.dyn.alpha)
plot(input.dyn.t,out.CN_LB)
