%% PROCESS AND PLOT POLAR
%%% Use polars to calibrate and/or check that your setup is ready 
%%% to rock n'roll!
% Date:     22/04/20
% File:     polarplotter.m
% By:       Sébastien Le Fouest

clear rawpath -date

% define you calipath if it doesn't exist
 if ~exist ('rawpath','var')
    rawpath = '../data/2020_SH/';
%     rawpath = '\\sti1raw.epfl.ch\unfold\smartH\2020_SLF\raw';
 date = 20200510; % set date for saving purposes in format yyyymmdd
 end
 
calibrationfolder = fullfile(rawpath,num2str(date),'calibration');

polarfolder = fullfile(rawpath,num2str(date),'loads');

%% INPUTS

% load your a_sweep data to plot the polar, specify whether you wish to
% determine alpha_0 from a sweep around your manual zeroing

% options
find_alpha0 = false; % do you wish to find alpha_0?
auto_home = false; % do you wish to use the result of the interpolation to go home?
saveplot = false; % do you wish to save a pdf of your figure?
savepolar = true; % do you wish to save the data to plot your polar?
if savepolar == 1
polarname = 'static_flatplate5'; % what is the file name for your polar data?
end

% name of files
polardata = 'ms121mpt005'; % data for the polar you wish to plot
forceoffsetdata = 'calibration_mpt'; % data with your non-hydrodynamic 
% force offset (LC offset, inertial forces, buoyancy etc.)

% wing data
param.chord = 0.15;         % chorch length in m    
param.span = 0.6;           % span in m

% fluid data
param.Uinf = 0.5;          % free stream velocity in m/s
param.rho = 1000;           % density in kg/m^3
param.mu =  10^-3;          % dynamic viscosity in N s/m^2            
param.nu = param.mu/param.rho;       % kinematic viscosity in m^2/s

% adim numbers %
param.Re = param.Uinf*param.chord/param.nu;
param.tc = param.chord./param.Uinf; % convective time

%LC angle 
LCangle = 90;

%% PROCESS ZERO FORCES

if ~exist('force_offset.mat','file')
load(fullfile(calibrationfolder,forceoffsetdata))

    %initialise forces
    n = length([out.alpha]);
    force_offset.alpha = [out.alpha];
    force_offset.F1 = NaN(n,1);
    force_offset.F2 = NaN(n,1);
    force_offset.FZ = NaN(n,1);
    force_offset.M1 = NaN(n,1);
    force_offset.M2 = NaN(n,1);
    force_offset.MZ = NaN(n,1);

    for i = 1:n
        force_offset.F1(i) = mean(out(i).forces(:,1));
        force_offset.F2(i) = mean(out(i).forces(:,2));
        force_offset.FZ(i) = mean(out(i).forces(:,3));
        force_offset.M1(i) = mean(out(i).forces(:,4));
        force_offset.M2(i) = mean(out(i).forces(:,5));
        force_offset.MZ(i) = mean(out(i).forces(:,6));
    end
    clear out
    save(fullfile(calibrationfolder,'force_offset'),'force_offset')
else 
    load ('force_offset.mat')
end
%% PROCESS POLAR
load(fullfile(polarfolder,polardata))
qA = 0.5.*param.rho.*(param.Uinf).^2.*param.chord.*param.span;
F1_0 = interp1([force_offset.alpha],[force_offset.F1],[out.alpha],'linear','extrap');
F2_0 = interp1([force_offset.alpha],[force_offset.F2],[out.alpha],'linear','extrap');
MZ_0 = interp1([force_offset.alpha],[force_offset.MZ],[out.alpha],'linear','extrap');
alpha = [out.alpha];
n = length(alpha);

F1 = NaN(size(alpha));
F2 = NaN(size(alpha));
MZ = NaN(size(alpha));

for i = 1:n
    F1(i) = mean(out(i).forces(:,1))-F1_0(i);
    F2(i) = mean(out(i).forces(:,2))-F2_0(i);
    MZ(i) = mean(out(i).forces(:,6))-MZ_0(i);
end
FX = (sind(LCangle).*F1 - cosd(LCangle).*F2);
FY = -(cosd(LCangle).*F1 + sind(LCangle).*F2);

Lift = - FX.*sind(alpha) + FY.*cosd(alpha);
Drag = FX.*cosd(alpha) + FY.*sind(alpha);
CLift = Lift./qA;
CDrag = Drag./qA;
CNormal = FY./qA;
if savepolar == 1
    polarforces.FX = FX;
    polarforces.FY = FY;
    polarforces.MZ = MZ;
    polarforces.Lift = Lift;
    polarforces.Drag = Drag;
    polarforces.CL = CLift;
    polarforces.CD = CDrag;
    polarforces.CN = CNormal;
    save(fullfile(calibrationfolder,polarname),'alpha')
    save(fullfile(calibrationfolder,polarname),'-struct','polarforces','CN','CL','CD','-append')
end

%% PLOT POLAR

scnsize = get(0,'screensize'); %To get nicer figures
figFullsize = [5 50 scnsize(3)*0.7 scnsize(4)*0.85]; % to get figures on full screen
whiteBckgnd = [1 1 1];

name = ('Lift and drag polar');
f2 = figure('name',name,'position', figFullsize , 'color', whiteBckgnd); %White Background
hax=axes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT 1 %%%%%%%%%%%%%%%%%%%%%%%
ax(1) = subplot(2,1,1);

hold on
ToPlot1 = 'CLift';
ToPlot2 = 'CDrag';

a = plot(alpha,CLift,'r','LineWidth',1.5,'LineStyle','-','Marker','+'); % Plot lift coef obtained from out own calculations
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
alpha2pi = -5:5;
a= plot(alpha2pi,2*pi*deg2rad(alpha2pi),'k--');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([-1 2])
xlabel('angle of attack') 
ylabel(ToPlot1)

% hl = legend('Location','SouthEast');%(p,DisplayName)
% set(hl, 'Interpreter','latex')
grid on
grid minor
% xlim([alpha(1) alpha(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT 2 %%%%%%%%%%%%%%%%%%%%%%%
ax(2) = subplot(2,1,2);

hold on
a = plot(alpha,CDrag,'r','LineWidth',1.5,'LineStyle','-','Marker','+');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
% hl = legend('Location','SouthEast');%(p,DisplayName)
% set(hl, 'Interpreter','latex')
grid on
grid minor
xlabel('angle of attack') 
ylabel(ToPlot2)
linkaxes(ax,'x')
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
if saveplot == 1
% print(gcf, '-dpdf', fullfile(calibrationfolder,'CL_CD_vsAlpha.pdf'),'-r0')
print(gcf, '-dpdf', 'CL_CD_vsAlpha.pdf','-r0')
end

%% FIND ALPHA_0

if find_alpha0
   alpha_0 = interp1(Lift,alpha,0);
   f = msgbox(sprintf('The angle offset is %.2f°', alpha_0));
    
   if auto_home
       g = ConnectGalil('192.168.255.201');
        if ~exist('flag_init','var')
            flag_init = false;
        end
        
        if flag_init == false %set everything once. Don't need to do it everytime
            %%%%%%%%%% Motor setup %%%%%%%%%%%%%%
            % for stepper motors decide :
            % NEMA17
            g.command('MO');    % Motor off
            g.command('MT-2.5');% set motor type -2.5 (reversed direction cw)
            g.command('LC2'); % when nothing happens, low current applied in motor => avoid overheating
            m(1).ms = 1280000/360; % Number of counts for 1 degree NEMA17 + Gear 1:25 (!) 1280000 counts for 360degrees
            m(1).n = 'A';
            m(1).es = 12500/9;
            m(1).RCN = 10; % Encoder sampling rate parameter (figuring out how it works, usually was at 10)

            flag_init = true;
        end
        
        g.command(['SH' AllMotNam(m)]); % Turn on inactive motor
        g.command('RC0'); % close eventual recording not aborted
        pause(0.1)

        % Get actual position
        actual_pos = str2double(g.command('TP'))./m.es; %[deg]
        curr_pos_motor = num2str(round(actual_pos .* m.ms)); % in motor step
        % Set position using encoder data
            pause(0.2)
            g.command(['DP' curr_pos_motor]);

        %Get the position of the angle where we want to go 
        goal_pos_motor = num2str(round(alpha_0 .* m.ms));

        % Get an absolute position from a position set by the number of step
        pause(0.2)
        g.command(['PA' goal_pos_motor]);
        g.command('BG'); %Turn little cute motor!

        while str2double(g.command('MG _BGA'))
                pause(10^-2) % wait the rotation is finished
        end
        
        g.command('DP0');   % set motor position 0
        g.command('DE0');   % set encoder position 
        f = msgbox(sprintf('The blade was brought to alpha_0, try to record data and check whether lift values are approximately 0.'));
    end 
end