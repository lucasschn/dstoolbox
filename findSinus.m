function [mean,amp,f_pts,phi]=findSinus(alpha,plotanalpha)
% findSinus(alpha) returns the mid-point mean and the amplitude of
% the oscillation amp_rad of a sinusoidal given signal alpha. These
% three variables can be in radians or in degrees. Additionally, findSinus
% returns also the oscillation frequency f_pts in number of turns per sampling interval and the phase-lag phi of the
% main harmonic.

N = length(alpha);

% factor to lower the sampling resolution in order to improve computing
% time. If m=1, the input data is unchanged.
m = 1;
alpha = alpha(1:m:N);
N = N/m;

fft_coeffs = fft(alpha,N); % freq : [freq 0, N/2-1 freq pos, N/2 freq neg]
P2 = abs(fft_coeffs/N); % two-fold power spectrum
P1 = P2(1:N/2-1);
P1(2:end-1) = 2*P1(2:end-1); % one-fold power spectrum
f = linspace(0,0.5,N/2-1); % frequency vector for one-fold power spectrum in rev/sampling interval

[pks, locs] = maxk(P1(f>0),1);

mean = P1(f==0);
amp = pks(1);
f_pts = f(locs(1)+1);
z = fft_coeffs(locs(1)+1);
phi = atan(real(z)/imag(z));

if nargin>1 && plotanalpha
    analpha = mean + amp*sin(2*pi*f_pts*(0:N)+phi);
    figure
    plot(alpha,'DisplayName','exp')
    hold on
    plot(analpha,'DisplayName','ideal')
    grid on
    legend show
end
end