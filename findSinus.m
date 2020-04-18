function [mean,amp,omega,phi]=findSinus(alpha,plotanalpha)
% findSinus(alpha) returns the mid-point mean and the amplitude of
% the oscillation amp_rad of a sinusoidal given signal alpha. These
% three variables can be in radians or in degrees.
m = 1;

N = length(alpha);
alpha = alpha(1:m:N);
N = N/m;

fft_coeffs = fft(alpha,N); % freq : [freq 0, N/2-1 freq pos, N/2 freq neg]
P2 = abs(fft_coeffs/N); % two-fold power spectrum
P1 = P2(1:N/2-1);
P1(2:end-1) = 2*P1(2:end-1); % one-fold power spectrum
f = linspace(0,0.5,N/2-1); % frequency vector for one-fold power spectrum in rev/sampling interval

[pks, locs] = maxk(P1(f>0),1);

alpha0 = P1(f==0);
alpha1 = pks(1);

z = fft_coeffs(locs(1)+1);
phi = atan(real(z)/imag(z));

if nargin>1 && plotanalpha
    analpha = alpha0 + alpha1*sin(2*pi*f(locs(1)+1)*(0:N)+phi);
    figure
    plot(alpha,'DisplayName','exp')
    hold on
    plot(analpha,'DisplayName','ideal')
    grid on
    legend show
end
end