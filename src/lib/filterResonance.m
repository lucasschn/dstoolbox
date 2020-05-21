function rampout = filterResonance(ramp)
N = length(ramp.CL);

% factor to lower the sampling resolution in order to improve computing
% time. If m=1, the input data is unchanged.
m = 1;
CL = ramp.CL(1:m:N);
N = N/m;

% if N is even
% freq : [freq 0, N/2-1 freq>0, N/2 freq<0]
% if N is odd 
% freq : [freq 0, N/2 freq>0, N/2 freq<0]

fft_coeffs = fft(CL,N); 
P2 = abs(fft_coeffs/N); % two-fold power spectrum

if mod(N,2)==0
    P1 = P2(1:N/2-1);
else 
    P1 = P2(1:N/2);
end

P1(2:end) = 2*P1(2:end); % one-fold power spectrum
f = reshape(linspace(0,0.5,(N-mod(N,2))/2)/ramp.Ts,size(P1)); % frequency vector for one-fold power spectrum in rev/sampling interval

K = lsqcurvefit(@(x,xdata) x./xdata,0.06,f(2:200),P1(2:200));

filtered_CL = movmean(ramp.CL,4);


figure
plot(f,P1,'DisplayName','spectrum')
hold on
plot(f,0.06./f,'DisplayName','init')
plot(f,K./f,'DisplayName','opt')
grid on 
legend('Location','NorthEast')

% filtered_CL = ifft(0.06./f);

figure
plot(ramp.t,ramp.CL)
hold on 
plot(ramp.t,filtered_CL)
grid on 
xlabel('t (s)')
ylabel('C_L')

end