function xfff = myFilter(x,fs)
% This function applies a Butterworth filter, a moving average and a
% Chebychev type-II filter to the input signal x sampled at sampling
% frequency fs, following the procedure from Keneth et al. 2011.

% Butterworth filter
fc = 35;
[b,a] = butter(5,fc/(fs/2));
xfiltered = filter(b,a,x);

% Moving average filter
%xff = movmean(xfiltered,30);
xff = xfiltered; % moving average deactivated

% Chebychev type-II filter
fp = 1/3;
[b,a] = cheby2(6,20,36*fp/(fs/2));
xfff = filter(b,a,xff);

end
