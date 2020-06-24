function x_filtered = myFilterTwice(x,fs,domovmean)
% Passes myFilter twice in both direction to cancel induced
% phase lag, following the procedure from Keneth et al. 2011.
if exist('domovmean','var')
    xf = myFilter(x,fs,domovmean);
    x_filtered = flip(myFilter(flip(xf),fs,domovmean));
else
    xf = myFilter(x,fs);
    x_filtered = flip(myFilter(flip(xf),fs));
end
end