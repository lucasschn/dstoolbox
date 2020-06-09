function x_filtered = myFilterTwice(x,fs)
% Passes myFilter twice in both direction to cancel induced
% phase lag, following the procedure from Keneth et al. 2011. 
xf = myFilter(x,fs);
x_filtered = flip(myFilter(flip(xf),fs));
end