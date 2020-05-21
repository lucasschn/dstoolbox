function Sheng(ramp,airfoil)
% Loads linfit variables from a previous execution of setLifit.

load(sprintf('linfit_%s.mat',airfoil.name));
ramp.computeAlphaLag(airfoil,Talpha),

alpha_ss = airfoil.steady.alpha_ss;
alpha_ds = ramp.alpha_CConset;

% definition of alpha_crit depending on r (only for Sheng)
ramp.findModelOnset(obj); % alpha_lagonset = alpha_lag_ds only if Talpha is correct
if ramp.r >= r0
    alpha_crit = alpha_ds0;
else % if r<r0
    alpha_crit = alpha_ss + (alpha_ds0 - alpha_ss)*ramp.r/r0;
end


if isempty(ramp.alpha)
    alpha_lag_ds = interp1(ramp.analpha,ramp.analpha_lag,alpha_ds);
elseif isempty(ramp.i_continuous_grow)
    alpha_lag_ds = interp1(ramp.alpha,ramp.alpha_lag,alpha_ds);
else % if alpha_continuous_grow is defined
    alpha_lag_ds = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),alpha_ds);
end

fprintf('alpha_lag_ds is %.2f while alpha_crit is %.2f',alpha_lag_ds,alpha_crit)
end
