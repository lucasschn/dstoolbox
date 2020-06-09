function [alpha_lag_ds,Talpha] = Expfit(airfoil,ramp)
% Convenience routine that computes the dynamic stall angle,
% and fits an exponential fitting to a bunch of ramp-up motions. It
% then computes the delay to be applied and computes the delayed AoA alpha_lag for each of
% them. The model specific critical pitch angle, alpha_crit, is
% then computed and all the created data is ploted on a reduced
% pitch rate vs pitch angle graph.
% initialization of the vectors
debug = false;
% abreviations
if isempty(ramp.r)
    error('Reduced pitch rate is not determined.')
elseif isempty(airfoil.steady.alpha_ss)
    error('Static stall angle is not determined.')
else
    r = ramp.r;
    alpha_ss = airfoil.steady.alpha_ss;
end

% define shape of exponential fit
alpha_ds_r = @(x,r) x(1)-(x(1)-alpha_ss)*exp(-x(2)*r);

% determine dynamic stall for this 
t0 = interp1(ramp.analpha,ramp.t,0);
if ~isempty(ramp.i_CConset)
    if ramp.alpha_CConset > alpha_ss
        t_ds = ramp.t(ramp.i_CConset)-t0; % time since start of ramp
        alpha_ds = ramp.alpha_CConset;
    elseif ramp.alpha_CLonset > alpha_ss
        warning('CL was used to define DS for r=%.3f',ramp.r)
        t_ds = ramp.t(ramp.i_CLonset)-t0;
        alpha_ds = ramp.alpha_CLonset;
    else
        error('The dynamic stall angle is lower than the static stall angle.')
    end
    if debug
        fprintf('From experiment : alpha_ds = %.2f°\n',alpha_ds)
        fprintf('From experiment : t_ds = %.2fs\n',t_ds)
    end
else % rely on the expfit
    load(sprintf('/Users/lucas/Documents/EPFL/PDM/expfit_%s',airfoil.name),'A','B')
    fprintf('A = %.2f, B = %.2f \n',A,B)
    alpha_ds = alpha_ds_r([A B],r);
    t_ds = alpha_ds/ramp.alphadot; % time since start of ramp
    if debug
        fprintf('From expfit : alpha_ds = %.2f° \n',alpha_ds)
        fprintf('From expfit : t_ds = %.2fs \n',t_ds)
    end
end

if isempty(alpha_ss)
    error('Static stall angle is not determined.')
elseif isempty(ramp.alphadot)
    error('Pitch rate is not determined.')
elseif isempty(t_ds)
    error('DS time is not determined.')
elseif isempty(ramp.V)
    error('Inflow velocity is not determined.')
else
% determination of Talpha so that alpha_lag(t_ds) = alpha_ss
syms tau % dimensional time constant
sol = solve(alpha_ss == ramp.alphadot*(t_ds - tau*(1-exp(-t_ds/tau))),tau,'Real',true,'IgnoreAnalyticConstraints',true);
Talpha = 2*ramp.V/airfoil.c*double(sol); % in adimensional time here
end

% compute alpha_lag using Talpha and finds alpha_lagonset
if ~isempty(Talpha)
    ramp.computeAlphaLag(airfoil,Talpha)
else 
    error('Talpha could not be found to satisfy the condition alpha_lag_ds=alpha_ss')
end

% looking for the value of alpha_lag(alpha) at the point alpha_ds
if isempty(ramp.alpha)
    alpha_lag_ds = interp1(ramp.analpha,ramp.analpha_lag,alpha_ds);
elseif isempty(ramp.i_continuous_grow)
    alpha_lag_ds = interp1(ramp.alpha,ramp.alpha_lag,alpha_ds);
else % if alpha_continuous_grow is defined
    alpha_lag_ds = interp1(ramp.alpha_continuous_grow,ramp.alpha_lag(ramp.i_continuous_grow),alpha_ds);
end

end