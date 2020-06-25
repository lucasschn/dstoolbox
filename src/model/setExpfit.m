function  fig = setExpfit(airfoil,varargin)
% In Sheng 2006, the curve fitting alpha_ds as a function of r is
% linear. Here it is exponential.


% abreviations
alpha_ss = airfoil.steady.alpha_ss;


% vector initialization
r = -ones(size(varargin));
alpha_ds = -ones(size(varargin));
t_ss = -ones(size(varargin));
t_ds = -ones(size(varargin));

for k=1:length(varargin)
    ramp = varargin{k};
    r(k) = ramp.r; 
    t_ss(k) = ramp.t(find(ramp.alpha>alpha_ss,1));
    if ramp.alpha_CConset > alpha_ss
        alpha_ds(k) = ramp.alpha_CConset;
        t_ds(k) = ramp.t(ramp.i_CConset);
    elseif ramp.alpha_CLonset > alpha_ss
        warning('DS was based on CL for r=%.3f',ramp.r)
        alpha_ds(k) = ramp.alpha_CLonset;
        t_ds(k) = ramp.t(ramp.i_CLonset);
    else
        error('The dynamic stall angle is lower than the static stall angle.')
    end
end

alpha_ds_r = @(x,r) x(1)-(x(1)-alpha_ss)*exp(-x(2)*r);

% compute the exponential fit
opts = optimset('Diagnostics','off','Display','off');
xopt = lsqcurvefit(alpha_ds_r,[alpha_ds(end) 1],r,alpha_ds,[0 0],[Inf Inf],opts);

A = xopt(1);
B = xopt(2);
save('../expfit_flatplate','A','B')

% compute the fitting error

fit_error = sum((alpha_ds - alpha_ds_r([A,B],r)).^2)/length(alpha_ds);
fprintf('The exponential fitting error is %.2f. \n',fit_error)

% plot alpha_ds(r) and exponential fitting

fig = figure;
subplot(311)
plot(r,alpha_ds,'d','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','\alpha_{ds} (exp)')
grid on
ylabel('\alpha_{ds} (°)','FontSize',20);
axis([0 0.12 0 30]);
ax = gca;
ax.FontSize = 20;
hold on
plot(r,alpha_ds_r(xopt,r),'LineWidth',2,'DisplayName','exponential fit')
subplot(313)
plot(r,t_ds-t_ss,'.','MarkerSize',20)
xlabel('r')
ylabel('\Delta t_{ds} (s)')
grid on
ax = gca;
ax.FontSize = 20;
end