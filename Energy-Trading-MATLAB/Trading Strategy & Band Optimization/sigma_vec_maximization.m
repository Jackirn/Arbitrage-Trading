function [d_grid, u_grid, mu_grid, f_grid] = sigma_vec_maximization(sigma_vec, C, l, f)
%SIGMA_VEC_MAXIMIZATION  Grid evaluation of optimal bands for vector of σ
%
% INPUTS
%   sigma_vec – Vector of stationary σ values (scalar or array)
%   C         – Transaction cost (log scale)
%   l         – Stop-loss level (negative, in σ-units)
%   f         – Leverage (numeric or 'opt')
%
% OUTPUTS
%   d_grid    – Lower band (abs value, σ-units, negative)
%   u_grid    – Upper band (σ-units, positive)
%   mu_grid   – Long-run log-return (in-sample)
%   f_grid    – Optimal leverage (only if f = 'opt')

len = length(sigma_vec);
d_grid  = zeros(len,1);
u_grid  = zeros(len,1);
mu_grid = zeros(len,1);
f_grid  = zeros(len,1);

for k = 1:len
    sig_k = sigma_vec(k);              % current σ
    c_k   = C / sig_k;                 % scaled transaction cost
    lb_k  = [l + 0.01, l + c_k];       % lower bounds for [d u]
    ub_k  = [0.6, 3];   % upper bounds for [d u]

    x0 = [-0.5 ,0.5];

    options = optimoptions('fmincon',...
        'Display','off','Algorithm','interior-point',...
        'MaxIterations', 500, 'OptimalityTolerance', 1e-8);

    obj = @(x) -long_return(x(1), x(2), c_k, l, sig_k, f);  % objective (negated)
    
    % Initial guess
    val = obj(x0);
    if ~isfinite(val) || ~isreal(val)
        x0 = (lb_k + ub_k) / 2;        % fallback: midpoint of bounds
        val = obj(x0);
        if ~isfinite(val) || ~isreal(val)
            error('Objective function still invalid at adjusted initial guess.');
        end
    end

    % Constrained optimization
    [x_k, negmu_k, flg] = fmincon(obj, x0, [], [], [], [], ...
                                  lb_k, ub_k, @(x) constraints(x, c_k, l), options);

    if flg > 0
        d_grid(k)  = abs(x_k(1));
        u_grid(k)  =      x_k(2);
        mu_grid(k) = -negmu_k;

        % Compute optimal leverage if requested
        if ~isnumeric(f)
            [~, f_grid(k)] = long_return(-d_grid(k), u_grid(k), c_k, l, sig_k, f);
            
            % Handle cases where optimal f = 0 (no trade zone)
            if f_grid(k) == 0
                u_grid(k) = NaN;
                d_grid(k) = NaN;
            end
        end
    else
        d_grid(k) = NaN;
        u_grid(k) = NaN;
        mu_grid(k) = NaN;
    end
end
end
