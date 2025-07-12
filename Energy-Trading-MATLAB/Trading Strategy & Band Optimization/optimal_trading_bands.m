function [R, grids] = optimal_trading_bands(M, l, f, P, C, alpha, return_grids)
%OPTIMAL_TRADING_BANDS  Optimal bands & CIs for an OU strategy.
%
%   INPUT ARGUMENTS
%   M            – # bootstrap replications (M = 1 -> point estimate only)
%   l            – Stop‑loss level (sigma‑units)
%   f            – Leverage (numeric)  or  'opt'  to optimise f*
%   P            – Struct with fields:
%                     .parameters       ← single‑fit theta,k,sigma,eta
%                     .bootstrapParams  ← vectors theta_i,k_i,sigma_i,eta_i (size M)
%   C            – Log‑scale transaction cost (not in sigma‑units)
%   alpha        – CI significance level (e.g. 0.05 -> 95 % two‑sided CI)
%   return_grids – Logical.  When TRUE, returns sigma‑grid arrays in *grids*
%
%   OUTPUT STRUCT R (one per l/f combination)
%   R.d_estimated        – Point estimate d* (sigma‑units, ‑ve magnitude)
%   R.u_estimated        – Point estimate u* (sigma‑units, +ve)
%   R.mu_estimated       – Point estimate mu* (log‑return / unit time)
%   R.d_CI, R.u_CI       – Two‑sided CIs for d*, u*
%   R.mu_CI              – CI for mu*
%   R.f_estimated        – Optimal f* (only if f == 'opt')
%   R.f_opt_CI           – CI for f*  (idem)
%   R.f_input            – Original f (numeric) OR NaN when optimised
%   R.profitable_sample  – Share of bootstrap draws where a trade exists
%
%   grids (optional)     – sigma‑grid arrays for diagnostic plotting
%                          .sigma  .d  .u  .f

    if nargin < 7
        return_grids = false; 
    end
    
    % Container for diagnostic grid
    grids = struct('sigma',[],'d',[],'u',[],'f',[],'label','');
    
    % Result struct – all numerical fields pre‑filled with NaN for safety
    R = struct( ...
        'd_estimated',  NaN, ...
        'u_estimated',  NaN, ...
        'mu_estimated', NaN, ...
        'd_CI',         [NaN NaN], ...
        'u_CI',         [NaN NaN], ...
        'mu_CI',        [NaN NaN], ...
        'f_estimated',  NaN, ...
        'f_opt_CI',     [NaN NaN], ...
        'f_input',      [], ...
        'profitable_sample', NaN ...
    );
    
    profit_perc = NaN;              % default when M == 1
    
    % Validate leverage argument
    
    if ~(isnumeric(f) || (ischar(f) && strcmpi(f,'opt')))
        error('f must be numeric or the string ''opt''.');
    end
    
    % Point estimates on full sample

    theta  = 1 / P.parameters.k;                                 % mean‑reversion scale
    sigma0 = P.parameters.sigma / sqrt(2 * P.parameters.k);      % stationary sigma
    
    [d,u,mu_no_theta,f_est] = sigma_vec_maximization(sigma0, C, l, f);
    mu = mu_no_theta / theta;                                    % long run return estimated
    
    % Bootstrap: uncertainty & discontinuity handling
    
    if M > 1
        % Prepare bootstrap vectors
    
        sigma_vec  = P.bootstrapParams.sigma ./ sqrt(2 * P.bootstrapParams.k);
        theta_vec  = 1 ./ P.bootstrapParams.k;
    
        % Build sigma‑grid for interpolation
    
        h = 1e-4;                                  
        n_grid = round((max(sigma_vec) - min(sigma_vec)) / h) + 1;
        sigma_grid = linspace(min(sigma_vec), max(sigma_vec), n_grid).';
    
        [d_grid, u_grid, mu_grid, f_grid] = sigma_vec_maximization(sigma_grid, C, l, f);
    
        if return_grids
            grids = struct('sigma', sigma_grid, 'd', d_grid, 'u', u_grid, ...
                            'f', f_grid, 'label', sprintf('f = %s',num2str(f)));
        end
    
        % Leverage jump (non interpolation samples)
           
        tol = 1e-12;
        
        if all(f_grid >= tol)
            % Entire sigma-grid is profitable -> Interpolate all
            sigma_noInt = [];
            theta_noInt = [];
            sigma_int   = sigma_vec;
            theta_int   = theta_vec;
        
        else
            % Jump detected -> find boundary for safe interpolation
            idx_jump = find(f_grid >= tol, 1, 'first');
        
            if isempty(idx_jump) || idx_jump <= 1
                % Either no profitable region or jump at beginning -> Interpolate all
                sigma_noInt = [];
                theta_noInt = [];
                sigma_int   = sigma_vec;
                theta_int   = theta_vec;
            else
                % Identify unsafe region near discontinuity
                jump_left  = sigma_grid(idx_jump - 1);
                jump_right = sigma_grid(idx_jump);
        
                sigma_mask   = (sigma_vec > jump_left) & (sigma_vec < jump_right);
                sigma_noInt  = sigma_vec(sigma_mask);
                theta_noInt  = theta_vec(sigma_mask);
                sigma_int    = sigma_vec(~sigma_mask);
                theta_int    = theta_vec(~sigma_mask);
            end
        end
    
        % Exact evaluation where interpolation would fail
    
        if ~isempty(sigma_noInt)
            [d_noInt,u_noInt,mu_noInt_noTheta,f_noInt] = sigma_vec_maximization(sigma_noInt,C,l,f);
            mu_noInt = mu_noInt_noTheta ./ theta_noInt;
        else
            d_noInt = []; u_noInt = []; mu_noInt = []; f_noInt = [];
        end
    
        % Interpolate elsewhere (if any remain)
    
        if ~isempty(sigma_int)
            d_boot  = interp1(sigma_grid,d_grid ,sigma_int);
            u_boot  = interp1(sigma_grid,u_grid ,sigma_int);
            mu_boot = interp1(sigma_grid,mu_grid,sigma_int) ./ theta_int;
        else
            d_boot = []; u_boot = []; mu_boot = [];
        end
    
        % Merge samples
    
        d_all  = [d_boot(:);  d_noInt(:)];
        u_all  = [u_boot(:);  u_noInt(:)];
        mu_all = [mu_boot(:); mu_noInt(:)];
    
        % Confidence intervals (filter NaNs)
    
        pct = [100*alpha/2, 100*(1-alpha/2)];
        R.d_CI  = prctile(d_all(~isnan(d_all)),  pct);
        R.u_CI  = prctile(u_all(~isnan(u_all)),  pct);
        R.mu_CI = prctile(mu_all(~isnan(mu_all)),pct);
    
        % Handle f* if optimised
    
        if ~isnumeric(f)
            if ~isempty(sigma_int)
                f_boot = interp1(sigma_grid,f_grid,sigma_int);
            else
                f_boot = [];
            end
            f_all = [f_boot(:); f_noInt(:)];
            R.f_opt_CI = prctile(f_all(~isnan(f_all)), pct);
        end
        
        profit_perc = sum(~isnan(d_all)) / M;
    
    end  
    
    % Populate estimation struct
    R.d_estimated         = d;
    R.u_estimated         = u;
    R.mu_estimated        = mu;
    R.profitable_sample   = profit_perc;
    
    if isnumeric(f)               % leverage fixed
        R.f_input      = f;
        R.f_estimated  = NaN;
        R.f_opt_CI     = [NaN NaN];
    else                          % leverage optimised
        R.f_input      = NaN;
        R.f_estimated  = f_est;
    end
end