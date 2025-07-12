function S = ouBootstrap(x, dt, M, alpha)
%OUBOOTSTRAP  Perform parametric bootstrap for OU parameter estimation.
%
%   S = OUBOOTSTRAP(x, dt, M, alpha)
%
%   Simulates M bootstrap paths of an Ornstein-Uhlenbeck (OU) process using
%   the estimated parameters from the input series x, and returns bootstrap
%   distributions, confidence intervals, and summaries.
%
%   INPUTS:
%     x      - Time series (vector) of observed values.
%     dt     - Time step between observations.
%     M      - Number of bootstrap samples to generate.
%     alpha  - Significance level (e.g., 0.05 for 95% confidence intervals).
%
%   OUTPUT:
%     S - Struct with fields:
%          parameters      : table with estimated (k, eta, sigma) from x.
%          bootstrapParams : M×3 table with bootstrap estimates.
%          CI              : 2×3 table with lower and upper percentile bounds.
%          median          : 1×3 table with median bootstrap values.
%
%   NOTES:
%     - Uses functions: ou_mle, ou_sim.
%     - CI is computed using empirical percentiles from bootstrap distribution.

    N = length(x) - 1;               % Number of time steps for simulation
    [k, eta, sigma] = ou_mle(x, dt); % Estimate OU parameters via MLE on original data
    
    P = zeros(M,3); % Preallocate matrix for bootstrap estimates
    for m = 1:M
        xSim = ou_sim(x(1), k, eta, sigma, dt, N); % Simulate OU path
        [k_m, eta_m, sigma_m] = ou_mle(xSim, dt);  % Estimate MLE on simulated path
        P(m,:) = [k_m, eta_m, sigma_m];            % Store bootstrap estimates
    end
    
    colNames = {'k', 'eta', 'sigma'}; % Parameter names
    T_hat    = array2table([k, eta, sigma], 'VariableNames', colNames); % Original parameter estimates
    T_sim    = array2table(P,               'VariableNames', colNames); % Bootstrap parameter estimates
    
    CI_mat = prctile(P, [alpha*100/2, 100*(1 - alpha/2)]); % Compute bootstrap confidence intervals
    CI_tab = array2table(CI_mat, ...
               'RowNames', {'Lower', 'Upper'}, ...
               'VariableNames', colNames); % Confidence interval table
    
    medTab = array2table(median(P), 'VariableNames', colNames); % Bootstrap medians
    
    % Final output struct
    S.parameters      = T_hat;  % MLE estimates from original data
    S.bootstrapParams = T_sim;  % Bootstrap samples
    S.CI              = CI_tab; % Bootstrap confidence intervals
    S.median          = medTab; % Bootstrap medians
end