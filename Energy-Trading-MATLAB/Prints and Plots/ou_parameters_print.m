function ou_parameters_print(P)
%OU_PARAMETERS_PRINT  Display MLE estimates and bootstrap confidence intervals.
%
%   ou_parameters_print(P)
%
%   Prints to the Command Window the estimated parameters of an
%   Ornstein-Uhlenbeck (OU) process, along with their 95% bootstrap
%   confidence intervals.
%
%   INPUT:
%     P - Struct with fields:
%           P.parameters : table with fields 'k', 'eta', 'sigma'
%           P.CI         : table with 2 rows ('Lower', 'Upper') and 
%                            columns 'k', 'eta', 'sigma'

    fprintf('------------------------------\n');
    fprintf('     MLE Parameter Estimates\n');
    fprintf('------------------------------\n');
    fprintf('  kappa_hat : %12.4f\n', P.parameters.k);
    fprintf('    eta_hat : %11.4f%%\n', 100*P.parameters.eta);
    fprintf('  sigma_hat : %11.2f%%\n', 100 * P.parameters.sigma);
    fprintf('------------------------------\n');

    fprintf('---------------------------------------------\n');
    fprintf('     Bootstrap Confidence Intervals (95%%)\n');
    fprintf('---------------------------------------------\n');
    fprintf('     kappa : [%10.5f , %10.5f]\n', P.CI.k(1), P.CI.k(2));
    fprintf('       eta : [%10.5f%% , %10.5f%%]\n', 100*P.CI.eta(1),   100*P.CI.eta(2));
    fprintf('     sigma : [%10.5f%% , %10.5f%%]\n', 100*P.CI.sigma(1), 100*P.CI.sigma(2));
    fprintf('---------------------------------------------\n');
end