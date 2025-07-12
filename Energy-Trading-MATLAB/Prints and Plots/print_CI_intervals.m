function print_CI_intervals(CI_u, CI_d, CI_mu)
%PRINT_CI_INTERVALS Stampa a schermo gli intervalli di confidenza al 95%
%
%   INPUT:
%     CI_u  - 1x2 vector [lower_u, upper_u]
%     CI_d  - 1x2 vector [lower_d, upper_d]
%     CI_mu - 1x2 vector [lower_mu, upper_mu]

    fprintf('------------------------------\n');
    fprintf('     95%% Confidence Intervals\n');
    fprintf('------------------------------\n');
    fprintf('  u*   : [ %8.4f , %8.4f ]\n', CI_u(1),  CI_u(2));
    fprintf('  d*   : [ %8.4f , %8.4f ]\n', -CI_d(2),  -CI_d(1));
    fprintf('  mu*  : [ %6.2f%% , %6.2f%% ]\n', 100*CI_mu(1), 100*CI_mu(2));
    fprintf('------------------------------\n');
end
