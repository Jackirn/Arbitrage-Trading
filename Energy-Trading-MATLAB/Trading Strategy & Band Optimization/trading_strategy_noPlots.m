function [optimal_bands, mu_OS, f_used_mat] = trading_strategy_noPlots(P, f_list, l_list, M, alpha, C, OS_data)

% TRADING_STRATEGY_NO_PLOT Evaluate optimal trading bands and compute OOS returns without plotting
%
% INPUTS and OUTPUTS same as original trading_strategy

    % Compute stationary parameters from OU process
    cap_sigma = P.parameters.sigma / sqrt(2 * P.parameters.k);
    cap_eta   = P.parameters.eta;

    N = numel(f_list); % Number of leverage levels
    L = numel(l_list); % Number of stop-loss levels

    mu_OS = zeros(N, L); % Matrix to store OS returns
    f_used_mat = zeros(N, L); % Leverage

    % Preallocate empty struct array
    template = struct( ...
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

    optimal_bands(N, L) = template;

    for j = 1:L

        l = l_list{j};

        fprintf('\n=============================================\n');
        fprintf(' Evaluating Optimal Trading Bands with Stop-Loss\n');
        fprintf('=============================================\n');
        fprintf('Using stop-loss level: l = %.3f (sigma units)\n\n', l);

        for i = 1:N

            f = f_list{i};

            % Compute optimal bands (plot flag set false)
            [R, ~] = optimal_trading_bands(M, l, f, P, C, alpha, false);
            optimal_bands(i,j) = R;

            % Determine leverage used
            if isequal(f, 'opt')
                f_used = optimal_bands(i,j).f_estimated;
            else
                f_used = f;
            end
            
            % Store leverage used
            f_used_mat(i,j) = f_used;
            % Compute OS return
            mu_OS(i,j) = OS_return(optimal_bands(i,j).u_estimated, ...
                                   optimal_bands(i,j).d_estimated, ...
                                   l, f_used, OS_data, C, cap_sigma, cap_eta, []);

            % Print results
            fprintf('------------------------------------------------\n');
            fprintf(' Results for f = %s (used f = %.4f)\n', num2str(f), f_used);
            fprintf('------------------------------------------------\n');

            if ~isnan(optimal_bands(i,j).u_estimated) && ~isnan(optimal_bands(i,j).d_estimated)
                fprintf('Estimated u*   : %.3f\n', optimal_bands(i,j).u_estimated);
                fprintf('Estimated d*   : %.3f\n', -optimal_bands(i,j).d_estimated);
            else
                fprintf(' No estimated bands since f optimal estimated is equal to 0\n\n')
            end

            fprintf('Estimated mu*  : %.3f\n\n', optimal_bands(i,j).mu_estimated);

            if ~isnan(optimal_bands(i,j).d_CI)
                print_CI_intervals(optimal_bands(i,j).u_CI, ...
                                   optimal_bands(i,j).d_CI, ...
                                   optimal_bands(i,j).mu_CI);
            else
                fprintf('No confidence interval founded since there are no profitable samples\n\n')
            end

            if isequal(f, 'opt') && optimal_bands(i,j).f_opt_CI(1) ~= 0.0 && optimal_bands(i,j).f_opt_CI(2) ~= 0.0
                fprintf('------------------------------\n');
                fprintf('  95%% CI for optimal f\n');
                fprintf('  f* : [%.4f , %.4f]\n', ...
                         optimal_bands(i,j).f_opt_CI(1), ...
                         optimal_bands(i,j).f_opt_CI(2));
                fprintf('------------------------------\n');
            end

            if optimal_bands(i,j).profitable_sample < 1
                fprintf('Profitable sample for CI computation: %.3f%%\n\n', ...
                    optimal_bands(i,j).profitable_sample * 100);
            end

            if ~isnan(mu_OS(i,j))
                fprintf('Out-of-sample return: %.3f\n\n', mu_OS(i,j));
            end
        end
    end

    % Convert results into tables for readability
    rowNames = cellfun(@num2str, f_list, 'UniformOutput', false);
    colNames = cellfun(@num2str, l_list, 'UniformOutput', false);

    optimal_bands = cell2table(mat2cell(optimal_bands, ones(1,N), ones(1,L)), ...
                                'RowNames', rowNames, 'VariableNames', colNames);

    mu_OS = array2table(mu_OS, 'RowNames', rowNames, 'VariableNames', colNames);

end
