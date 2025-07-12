function results = pair_trading_analysis_noPlots(T, splitMonths, M, alpha, l_list, f_list)
    % Step 1: Split dataset into IS and OS
    [data_pre_out, data_OS_pre_out] = trimAndSplitPriceTable(T, [], [], [], [], splitMonths);

    % Step 2: Remove outliers
    [data_IS, outliers_IS, isOutlierMask_IS] = removeOutliersIS(data_pre_out);
    [data_OS, outliers_OS, isOutlierMask_OS] = removeOutliersIS(data_OS_pre_out);

    % Step 3: (Plot calls removed)

    % Step 4: Calculate dt
    N = height(data_IS);
    t = yearfrac(data_IS.Time(1), data_IS.Time(end), 1);
    dt = t / N;

    % Step 5: Bootstrap OU parameters
    P = ouBootstrap(data_IS.Rt, dt, M, alpha);

    % Step 6: Filter negative k samples
    k_positive_mask = P.bootstrapParams.k > 0;
    num_removed = sum(~k_positive_mask);
    if num_removed > 0
        P.bootstrapParams = P.bootstrapParams(k_positive_mask, :);
        fprintf("With a percentage of negative k removed from the bootstrap of %.4f%%\n", ...
            num_removed / M * 100);
        fprintf("Total amount of sample removed is %d\n", num_removed);
    end

    % Step 7: Print OU parameter estimates
    ou_parameters_print(P);

    % Step 8: Compute transaction costs
    C = transaction_cost(data_IS);

    % Step 9: Run trading strategy optimization and evaluation on OS data
    [optimal_bands, mu_OS, f_used_mat] = trading_strategy_noPlots(P, f_list, l_list, M, alpha, C, data_OS);

    % Step 10: (Plot call removed)

    % Collect results
    results = struct();
    results.data_IS = data_IS;
    results.data_OS = data_OS;
    results.outliers_IS = outliers_IS;
    results.outliers_OS = outliers_OS;
    results.bootstrapParams = P.bootstrapParams;
    results.OU_parameters = P.parameters;
    results.CI = P.CI;
    results.optimal_bands = optimal_bands;
    results.mu_OS = mu_OS;
    results.C = C;
    results.f_used = f_used_mat;
end
