function results = pair_trading_analysis(T, splitMonths, M, alpha, l_list, f_list, f_up)
%PAIR_TRADING_ANALYSIS Perform pair trading analysis with bootstrap OU estimation
%
% Inputs:
%   T           - Table containing price/time data
%   splitMonths - Month index to split data into In-Sample (IS) and Out-of-Sample (OS)
%   M           - Number of bootstrap samples for parameter estimation
%   alpha       - Significance level for confidence intervals (e.g., 0.05 for 95%)
%   l_list      - List of stop-loss thresholds for the trading strategy
%   f_list      - List of leverage factors or optimization modes
%
% Steps:
%   1. Trim and split the dataset into IS and OS based on splitMonths
%   2. Remove outliers from both IS and OS datasets
%   3. Plot and print identified outliers
%   4. Calculate time step dt for the OU parameter estimation
%   5. Bootstrap Ornstein-Uhlenbeck (OU) process parameters from IS returns
%   6. Filter bootstrap samples to remove any negative mean-reversion rates k
%   7. Print summary of OU parameter estimation and bootstrap filtering
%   8. Compute transaction costs based on the dataset
%   9. Run trading strategy optimization and evaluation on OS data
%   10. Plot return vs leverage for various stop-loss thresholds

    % Step 1: Split the dataset into In-Sample (IS) and Out-of-Sample (OS) segments
    [data_pre_out, data_OS_pre_out] = trimAndSplitPriceTable(T, [], [], [], [], splitMonths);

    % Step 2: Remove outliers from IS and OS datasets
    [data_IS, outliers_IS, isOutlierMask_IS] = removeOutliersIS(data_pre_out);
    [data_OS, outliers_OS, isOutlierMask_OS] = removeOutliersIS(data_OS_pre_out);

    % Step 3: Visualize and print outliers for both IS and OS
    plotandprintOutliers(data_pre_out, isOutlierMask_IS, 'Outliers from IS Dataset');
    plotandprintOutliers(data_OS_pre_out, isOutlierMask_OS, 'Outliers from OS Dataset');

    % Step 4: Calculate number of samples and time increment (dt) for IS data
    N = height(data_IS);  % Number of observations in IS dataset
    t = yearfrac(data_IS.Time(1), data_IS.Time(end), 1);  % Total time in years (using actual calendar)
    dt = t / N;  % Time step between observations

    % Step 5: Bootstrap Ornstein-Uhlenbeck (OU) parameters using IS returns
    P = ouBootstrap(data_IS.Rt, dt, M, alpha);

    % Step 6: Filter out bootstrap samples with negative mean reversion rate 'k'
    k_positive_mask = P.bootstrapParams.k > 0; % Logical index where k is positive
    num_removed = sum(~k_positive_mask);  % Count how many samples were removed
    if num_removed > 0

        P.bootstrapParams = P.bootstrapParams(k_positive_mask, :);  % Keep only positive k samples
        fprintf("With a percentage of negative k removed from the bootstrap of %.4f%%\n", ...
            num_removed / M * 100);
        fprintf("Total amount of sample removed is %d\n", num_removed);
        
    end

    % Step 7: Print estimated OU parameters and removal summary
    ou_parameters_print(P);

    % Step 8: Compute transaction costs based on input data (implementation dependent)
    C = transaction_cost(data_IS);

    % Step 9: Run trading strategy optimization and evaluation on OS data
    [optimal_bands, mu_OS] = trading_strategy(P, f_list, l_list, M, alpha, C, data_OS);

    % Step 10: Plot returns versus leverage for stop-loss levels in l_list
    mu_leverage_plot(l_list, f_up, P, optimal_bands, C, alpha); 

    % Collect results into a structured output
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
end
