function [clean_IS, outliers, isOutlier] = removeOutliersIS(data_IS)
%REMOVEOUTLIERSIS  Removes log-spread and antipersistent outliers from IS data.
%
%   [clean_IS, outliers, isOutlier] = REMOVEOUTLIERSIS(data_IS)
%
%   This function filters in-sample (IS) price data by applying two 
%   sequential outlier detection procedures:
%
%     1. Log-spread outliers: Observations where the 'Rt' (log-spread)
%        deviates significantly from the interquartile range (IQR).
%
%     2. Antipersistent outliers: Observations that violate statistical 
%        properties of mean-reversion or consistency over time.
%
%   INPUT:
%     data_IS    - A table containing IS data with a 'Rt' column (log-spread).
%
%   OUTPUT:
%     clean_IS   - The filtered table with outliers removed.
%     outliers   - A table containing all detected outliers.
%     isOutlier  - A logical array flagging the outlier rows in the original data.

    % Step 1: Remove log spread IQR outliers
    [data_no_spread_outliers, isOutlier_log, ~] = filterLogSpreadOutliers(data_IS);

    % Step 2: Remove antipersistent outliers
    [clean_IS, isOutlier_anti, ~] = filterAntipersistentOutliers(data_no_spread_outliers);
    
    % Reconstruct full outlier mask
    isOutlier = [isOutlier_log; isOutlier_anti];
    
    % Create outlier table
    outliers = data_IS(isOutlier, :);
end
