function [clean_data, outlier_indices, outlier_values] = filterAntipersistentOutliers(data)
%FILTERANTIPERSISTENTOUTLIERS  Removes antipersistent outliers from Rt time series.
%
%   [clean_data, outlier_indices, outlier_values] = FILTERANTIPERSISTENTOUTLIERS(data)
%
%   Identifies and removes local outliers in the log-spread ('Rt') column that
%   exhibit antipersistent behavior — i.e., large jumps in both backward and 
%   forward directions compared to the IQR.
%
%   INPUT:
%     data             - A table containing at least the 'Rt' column (log-spread).
%
%   OUTPUT:
%     clean_data       - Table with antipersistent outliers removed.
%     outlier_indices  - Indices of the detected outlier rows.
%     outlier_values   - Table containing the outlier rows.
%
%   METHOD:
%     For each interior time point t (2 ≤ t ≤ n−1), the function checks:
%
%         |Rt(t) - Rt(t−1)| > IQR   AND   |Rt(t+1) - Rt(t)| > 0.95 × IQR
%
%     If both conditions are met, Rt(t) is flagged as an outlier.
%
%   This heuristic captures abrupt, transient deviations in Rt inconsistent 
%   with a mean-reverting (persistent) time series.

    Rt = data.Rt;
    n = height(data);
    
    % Compute IQR
    IQR_val = iqr(Rt);
    
    % Initialize logical index for outliers
    isOutlier = false(n, 1);
    
    % Apply the 3-point rule
    for t = 2:(n-1)
        delta_prev = abs(Rt(t) - Rt(t-1));
        delta_next = abs(Rt(t+1) - Rt(t));
        
        if delta_prev > IQR_val && delta_next > 0.95 * IQR_val
            isOutlier(t) = true;
        end
    end
    
    % Find outlier indices and values
    outlier_indices = find(isOutlier);
    outlier_values = data(outlier_indices, :);

    % Remove them from data
    clean_data = data(~isOutlier, :);
end
