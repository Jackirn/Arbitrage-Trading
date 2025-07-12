function [cleanData, outlier_indices, outlier_values] = filterLogSpreadOutliers(data)
%FILTERLOGSPREADOUTLIERS  Remove log-spread (Rt) outliers from a price table.
%
%   [cleanData, outlier_indices, outlier_values] = FILTERLOGSPREADOUTLIERS(data)
%
%   Identifies and removes outliers in the log-spread ('Rt') column of the
%   input table using a 3×IQR rule.
%
%   INPUT:
%     data             - A table containing at least the 'Rt' column (log-spread).
%
%   OUTPUT:
%     cleanData        - Table with outlier rows removed.
%     outlier_indices  - Indices of the detected outlier rows in the input table.
%     outlier_values   - Table containing the outlier rows.
%
%   METHOD:
%     An observation Rt(i) is considered an outlier if it lies outside the
%     interval:
%        [Q1 - 3×IQR, Q3 + 3×IQR]
%     where Q1 and Q3 are the 25th and 75th percentiles of Rt, and
%     IQR = Q3 - Q1.

    Rt = data.Rt;
    
    Q1       = quantile(Rt, 0.25);
    Q3       = quantile(Rt, 0.75);
    iqrVal   = Q3 - Q1;
    loBound  = Q1 - 3 * iqrVal;
    hiBound  = Q3 + 3 * iqrVal;
    isOutlier = Rt < loBound | Rt > hiBound;
        
    outlier_indices = find(isOutlier);
    outlier_values = data(outlier_indices, :);

    cleanData   = data(~isOutlier, :);
end
