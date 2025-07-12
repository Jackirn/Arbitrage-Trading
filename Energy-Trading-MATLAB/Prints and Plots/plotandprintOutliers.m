function plotandprintOutliers(data, outlier_indices, titleStr)
% PLOTANDPRINTOUTLIERS Plots Rt over time, highlights and prints outliers.
%
% data           : table with 'Time' (datetime) and 'Rt' columns
% outlier_indices: logical vector or indices of outliers
% titleStr       : title for the plot
%
% Prints summary of outliers to Command Window.

    % Styling parameters
    fsTitle = 16;     % Title font size
    fsLabel = 14;     % Axis label font size
    fsTicks = 12;     % Tick/Legend font size
    lw      = 1.8;    % Line width

    % Create figure
    figure('Color', 'w', 'Position', [100 100 850 400]);
    hold on;
    box on;
    grid on;
    set(gca, 'FontSize', fsTicks, 'LineWidth', 1);

    % Plot Rt
    plot(data.Time, data.Rt, '-', ...
         'Color', [0.1 0.4 0.8], ...
         'LineWidth', lw, ...
         'DisplayName', '$R_t$');

    % Create outlier mask
    if islogical(outlier_indices)
        mask = outlier_indices;
    else
        mask = false(height(data),1);
        mask(outlier_indices) = true;
    end

    % Plot outliers
    plot(data.Time(mask), data.Rt(mask), 'd', ...
         'MarkerSize', 10, ...
         'MarkerEdgeColor', 'r', ...
         'MarkerFaceColor', [1 0.6 0.6], ...
         'LineWidth', 1.5, ...
         'DisplayName', 'Outliers');

    % Set y-axis limits with margin
    y_margin = 0.05 * range(data.Rt);
    ylim([min(data.Rt)-y_margin, max(data.Rt)+y_margin]);

    % Set x-axis limits
    if isdatetime(data.Time)
        xlim([min(data.Time), max(data.Time)]);
        xtickformat('yyyy-MM-dd');
    end

    % Axis labels and title
    xlabel('Time', 'FontSize', fsLabel, 'Interpreter', 'latex', 'FontWeight', 'bold');
    ylabel('$R_t$ (log-spread)', 'FontSize', fsLabel, 'Interpreter', 'latex', 'FontWeight', 'bold');
    title(titleStr, 'FontSize', fsTitle, 'Interpreter', 'latex', 'FontWeight', 'bold');

    % Legend
    legend('Location', 'best', ...
           'FontSize', fsTicks, ...
           'Interpreter', 'latex', ...
           'Box', 'off');

    % Print outlier summary
    n_total = height(data);
    n_outliers = nnz(mask);
    perc = 100 * n_outliers / n_total;
    outlierTab = data(mask, :);

    fprintf('---------------------------------------------\n');
    fprintf('           %s: %d out of %d (%.2f%%)\n', titleStr, n_outliers, n_total, perc);
    fprintf('---------------------------------------------\n');
    if n_outliers > 0
        disp(outlierTab);
    else
        fprintf('No outliers detected.\n');
    end
    fprintf('---------------------------------------------\n\n');
end
