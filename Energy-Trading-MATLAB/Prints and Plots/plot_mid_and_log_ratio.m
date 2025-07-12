function plot_mid_and_log_ratio(T_g)
%PLOT_MID_AND_LOG_RATIO  Plot mid prices and log-price ratio over time.
%
%   plot_mid_and_log_ratio(T_g)
%
%   Visualizes the time series of mid prices for two assets and their
%   log-price ratio, using a two-panel plot.
%
%   INPUT:
%     T_g - Table containing the following columns:
%           • Time  : datetime or numeric vector
%           • Mid1  : Mid price of asset 1
%           • Mid2  : Mid price of asset 2
%           • Rt    : Log-price ratio, typically log(Mid1 / Mid2)
%
%   OUTPUT:
%     A figure with:
%       - Top subplot: Mid1 and Mid2 time series
%       - Bottom subplot: Log-price ratio series (Rt)

    % Extract time range
    startTime = T_g.Time(1);
    endTime   = T_g.Time(end);

    % Create figure
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);

    % Fonts and colors
    fsTitle = 15;
    fsLabel = 13;
    fsAxis  = 11;
    lw = 1.8;

    % Subplot 1: Mid Prices
    subplot(2,1,1)
    hold on; grid on; box on;
    stairs(T_g.Time, T_g.Mid1, 'Color', [0.1 0.4 0.8], 'LineWidth', lw);
    stairs(T_g.Time, T_g.Mid2, 'Color', [0.85 0.33 0.1], 'LineWidth', lw);
    xlim([startTime, endTime]);

    xlabel('Time', 'FontSize', fsLabel, 'Interpreter', 'latex');
    ylabel('Mid Prices', 'FontSize', fsLabel, 'Interpreter', 'latex');
    title('Mid Prices of Asset 1 and Asset 2', ...
          'FontSize', fsTitle, 'FontWeight', 'bold', 'Interpreter', 'latex');
    legend({'Mid1', 'Mid2'}, 'FontSize', fsAxis, ...
           'Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
    set(gca, 'FontSize', fsAxis, 'LineWidth', 1);

    % Subplot 2: Log-Price Ratio
    subplot(2,1,2)
    hold on; grid on; box on;
    stairs(T_g.Time, T_g.Rt, 'Color', [0.2 0.6 0.3], 'LineWidth', lw);
    xlim([startTime, endTime]);

    xlabel('Time', 'FontSize', fsLabel, 'Interpreter', 'latex');
    ylabel('Log-Price Ratio', 'FontSize', fsLabel, 'Interpreter', 'latex');
    title('Log of Price Ratio ($R_t$)', ...
          'FontSize', fsTitle, 'FontWeight', 'bold', 'Interpreter', 'latex');
    legend({'$R_t$'}, 'FontSize', fsAxis, ...
           'Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
    set(gca, 'FontSize', fsAxis, 'LineWidth', 1);
end
