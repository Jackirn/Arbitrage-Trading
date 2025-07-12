function plot_bands(d_opt, u_opt, c_values)
%PLOT_BANDS Plots optimal entry and exit bands vs transaction cost.
%   plot_bands(d_opt, u_opt, c_values) plots:
%     - |d| (entry band) in red
%     - u (exit band) in blue
%   with a vertical line at c^* = 0.76.
%
%   Inputs:
%     d_opt     - Vector of optimal d values (entry band)
%     u_opt     - Vector of optimal u values (exit band)
%     c_values  - Corresponding transaction cost values

    figure;
    plot(c_values, d_opt, 'r-', 'LineWidth', 2); hold on;
    plot(c_values, u_opt, 'b--', 'LineWidth', 2);
    xline(0.76, 'k:', 'c^*', 'LabelVerticalAlignment', 'bottom');
    xlabel('Transaction cost c (in S units)');
    ylabel('Optimal trading bands');
    title('Optimal d and u vs Transaction Cost (l = -1.960)');
    legend('|d| (Entry)', 'u (Exit)', 'Location', 'northwest');
    grid on;
    %exportgraphics(gcf, 'Plots/Optimal Bands (figure 3).pdf', 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', 300);
end