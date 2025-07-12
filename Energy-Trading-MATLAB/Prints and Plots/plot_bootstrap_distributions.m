function plot_bootstrap_distributions(S, M)
%PLOT_BOOTSTRAP_DISTRIBUTIONS  Plot bootstrap distributions of OU parameters.
%
%   plot_bootstrap_distributions(S, M)
%
%   INPUTS:
%     S - struct returned by ouBootstrap, containing:
%          S.bootstrapParams: table with M rows and columns {'k', 'eta', 'sigma'}
%          S.CI: table with rows {'Lower', 'Upper'} and same columns
%     M - number of bootstrap samples
%
%   OUTPUT:
%     A 3-panel histogram figure with CI and mean lines.

    if nargin < 2
        M = height(S.bootstrapParams);
    end

    % Ensure export directory exists
    if ~exist('Plots', 'dir')
        mkdir('Plots');
    end

    figure('Position', [100 100 700 900], 'Color', 'w');

    % Labels
    paramNames = {'k', 'eta', 'sigma'};
    latexNames = {'\kappa', '\eta', '\sigma'};
    colors     = {[0.2 0.4 0.6], [0.2 0.6 0.4], [0.6 0.4 0.2]};

    for i = 1:3
        param = paramNames{i};
        values = S.bootstrapParams.(param);
        CI = [S.CI{1, param}, S.CI{2, param}];
        mu = mean(values);

        subplot(3,1,i);
        histogram(values, 30, ...
            'FaceColor', colors{i}, ...
            'EdgeColor', 'none', ...
            'Normalization', 'pdf');
        hold on;

        % CI lines
        xline(CI(1), '--r', 'LineWidth', 1.5, ...
              'Label', 'CI Lower', ...
              'LabelOrientation', 'horizontal', ...
              'LabelHorizontalAlignment', 'left');
        xline(CI(2), '--r', 'LineWidth', 1.5, ...
              'Label', 'CI Upper', ...
              'LabelOrientation', 'horizontal', ...
              'LabelHorizontalAlignment', 'right');

        % Mean line
        xline(mu, '-', 'Color', 'k', 'LineWidth', 1.5, ...
              'Label', 'Mean', ...
              'LabelOrientation', 'horizontal', ...
              'LabelHorizontalAlignment', 'center');

        % Axes and labels
        title(sprintf('Bootstrap Distribution of $%s$ (%d samples)', latexNames{i}, M), ...
              'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
        xlabel(['$' latexNames{i} '$'], 'Interpreter', 'latex', 'FontSize', 11);
        ylabel('Density', 'FontSize', 11);
        grid on; box on;
        hold off;
    end

    sgtitle('Parametric Bootstrap: OU Parameter Estimates and 95\% CI', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

    % Optional: save the plot
    % exportgraphics(gcf, 'Plots/MLE_Parameters.pdf', 'ContentType', 'vector', ...
    %     'BackgroundColor', 'white', 'Resolution', 300);
end