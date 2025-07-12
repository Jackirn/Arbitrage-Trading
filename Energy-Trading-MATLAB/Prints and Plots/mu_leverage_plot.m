function mu_leverage_plot(l_list, f_up, P, optimal_bands, C, alpha)
%MU_LEVERAGE_PLOT  Plot long-run return mu(f) as a function of leverage f.
%
%   mu_leverage_plot(l_list, f_up, P, optimal_bands, C, alpha)
%
%   For each stop-loss level `l` in `l_list`, this function:
%     - Plots the curve mu(f), the estimated long-run return as a function of leverage f.
%     - Marks the estimated optimal leverage f* obtained from `optimal_bands`.
%     - Organizes the plots in a compact tiled layout.
%     - Displays a shared legend for all subplots.
%
%   INPUTS:
%     l_list        - Cell array of stop-loss levels (in sigma-units, negative)
%     f_up          - Maximum leverage to consider in the plot (upper bound for f)
%     P             - OU process parameters (struct from ouBootstrap)
%     optimal_bands - Table of optimal band results (from trading_strategy)
%     C             - Transaction cost (log-scale)
%     alpha         - Significance level for CI (used internally by optimal_trading_bands)
%
%   NOTE:
%     - The function assumes 'opt' is the row name in optimal_bands for f = f_opt.
%     - Uses optimal_trading_bands with M = 1 for fast evaluation of mu(f).

    f = 'opt';         % Row name for optimal leverage results
    N = numel(l_list); % Number of stop-loss levels
    
    % Dynamic tiled layout size 
    if N <= 3                % Vertical layout
        rows = N;  cols = 1;
    elseif N == 4            % Perfect square layout
        rows = 2;  cols = 2;
    else                     % Nearly square grid
        rows = ceil(sqrt(N));
        cols = ceil(N / rows);
    end

    fig = figure;
    tlo = tiledlayout(rows, cols, 'TileSpacing','compact','Padding','compact');
    
    allLines = gobjects(0); % Collect handles for global legend
    
    for i = 1:N
        ax = nexttile(tlo); % Draw in the next tile
        l  = l_list{i};     % Current stop-loss level
    
        % Compute mu(f) curve 
        f_vals  = linspace(0, f_up, 100); % Leverage values
        mu_vals = arrayfun(@(f_i) ...     % Compute mu(f_i) for each f
                    optimal_trading_bands(1, l, f_i, P, C, alpha).mu_estimated, ...
                    f_vals);
    
        % Plot mu(f) curve
        h1 = plot(ax, f_vals, mu_vals, 'LineWidth',1.5,'DisplayName','\mu(f)');
        hold(ax,'on')
    
        % Plot optimal f point 
        f_opt = optimal_bands{f, num2str(l)}.f_estimated; % Optimal leverage
        y_opt = interp1(f_vals, mu_vals, f_opt);          % Interpolated mu(f*)
        h2 = plot(ax, f_opt, y_opt, 'ro', 'MarkerSize',8, ...
                  'LineWidth',2, 'DisplayName','Optimal f');
    
        % Axis styling and labels
        xlabel(ax,'Leverage f'); 
        ylabel(ax,'Long-run return \mu');
        title(ax,sprintf('l = %.3f', l)); 
        grid(ax,'on'); 
        hold(ax,'off');
    
        % Store handles for global legend
        allLines = [allLines; h1; h2];
    end
    
    % Shared legend across all subplots 

    % Filter valid handles with unique DisplayNames
    valid = allLines(ishandle(allLines) & ...
             arrayfun(@(h) ~isempty(h.DisplayName), allLines));
    names = arrayfun(@(h) string(h.DisplayName), valid, 'UniformOutput', false);
    [~, ia] = unique([names{:}], 'stable'); % Keep first occurrence
    
    uniqueLines = valid(ia);                % Unique line handles
    
    % Create shared legend in the bottom tile
    lgd = legend(uniqueLines, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';

end
