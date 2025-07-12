function [optimal_bands, mu_OS] = trading_strategy(P, f_list, l_list, M, alpha, C, OS_data)
% TRADING_STRATEGY Evaluate optimal trading bands and compute OOS returns
%
% INPUTS
%   P         – OU process parameters (struct with .parameters and .bootstrapParams)
%   f_list    – Cell array of leverage values (numeric or 'opt')
%   l_list    – Stop-loss levels (in sigma-units, negative)
%   M         – Number of bootstrap replications
%   alpha     – Confidence level for CI
%   C         – Transaction cost (log-scale)
%   OS_data   – Struct or matrix containing out-of-sample price data
%
% OUTPUTS
%   optimal_bands – Struct array with estimated trading bands and CIs
%   mu_OS         – Out-of-sample (OOS) returns for each leverage strategy

    % Compute stationary parameters from OU process
    
    cap_sigma = P.parameters.sigma / sqrt(2 * P.parameters.k); % Long-run standard deviation
    cap_eta   = P.parameters.eta;                              % Long-run mean
    
    % Initialize output containers
    
    N = numel(f_list); % Number of leverage levels
    L = numel(l_list); % Number of stop-loss levels
    
    mu_OS = zeros(N, L); % Matrix to store OS returns
    
    % Preallocate empty struct array
    template = struct( ...
        'd_estimated',  NaN, ...
        'u_estimated',  NaN, ...
        'mu_estimated', NaN, ...
        'd_CI',         [NaN NaN], ...
        'u_CI',         [NaN NaN], ...
        'mu_CI',        [NaN NaN], ...
        'f_estimated',  NaN, ...
        'f_opt_CI',     [NaN NaN], ...
        'f_input',      [], ...
        'profitable_sample', NaN ...
    );
    
    optimal_bands(N, L) = template;
    grid_data = cell(N, L);
    
    % Main loop over stop-loss levels (each l gets its own figure)
    
    for j = 1:L
    
        iPlot = 0;
        l = l_list{j};
    
        % Create figure layout for OS returns
        fig = figure;
        tlo = tiledlayout(N, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(sprintf('Out-of-Sample trading bands (l = %.3f sigma)', l), 'FontWeight', 'bold');
        
    
        fprintf('\n=============================================\n');
        fprintf(' Evaluating Optimal Trading Bands with Stop-Loss\n');
        fprintf('=============================================\n');
        fprintf('Using stop-loss level: l = %.3f (sigma units)\n\n', l);
    
        for i = 1:N
    
            f = f_list{i};
    
            % Compute optimal bands
            [R, grid_data{i}] = optimal_trading_bands(M, l, f, P, C, alpha, true);
            optimal_bands(i,j) = R;
    
            % Determine leverage used
            if isequal(f, 'opt')
                f_used = optimal_bands(i,j).f_estimated;
            else
                f_used = f;
            end
    
            % Compute OS return
            mu_OS(i,j) = OS_return(optimal_bands(i,j).u_estimated, ...
                                   optimal_bands(i,j).d_estimated, ...
                                   l, f_used, OS_data, C, cap_sigma, cap_eta, tlo);
    
            % Print results
            fprintf('------------------------------------------------\n');
            fprintf(' Results for f = %s (used f = %.4f)\n', num2str(f), f_used);
            fprintf('------------------------------------------------\n');
    
            if ~isnan(optimal_bands(i,j).u_estimated) && ~isnan(optimal_bands(i,j).d_estimated)
                fprintf('Estimated u*   : %.3f\n', optimal_bands(i,j).u_estimated);
                fprintf('Estimated d*   : %.3f\n', -optimal_bands(i,j).d_estimated);
            else
                fprintf(' No estimated bands since f optimal estimated is equal to 0\n\n')
            end
    
            fprintf('Estimated mu*  : %.3f\n\n', optimal_bands(i,j).mu_estimated);
    
            if ~isnan(optimal_bands(i,j).d_CI)
    
                print_CI_intervals(optimal_bands(i,j).u_CI, ...
                               optimal_bands(i,j).d_CI, ...
                               optimal_bands(i,j).mu_CI);
            else
                fprintf('No confidence interval founded since there are no profitable samples\n\n')
            end
            if isequal(f, 'opt') && optimal_bands(i,j).f_opt_CI(1)~=0.0 && optimal_bands(i,j).f_opt_CI(2)~=0.0
                fprintf('------------------------------\n');
                fprintf('  95%% CI for optimal f\n');
                fprintf('  f* : [%.4f , %.4f]\n', ...
                         optimal_bands(i,j).f_opt_CI(1), ...
                         optimal_bands(i,j).f_opt_CI(2));
                fprintf('------------------------------\n');
            end
    
            if optimal_bands(i,j).profitable_sample < 1
                fprintf('Profitable sample for CI computation: %.3f%%\n\n', ...
                    optimal_bands(i,j).profitable_sample * 100);
            end
    
            if ~isnan(mu_OS(i,j))
                fprintf('Out-of-sample return: %.3f\n\n', mu_OS(i,j));
            end
        end
    
        % Find all lines and constant lines with DisplayName inside your figure
        linesWithNames = findall(fig, 'Type', 'Line', '-or', 'Type', 'ConstantLine');
        linesWithNames = flipud(linesWithNames); % Optional: flip to control order
        
        % Get unique handles by DisplayName to avoid duplicates
        [~, ia] = unique(arrayfun(@(h) string(get(h, 'DisplayName')), linesWithNames), 'stable');
        linesUnique = linesWithNames(ia);
        
        % Create legend on one of the axes or the figure
        lgd = legend(linesUnique, 'Location', 'southoutside', 'Orientation', 'horizontal');
        
        % Attach the legend to the bottom tile of your tiledlayout
        lgd.Layout.Tile = 'south';
    
    
    % Figure: c-grid trading bands – one subplot for each f


    % - Determine which f-entries have at least one valid d value 
    valid_idx = false(1,N);
    for ii = 1:N
        valid_idx(ii) = ~all( isnan( grid_data{ii}.d ) );
    end
    plot_idx   = find(valid_idx);    % indices that will be plotted
    nPlots     = numel(plot_idx);
    
    if nPlots == 0
        warning('No plottable d-grid curves – figure skipped.');
    else
        % - Create layout sized to valid plots
        nCols  = ceil( sqrt(nPlots) );
        nRows  = ceil( nPlots / nCols );
    
        figGrid = figure;
        tloGrid = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');
        sgtitle(tloGrid, sprintf('Trading Bands  (l = %.3f sigma)', l), 'FontWeight','bold');
    
        sg     = grid_data{plot_idx(1)}.sigma; % same sigma grid for all
        c_grid = C ./ sg;
        xMin   = min(c_grid) - 0.01;
        xMax   = max(c_grid) + 0.01;
    
        % ─ Loop over plottable indices only
        for k = 1:nPlots
            i  = plot_idx(k); % original leverage index
    
            ug  = grid_data{i}.u;
            dg  = grid_data{i}.d;
            dCI = optimal_bands(i,j).d_CI;
            uCI = optimal_bands(i,j).u_CI;
    
            nexttile(tloGrid);
            hold on;  box on;  grid on;
            xlim([xMin xMax]);
    
            % curves
            plot(c_grid, ug, 'b-', 'LineWidth',1.4, 'DisplayName','u(sigma)');
            plot(c_grid, dg, 'r-', 'LineWidth',1.4, 'DisplayName','|d(sigma)|');
    
            % CI lines
            yline(dCI(1),'--r','DisplayName','CI lower (d)');
            yline(dCI(2),'--r','DisplayName','CI upper (d)');
            yline(uCI(1),'--b','DisplayName','CI lower (u)');
            yline(uCI(2),'--b','DisplayName','CI upper (u)');
    
            % point & vertical line at chosen cost
            xline(C/cap_sigma,'--k','DisplayName','Estimated bands');
            plot(C/cap_sigma, optimal_bands(i,j).d_estimated,'rx','MarkerSize',10,'LineWidth',1.5, 'DisplayName','|d| estimated');
            plot(C/cap_sigma, optimal_bands(i,j).u_estimated,'bx','MarkerSize',10,'LineWidth',1.5, 'DisplayName', 'u estimated');
    
            % shaded CI areas (not shown in legend)
            xPatch = [xMin xMax xMax xMin];
            fill(xPatch,[dCI(1) dCI(1) dCI(2) dCI(2)], [1 .8 .8], 'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
            fill(xPatch,[uCI(1) uCI(1) uCI(2) uCI(2)], [.8 .9 1], 'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off');
    
            title(sprintf('f = %s', num2str(f_list{i})));
            xlabel('Transaction cost c  (sigma-units)');
            ylabel('Band level  (sigma-units)');
            
            if ~isnumeric(f_list{i})
                iPlot = i;
            end
        end
    
        % ─ Shared legend (lines only, unique names)
        linesWithNames = findall(figGrid,'Type','Line','-or','Type','ConstantLine');
        linesWithNames = flipud(linesWithNames);
        [~, ia] = unique(arrayfun(@(h) string(get(h,'DisplayName')), linesWithNames),'stable');
        lgd = legend(linesWithNames(ia), 'Location','southoutside', 'Orientation','horizontal');
        lgd.Layout.Tile = 'south';  % Correct way to attach legend to bottom tile in a tiledlayout
    
    end
    
        % Optional: plot step function for optimal f (if estimated)

        if iPlot > 0
    
            f_vec = grid_data{iPlot}.f(:).';
            f_zero = f_vec;      
            f_zero(f_vec ~= 0) = NaN;
            f_nonzero = f_vec;   
            f_nonzero(f_vec == 0) = NaN;
    
            figF = figure;
            hold on; box on; grid on;
    
            xlim([xMin xMax]);
    
            plot(c_grid, f_zero,    'k-', 'LineWidth',1.5, 'DisplayName','f = 0');
            plot(c_grid, f_nonzero, 'b-', 'LineWidth',1.5, 'DisplayName','f ≠ 0');
    
            y1 = optimal_bands(iPlot,j).f_opt_CI(1);
            y2 = optimal_bands(iPlot,j).f_opt_CI(2);
    
            x_patch = [xMin, xMax, xMax, xMin];
            y_patch = [y1 y1 y2 y2];
    
            fill(x_patch, y_patch, [0.8 0.9 1],  'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
            yline(y1, '--b', 'DisplayName','CI lower (f)');
            yline(y2, '--b', 'DisplayName','CI upper (f)');
            xline(C/cap_sigma, '--k', 'DisplayName',' Estimated transaction cost');
    
            plot(C/cap_sigma, optimal_bands(i,j).f_estimated, 'bx', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'f estimated');
    
            xlabel('Transaction cost c  (\sigma-units)');
            ylabel('f_{opt}(c)');
            title(sprintf('Stepwise plot of optimal f using stop-loss level: l = %.3f (sigma units)\n\n', l));
            legend('Location','best');
        end
    end
    

    % Convert results into tables for readability


    rowNames = cellfun(@num2str, f_list, 'UniformOutput', false);
    colNames = cellfun(@num2str, l_list, 'UniformOutput', false);
    
    optimal_bands = cell2table(mat2cell(optimal_bands, ones(1,N), ones(1,L)), ...
                                'RowNames', rowNames, 'VariableNames', colNames);
    
    mu_OS = array2table(mu_OS, 'RowNames', rowNames, 'VariableNames', colNames);

end
