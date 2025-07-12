function mu_OS = OS_return(u, d, l, f, Data_OS, C, sigma, eta, tlo)
%OS_RETURN  Compute out-of-sample return from a band-based trading strategy.
%
%   mu_OS = OS_RETURN(u, d, l, f, Data_OS, C, sigma, eta, tlo)
%
%   Computes the annualized out-of-sample return from a long/short mean-reverting 
%   trading strategy based on thresholds derived from an OU process.
%
%   INPUTS:
%     u        - Upper exit band (in sigma-units)
%     d        - Lower entry band (in sigma-units)
%     l        - Stop-loss level (in sigma-units, negative)
%     f        - Leverage factor
%     Data_OS  - Table with out-of-sample data, must contain 'Rt' and 'Time'
%     C        - Log transaction cost
%     sigma    - Long-run standard deviation of the process
%     eta      - Long-run mean of the process
%     tlo      - (Optional) tiledlayout handle for plotting. If not provided, 
%                a new figure with two subplots is created
%
%   OUTPUT:
%     mu_OS    - Annualized out-of-sample return of the trading strategy
%
%   NOTES:
%     - The strategy combines both long and short symmetric band rules.
%     - Completed trades are counted using a band-crossing detection function.
%     - Visualizes entry/exit/stop-loss bands for long and short strategies.

    if nargin < 9, tlo = []; end     
   
    % Center the spread around long-run mean
    X_t = Data_OS.Rt - eta;
    
    % Define long strategy thresholds
    D = -d * sigma;
    U =  u * sigma;
    L =  l * sigma;
    
    % Count completed trades for long strategy
    [count_plus, count_minus] = count_d_to_u_l(X_t, D, U, L);
    
    % Define short strategy thresholds (symmetric)
    Dshort =  d * sigma;
    Ushort = -u * sigma;
    Lshort = -l * sigma;
    
    % Count completed trades for short strategy
    [count_minus_short, count_plus_short] = count_d_to_u_l(X_t, ...
                                              Dshort, Lshort, Ushort);
    
    % Compute annualized out-of-sample return (long + short strategies)
    mu_OS = ( (count_plus + count_plus_short)  * ...
              log(1 + f*(exp(U - D - C) - 1))   + ...
              (count_minus + count_minus_short) * ...
              log(1 + f*(exp(L - D - C) - 1)) ) ...
               / yearfrac(Data_OS.Time(1), Data_OS.Time(end), 1);
    
    % Plotting
    
    % Create axes: from layout if provided, otherwise new figure
    if isempty(tlo)
        fig = figure;
        ax_long  = subplot(1,2,1, 'Parent', fig);
        ax_short = subplot(1,2,2, 'Parent', fig);
    else
        ax_long  = nexttile(tlo); % LONG  (column 1)
        ax_short = nexttile(tlo); % SHORT (column 2)
    end
    
    % LONG PLOT 
    axes(ax_long); 
    plot(Data_OS.Time, X_t, 'k-', 'LineWidth', 1.2, 'DisplayName','OS Dataset'); hold on
    yline(D, 'b-', 'LineWidth', 1.2, 'DisplayName','Entry Band Level');
    yline(U, 'g-', 'LineWidth', 1.2,'DisplayName','Exit Band Level');
    yline(L, 'r-', 'LineWidth', 1.2,'DisplayName','Stop-Loss Level'); hold off
    title(sprintf('Long Bands with  f = %s', num2str(f)))
    xlabel('Time'); ylabel('X_t'); grid on
    
    % SHORT PLOT 
    axes(ax_short); 
    plot(Data_OS.Time, X_t, 'k-', 'LineWidth', 1.2,'DisplayName','OS Dataset'); hold on
    yline(Dshort, 'b-', 'LineWidth', 1.2, 'DisplayName','Entry Band Level');
    yline(Ushort, 'g-', 'LineWidth', 1.2,'DisplayName','Exit Band Level');
    yline(Lshort, 'r-', 'LineWidth', 1.2,'DisplayName','Stop-Loss Level'); hold off
    title(sprintf('Short Bands with f = %s', num2str(f)))
    xlabel('Time'); ylabel('X_t'); grid on
end
