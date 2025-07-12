import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def ou_mle(x, dt):
    """
    Closed-form MLE for the Ornstein-Uhlenbeck process

    Parameters
    ----------
    x : array-like
        Time series observations (1D array: [X_0, ..., X_N])
    dt : float
        Time step Δt

    Returns
    -------
    k : float
        Mean reversion speed
    eta : float
        Long-term mean
    sigma : float
        Volatility
    """

    x = np.asarray(x).flatten()  # ensure 1D array
    x_minus = x[:-1]
    x_plus = x[1:]
    N = len(x_minus)

    Y_m = np.mean(x_minus)
    Y_p = np.mean(x_plus)
    Y_mm = np.mean(x_minus ** 2)
    Y_pp = np.mean(x_plus ** 2)
    Y_pm = np.mean(x_minus * x_plus)

    rho = (Y_pm - Y_m * Y_p) / (Y_mm - Y_m ** 2)
    k = -np.log(rho) / dt

    eta = Y_p + ((x[-1] - x[0]) / N) * \
          (Y_pm - Y_m * Y_p) / ((Y_mm - Y_m ** 2) - (Y_pm - Y_m * Y_p))

    sigma2 = Y_pp - Y_p ** 2 - ((Y_pm - Y_m * Y_p) ** 2) / (Y_mm - Y_m ** 2)
    sigma = np.sqrt((2 * k * sigma2) / (1 - np.exp(-2 * k * dt)))

    return k, eta, sigma



def ou_sim(x0, k, eta, sigma, dt, N, rng=None):
    """
    Simulate a single path of an Ornstein-Uhlenbeck (OU) process.

    Exact one-step solution on an equally spaced grid Δt

    Parameters
    ----------
    x0 : float
        Initial state X₀.
    k : float
        Mean-reversion speed (k > 0).
    eta : float
        Long-run mean eta.
    sigma : float
        Volatility parameter sigma > 0.
    dt : float
        Time-step length Δt > 0.
    N : int
        Number of steps (path length will be N+1).
    rng : np.random.Generator, optional
        NumPy random generator for reproducibility.

    Returns
    -------
    x : np.ndarray
        1-D array of length N+1 with the simulated path (x[0] = x0).
    """
    if rng is None:
        rng = np.random.default_rng()

    # Pre-allocate output array
    x = np.empty(N + 1)
    x[0] = x0

    # Pre-compute constants
    a  = np.exp(-k * dt)                     
    b  = eta * (1.0 - a)                    
    sd = sigma * np.sqrt((1.0 - a**2) / (2.0 * k))  # conditional std-dev

    # Draw N standard normal variates
    z = rng.standard_normal(N)

    # Exact simulation loop
    for i in range(N):
        x[i + 1] = a * x[i] + b + sd * z[i]

    return x


def ou_bootstrap(clean_data, M=1_000, alpha=0.05, rng=None):
    """
    Parametric bootstrap for an Ornstein-Uhlenbeck process.

    Parameters
    ----------
    x : array-like
        Observed OU path  (length N+1).
    dt : float
        Sampling step Δt.
    M : int, default 1000
        Number of bootstrap replications.
    alpha : float, default 0.05
        Two-sided confidence level (e.g. 0.05 → 95 % CI).
    rng : np.random.Generator, optional
        NumPy random generator for reproducibility.

    Returns
    -------
    result : dict
        {
          'parameters'       : DataFrame with MLE (k, eta, sigma),
          'bootstrap_params' : DataFrame with M bootstrap estimates,
          'CI'               : DataFrame with lower/upper percentile CI,
          'median'           : DataFrame with bootstrap medians
        }
    """
    x = clean_data['Rt']

    t = pd.to_datetime(clean_data['Time'])
    t_int = t.iloc[-1] - t.iloc[0]
    n = len(clean_data)
    dt = t_int.total_seconds() / (n * 60 * 60 * 24 * 365) # yearly constant dt

    if rng is None:
        rng = np.random.default_rng()

    x = np.asarray(x).flatten()
    N = len(x) - 1 # number of simulation steps

    k_hat, eta_hat, sigma_hat = ou_mle(x, dt)

    boot_estimates = np.zeros((M, 3))
    for m in range(M):
        x_sim = ou_sim(x0=x[0],
                       k=k_hat,
                       eta=eta_hat,
                       sigma=sigma_hat,
                       dt=dt,
                       N=N,
                       rng=rng)
        boot_estimates[m] = ou_mle(x_sim, dt)

    col_names = ['k', 'eta', 'sigma']

    parameters       = pd.DataFrame([ [k_hat, eta_hat, sigma_hat] ],
                                    columns=col_names)

    bootstrap_params = pd.DataFrame(boot_estimates, columns=col_names)

    lower_p, upper_p = alpha * 50, 100 - alpha * 50
    CI_values = np.percentile(boot_estimates, [lower_p, upper_p], axis=0)

    CI = pd.DataFrame(CI_values,
                      index=['Lower', 'Upper'],
                      columns=col_names)

    median_vals = pd.DataFrame(bootstrap_params.median().to_frame().T,
                               columns=col_names)

    return {
        'parameters'      : parameters,
        'bootstrap_params': bootstrap_params,
        'CI'              : CI,
        'median'          : median_vals
    }


def print_ou_estimates(results):
    """
    Pretty print OU parameter estimates and confidence intervals from a dict.

    Parameters:
    -----------
    results : dict
        Dictionary with keys 'parameters' and 'CI'.
        'parameters': dict or pd.Series with ['k', 'eta', 'sigma']
        'CI': pd.DataFrame with rows ['Lower', 'Upper'] and columns ['k', 'eta', 'sigma'].
    """
    params = results['parameters']
    ci = results['CI']

    print("Ornstein-Uhlenbeck Parameter Estimates\n" + "-" * 45)

    for param in ['k', 'eta', 'sigma']:
        # Safely extract values
        est = params[param]
        lower = ci.loc['Lower', param]
        upper = ci.loc['Upper', param]

        # Convert from Series/arrays to scalars if needed
        if hasattr(est, 'item'):
            est = est.item()
        if hasattr(lower, 'item'):
            lower = lower.item()
        if hasattr(upper, 'item'):
            upper = upper.item()

        print(f"{param:6s}: Estimate = {est:.6f}, 95% CI = [{lower:.6f}, {upper:.6f}]")


def plot_bootstrap_distributions(results, M, filename='Plots/MLE_parameters.png'):
    """
    Plot bootstrap distributions of OU parameters with CI and mean lines.

    Parameters
    ----------
    results : dict
        Output from ou_bootstrap containing bootstrap estimates and confidence intervals.
    M : int
        Number of bootstrap samples.
    filename : str, optional
        Path to save the plot.
    """
    import os

    # Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    params = ['k', 'eta', 'sigma']
    colors = [(0.2, 0.4, 0.6), (0.2, 0.6, 0.4), (0.6, 0.4, 0.2)]
    bootstrap_params = results['bootstrap_params']
    CI = results['CI']

    fig, axes = plt.subplots(3, 1, figsize=(8, 10))
    fig.subplots_adjust(hspace=0.4)
    fig.patch.set_facecolor('white')

    for i, param in enumerate(params):
        ax = axes[i]
        data = bootstrap_params[param].values
        lower, upper = CI.loc['Lower', param], CI.loc['Upper', param]
        mu = np.mean(data)

        ax.hist(data, bins=30, density=True, color=colors[i], edgecolor='none')
        ax.axvline(lower, color='red', linestyle='--', linewidth=1.5, label='CI Lower')
        ax.axvline(upper, color='red', linestyle='--', linewidth=1.5, label='CI Upper')
        ax.axvline(mu, color='black', linestyle='-', linewidth=1.5, label='Mean')

        ax.set_title(f"Bootstrap Distribution of {param} ({M} samples)", fontsize=12, fontweight='bold')
        ax.set_xlabel(param, fontsize=11)
        ax.set_ylabel("Density", fontsize=11)
        ax.grid(True)
        ax.legend()

    fig.suptitle("Parametric Bootstrap: OU Parameter Estimates and 95% CI", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename, dpi=300)
    plt.show()

def plot_optimal_trading_bands(c_values, d_opt, u_opt, l,
                               k_val=1.0, sigma_val=1.0,
                               c_sigma_max=0.76,
                               figsize=(8, 5)):
    """
    Plot optimal entry (|d*|) and exit (u*) bands as a function of transaction cost.

    Parameters
    ----------
    c_values : array-like
        Transaction costs in “natural” units (i.e. not sigma-units).
    d_opt : array-like
        Optimal lower band estimates (negative values; will be plotted as positive |d*|).
    u_opt : array-like
        Optimal upper band estimates.
    l : float
        Stop-loss level in sigma-units (negative).
    k_val : float, optional
        OU speed parameter, default 1.0.
    sigma_val : float, optional
        OU long-run volatility, default 1.0.
    c_sigma_max : float, optional
        Maximum transaction cost (in sigma-units) for the x-axis limit.
    figsize : tuple, optional
        Figure size passed to plt.figure().

    Returns
    -------
    None
    """
    # Convert to sigma-units for the x-axis
    c_plot = np.asarray(c_values) * np.sqrt(2 * k_val) / sigma_val

    # Mask out any invalid results
    mask = (~np.isnan(d_opt)) & (~np.isnan(u_opt))

    # Prepare the figure
    plt.figure(figsize=figsize)
    plt.plot(c_plot[mask], np.abs(d_opt[mask]), '--', lw=2,
             label='|d*| (entry)')
    plt.plot(c_plot[mask], u_opt[mask],        '-' , lw=2,
             label='u*   (exit)')

    plt.xlabel('Transaction Cost c (in sigma units)')
    plt.ylabel('Optimal trading bands')
    plt.title(f'Optimal d and u vs Transaction Cost (l = {l:.3f})')
    plt.grid(True)
    plt.legend(loc='upper left')

    plt.xlim(0, c_sigma_max)
    ymax = np.nanmax(u_opt) * 1.05
    plt.ylim(0, ymax)

    plt.tight_layout()
    plt.show()

def plot_CI_mu(df_plot, alpha=0.5):
    """
    Plot the long-run return μ with confidence intervals against stop-loss levels.
    """
    import matplotlib.pyplot as plt

    x = df_plot['Stop-loss']
    y = df_plot['μ']
    lower = df_plot['mu_CI_low']
    upper = df_plot['mu_CI_high']

    plt.figure(figsize=(10, 6))
    
    # Confidence interval band
    plt.fill_between(x, lower, upper, color='lightsteelblue', alpha=alpha, label='95% CI')

    # μ line
    plt.plot(x, y, color='navy', marker='o', linewidth=2.5, label='μ')

    # Dashed border lines for CI (optional for clarity)
    plt.plot(x, lower, color='lightsteelblue', linestyle='--', linewidth=1)
    plt.plot(x, upper, color='lightsteelblue', linestyle='--', linewidth=1)

    plt.xlabel("Stop-loss level", fontsize=12)
    plt.ylabel("μ (long-run return)", fontsize=12)
    plt.title("μ vs Stop-loss Level with Confidence Band", fontsize=14)
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
