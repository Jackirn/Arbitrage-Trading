import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from math import pi
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.integrate import quad

def constraints_fun(x, c, l):
    """
    Inequality constraints for the optimization:
    c1: c + buffer - (u - d) <= 0
    c2: l - d <= 0
    c3: d - u <= 0
    """
    d, u = x
    buffer = 0.1
    return [
        c + buffer - (u - d),  # c1
        l - d,                 # c2
        d - u                  # c3
    ]

def optimal_trading_bands(M, l, f, P, C, alpha, grid):
    """
    Estimate optimal OU trading bands and confidence intervals.
    """
    R = {
        'd_estimated':  np.nan,
        'u_estimated':  np.nan,
        'mu_estimated': np.nan,
        'd_CI':         [np.nan, np.nan],
        'u_CI':         [np.nan, np.nan],
        'mu_CI':        [np.nan, np.nan],
        'f_estimated':  np.nan,
        'f_opt_CI':     [np.nan, np.nan],
        'f_input':      f if isinstance(f, (int, float)) else np.nan
    }

    # Extract OU parameters
    k_hat = P['parameters']['k'].values[0]
    sigma_hat = P['parameters']['sigma'].values[0]
    theta = 1 / k_hat
    sigma_stat = sigma_hat / np.sqrt(2 * k_hat)
    c = C / sigma_stat

    # Objective: maximize long-run return
    def obj_fun(x):
        mu, _ = long_return(x[0], x[1], c, l, sigma_stat, f)  
        return -mu

    d_min = l + 0.01      
    d_max = 0.6            
    u_min = l + C           
    u_max = 3.0              

    bounds = [(d_min, d_max), (u_min, u_max)]
    x0     = [-0.5, 0.5]     

    # Constraints 
    cons = {'type': 'ineq', 'fun': lambda x: -np.array(constraints_fun(x, c, l))}

    res = minimize(obj_fun, x0, method='SLSQP', bounds=bounds, constraints=cons, options={'ftol': 1e-8, 'disp': False, 'maxiter': 500})

    if res.success:
        d = abs(res.x[0])
        u = res.x[1]
        mu, f_star = long_return(-d, u, c, l, sigma_stat, f)
        R['d_estimated'] = d
        R['u_estimated'] = u
        R['mu_estimated'] = mu / theta
        if not isinstance(f, (int, float)):
            R['f_estimated'] = f_star
    else:
        return R  # Optimization failed

    # Bootstrap confidence intervals if M > 1
    if M > 1:
        sigma_vec = P['bootstrap_params']['sigma'].values / np.sqrt(2 * P['bootstrap_params']['k'].values)
        theta_vec = 1 / P['bootstrap_params']['k'].values
        sigma_grid = np.linspace(sigma_vec.min(), sigma_vec.max(), grid)

        d_grid, u_grid, mu_grid, f_grid = [], [], [], []

        for sig_k in sigma_grid:
            c_k = C / sig_k
            lb_k = [l + 0.01, l + c_k]
            ub_k = [0.6, 3.0]
            bounds_k = [(lb_k[0], ub_k[0]), (lb_k[1], ub_k[1])]
            cons_k = {'type': 'ineq', 'fun': lambda x: -np.array(constraints_fun(x, c_k, l))}

            def obj_k(x):
                mu_k, _ = long_return(x[0], x[1], c_k, l, sig_k, f)
                return -mu_k

            res_k = minimize(obj_k, x0, method='SLSQP', bounds=bounds_k, constraints=cons_k, options={'ftol': 1e-8, 'disp': False})
            if res_k.success:
                d_k = abs(res_k.x[0])
                u_k = res_k.x[1]
                mu_k, f_k = long_return(-d_k, u_k, c_k, l, sig_k, f)
                d_grid.append(d_k)
                u_grid.append(u_k)
                mu_grid.append(mu_k)
                f_grid.append(f_k if not isinstance(f, (int, float)) else np.nan)
            else:
                d_grid.append(np.nan)
                u_grid.append(np.nan)
                mu_grid.append(np.nan)
                f_grid.append(np.nan)

        interp_d = interp1d(sigma_grid, d_grid, bounds_error=False, fill_value="extrapolate")
        interp_u = interp1d(sigma_grid, u_grid, bounds_error=False, fill_value="extrapolate")
        interp_mu = interp1d(sigma_grid, mu_grid, bounds_error=False, fill_value="extrapolate")
        interp_f = interp1d(sigma_grid, f_grid, bounds_error=False, fill_value="extrapolate") if not isinstance(f, (int, float)) else None

        d_boot = interp_d(sigma_vec)
        u_boot = interp_u(sigma_vec)
        mu_boot = interp_mu(sigma_vec) / theta_vec
        if interp_f:
            f_boot = interp_f(sigma_vec)

        R['d_CI'] = list(np.percentile(d_boot, [100 * alpha / 2, 100 * (1 - alpha / 2)]))
        R['u_CI'] = list(np.percentile(u_boot, [100 * alpha / 2, 100 * (1 - alpha / 2)]))
        R['mu_CI'] = list(np.percentile(mu_boot, [100 * alpha / 2, 100 * (1 - alpha / 2)]))
        if interp_f:
            R['f_opt_CI'] = list(np.percentile(f_boot, [100 * alpha / 2, 100 * (1 - alpha / 2)]))
            R['f_input'] = np.nan

    return R

def erfid_matlab(x, y):
    """
    Riproduce la funzione MATLAB:
        f(t) = exp(+t^2/2)
        val   = sqrt(2/pi) * ∫_y^x f(t) dt
    """
    integrand = lambda t: np.exp(t**2 / 2)
    integral_value, _ = quad(integrand, y, x)
    return np.sqrt(2/np.pi) * integral_value


def long_return(d, u, c, l, sigma, f):
    """
    Compute long-run return mu of OU trading strategy.

    Parameters
    ----------
    d : float
        Lower entry level (negative, in sigma-units)
    u : float
        Upper exit level (positive, in sigma-units)
    c : float
        Transaction cost (in sigma-units)
    l : float
        Stop-loss level (in sigma-units, negative)
    sigma : float
        Stationary sigma of OU process
    f : float or str
        Fixed leverage or 'opt'

    Returns
    -------
    mu : float
        Long-run log-return
    fStar : float
        Optimal leverage (if f == 'opt') or input leverage
    """
    if (u - d <= c) or (d <= l) or (u <= d):
        return -np.inf, np.nan

    expoUD = np.exp(sigma * (u - d - c)) - 1
    expoLD = np.exp(sigma * (l - d - c)) - 1

    if isinstance(f, (int, float)):  # Fixed leverage
        fStar = f
    else:
        denomUL = erfid_matlab(u, l)
        fStar = -erfid_matlab(d, l) / (expoLD * denomUL) - erfid_matlab(u, d) / (expoUD * denomUL)

    try:
        mu = (2 / pi) * (
            np.log(1 + fStar * expoUD) / erfid_matlab(u, d) +
            np.log(1 + fStar * expoLD) / erfid_matlab(d, l)
        )
    except:
        mu = -np.inf

    return mu, fStar
