function [k,eta,sigma] = ou_mle(x,dt)
%OU_MLE  Maximum Likelihood Estimation for the Ornstein–Uhlenbeck process.
%
%   [k, eta, sigma] = OU_MLE(x, dt)
%
%   Estimates the parameters of the continuous-time Ornstein–Uhlenbeck (OU)
%   process defined by the stochastic differential equation:
%
%       dX_t = k (eta − X_t) dt + sigma dW_t
%
%   using closed-form maximum likelihood estimators based on discrete-time
%   observations.
%
%   INPUT:
%     x   - Column vector of observations (X_0, X_1, ..., X_N)
%     dt  - Time increment between consecutive observations
%
%   OUTPUT:
%     k      - Mean-reversion rate
%     eta    - Long-term mean level
%     sigma  - Volatility coefficient
%
%   NOTE:
%     This implementation assumes constant time steps and is based on
%     moment matching and exact transition density of the OU process.

    x_plus  = x(2:end);        
    x_minus = x(1:end-1);      
    N       = length(x_minus);
    
    Y_m  = mean(x_minus);
    Y_p  = mean(x_plus);
    Y_mm = mean(x_minus.^2);
    Y_pp = mean(x_plus.^2);
    Y_pm = mean(x_minus .* x_plus);
    
    
    rho = (Y_pm - Y_m*Y_p) / (Y_mm - Y_m^2);
    
    k   = -log(rho) / dt;
    
    eta = Y_p + ((x(end)-x(1))/N) * ...
          (Y_pm - Y_m*Y_p) / ((Y_mm - Y_m^2) - (Y_pm - Y_m*Y_p));
    
    sigma2 = Y_pp - Y_p^2 - (Y_pm - Y_m*Y_p)^2 / (Y_mm - Y_m^2);
    
    sigma  = sqrt( (2*k*sigma2) / (1 - exp(-2*k*dt)) ); 

end