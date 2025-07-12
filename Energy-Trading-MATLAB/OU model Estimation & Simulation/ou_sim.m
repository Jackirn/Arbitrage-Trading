function x = ou_sim(x0,k,eta,sigma,dt,N)
%OU_SIM  Simulate a path from the Ornstein-Uhlenbeck (OU) process.
%
%   x = OU_SIM(x0, k, eta, sigma, dt, N)
%
%   Simulates a single path of the Ornstein-Uhlenbeck process defined by:
%
%       dX_t = k * (eta - X_t) * dt + sigma * dW_t
%
%   using the exact solution over a uniform time grid with step dt.
%
%   Exact one-step update formula:
%
%       X_{t+dt} = eta + (X_t - eta) * exp(-k * dt)
%                  + sigma * sqrt((1 - exp(-2k * dt)) / (2k)) * Z,
%         where Z ~ N(0,1)
%
%   INPUTS:
%     x0     - Initial value of the process
%     k      - Mean-reversion speed (k > 0)
%     eta    - Long-run mean level
%     sigma  - Volatility coefficient (sigma > 0)
%     dt     - Time step between observations (dt > 0)
%     N      - Number of time steps (returns a path of length N+1)
%
%   OUTPUT:
%     x      - (N+1)×1 vector representing the simulated OU path,
%              with x(1) = x0

    x          = zeros(N+1,1);             % pre-allocate output
    x(1)       = x0;                       % set starting value
    
    z          = randn(N,1);               % 
    a  = exp(-k*dt);                       % e^{-kΔt}
    b  = eta * (1 - a);                    % drift term  eta(1-e^{-kΔt})
    sd = sigma * sqrt((1 - a^2) / (2*k));  % conditional std-dev
    
    for i = 1:N
        % OU exact simulation
        x(i+1) = a * x(i) + b + sd * z(i);
    end
end