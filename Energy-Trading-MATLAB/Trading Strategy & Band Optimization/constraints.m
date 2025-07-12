function [cineq, ceq] = constraints(x, c, l)
% CONSTRAINTS Define inequality constraints for optimal trading bands
%
% Inputs:
%   x = [d, u] : entry and exit levels in standard units
%   c         : transaction cost (in same units)
%   l         : stop-loss level (fixed, in standard units)
%
% Outputs:
%   cineq : inequality constraints (must be <= 0)
%   ceq   : equality constraints (none used here)

    d = x(1);
    u = x(2);

    buffer = 0.1;  % Minimum additional gap between bands and transaction cost

    % Constraint 1: Ensure u - d > c + buffer
    % Rewritten as: c + buffer - (u - d) <= 0
    c1 = c + buffer - (u - d);

    % Constraint 2: Ensure d > l
    % Rewritten as: l - d <= 0
    c2 = l - d;

    % Constraint 3: Ensure u > d
    % Rewritten as: d - u <= 0
    c3 = d - u;

    % Return constraints as column vector
    cineq = [c1; c2; c3];
    ceq = [];  % No equality constraints
end