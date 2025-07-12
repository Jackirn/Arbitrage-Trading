function [mu, fStar] = long_return(d, u, c, l, sigma, f)
%LONG_RETURN  In‐sample long‐run return for OU bands (with optional leverage)
%
%   [mu, fStar] = long_return(d,u,c,l,sigma,f)
%
%   Inputs:
%     d, u    - lower/upper band (σ-units)
%     c       - cost (σ-units)
%     l       - stop-loss level (negative σ)
%     sigma   - stationary volatility
%     f       - numeric leverage or 'opt'
%
%   Outputs:
%     mu      - long‐run return (zero or negative means no trade)
%     fStar   - optimal f if 'opt', otherwise same as input f

    % 1) Feasibility check
    if (u - d <= c) || (d <= l) || (u <= d)
        mu    = -Inf;
        fStar = NaN;
        return;
    end

    c_bar = (l - d) ...
     + log( 1 + Erfid(d, l)/Erfid(u, l) .* ( exp( sigma * (u - l) ) - 1 ) )/sigma;
    

    % 2) Exponential terms
    expoUD = exp( sigma .* (u - d - c) ) - 1;   % scalar
    expoLD = exp( sigma .* (l - d - c) ) - 1;   % scalar

    % 3) Leverage
    if isnumeric(f)
        fStar = f;
    else

        % derive optimal f via closed‐form (all scalars here)
        denomUL = Erfid(u, l);
        fStar   = - Erfid(d, l) ./ ( expoLD .* denomUL ) ...
                  - Erfid(u, d) ./ ( expoUD .* denomUL );
        if c>c_bar
            fStar = 0;
         end


    end


    % 4) Long‐run return (scalar)
    term1 = log( 1 + fStar .* expoUD ) ./ Erfid(u, d);
    term2 = log( 1 + fStar .* expoLD ) ./ Erfid(d, l);

    mu = (2/pi) .* ( term1 + term2 );
end
