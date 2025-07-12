function [count_plus, count_minus]= count_d_to_u_l(data, d, u,l)
%COUNT_D_TO_U_L  Count completed trades in a band-crossing strategy.
%
%   [count_plus, count_minus] = COUNT_D_TO_U_L(data, d, u, l)
%
%   Counts the number of completed trades based on a threshold-crossing logic:
%
%     A "positive" (profit) trade is counted when the process moves from
%       below d to above u (i.e., d → u).
%
%     A "negative" (loss) trade is counted when the process moves from
%       above d to below l (i.e., d → l), acting as a stop-loss.
%
%   The logic simulates entry/exit rules for mean-reversion strategies with
%   both profit-taking and stop-loss exits.
%
%   INPUTS:
%     data  - Vector of de-meaned spread observations (X_t - eta)
%     d     - Entry level (in process units)
%     u     - Profit-taking level (must be > d)
%     l     - Stop-loss level (must be < d)
%
%   OUTPUTS:
%     count_plus  - Number of completed profitable trades (d → u)
%     count_minus - Number of completed loss trades (d → l)
%
%   NOTES:
%     - States reset when a trade completes or is aborted by re-touching d.
%     - Assumes process is observed at discrete time steps.

    count_plus = 0;
    count_minus = 0;
    state_p = 'u->d';  % 'u->d' or 'd->u'
    state_m = 'l->d';
    for i = 1:length(data)
        val = data(i);

        switch state_p
            case 'u->d'
                if val <= d
                    state_p = 'd->u';
                end

            case 'd->u'
                if val >= u
                    count_plus = count_plus + 1;
                    state_p = 'u->d';
                elseif val <= d
                    % reset the d touch if it happens again
                    state_p = 'd->u';
                end
        end
        
         switch state_m
            case 'l->d'
                if val >= d
                    state_m = 'd->l';
                end

            case 'd->l'
                if val <= l
                    count_minus = count_minus + 1;
                    state_m = 'l->d';
                elseif val >= d
                    % reset the d touch if it happens again
                    state_m = 'd->l';
                end
        end
    end
end