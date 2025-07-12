function C=transaction_cost(data)
%TRANSACTION_COST  Compute average log transaction cost for a two-asset pair.
%
%   C = TRANSACTION_COST(data)
%
%   Computes the mean of the total log transaction cost across time, based on 
%   bid-ask spreads of two traded instruments.
%
%   INPUT:
%     data - Table containing the columns:
%              Bid1, Ask1 : bid and ask prices of product 1
%              Bid2, Ask2 : bid and ask prices of product 2
%
%   OUTPUT:
%     C - Average total log transaction cost over the sample
%
%   FORMULA:
%     C_t = log(Ask1 / Bid1) + log(Ask2 / Bid2)
%     C   = mean(C_t)

    C_t = log(data.Ask1 ./ data.Bid1) + log(data.Ask2 ./ data.Bid2);
    C = mean(C_t); % Average log transaction cost

end