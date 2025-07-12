function T = buildPriceTable( ...
        data           , ...
        timeCol        , ...
        bid1Col        , ask1Col , mid1Col , tick1 , conv1 , ...
        bid2Col        , ask2Col , mid2Col , tick2 , conv2 , ...
        startDate      , endDate)
%BUILDPRICETABLE  Preprocess and standardize price data for two assets.
%
%   T = BUILDPRICETABLE(data, timeCol, ...
%                       bid1Col, ask1Col, mid1Col, tick1, conv1, ...
%                       bid2Col, ask2Col, mid2Col, tick2, conv2, ...
%                       startDate, endDate)
%
%   Returns a table T containing the time series of bid, ask, and mid
%   prices for two products, along with their log-price spread.
%
%   INPUTS:
%     data      - A table containing raw price data.
%     timeCol   - Name of the datetime column.
%
%     bid1Col   - Name of the bid price column for product 1 (can be empty).
%     ask1Col   - Name of the ask price column for product 1 (can be empty).
%     mid1Col   - Name of the mid price column for product 1 (used if bid/ask missing).
%     tick1     - Minimum price increment for product 1.
%     conv1     - Price scaling factor for product 1.
%
%     bid2Col   - Name of the bid price column for product 2 (can be empty).
%     ask2Col   - Name of the ask price column for product 2 (can be empty).
%     mid2Col   - Name of the mid price column for product 2 (used if bid/ask missing).
%     tick2     - Minimum price increment for product 2.
%     conv2     - Price scaling factor for product 2.
%
%     startDate - (Optional) Lower bound for time filtering (inclusive).
%     endDate   - (Optional) Upper bound for time filtering (exclusive).
%
%   OUTPUT:
%     T         - A table with the following columns:
%                  - Time : datetime
%                  - Bid1, Ask1, Mid1 : bid/ask/mid prices for product 1
%                  - Bid2, Ask2, Mid2 : bid/ask/mid prices for product 2
%                  - Rt : log spread between mid prices (log(Mid1 / Mid2))
%
%   If bid and ask prices are not available for a product, the function
%   estimates them using the mid price and tick size.

    % Trim by date if bounds are provided
    if exist('startDate','var') && exist('endDate','var') && ...
       ~isempty(startDate) && ~isempty(endDate)
        data = data(data.(timeCol) >= startDate & data.(timeCol) < endDate, :);
    end
    
    % Extract / prepare time
    Time = data.(timeCol);
    
    % Product 1: Bid / Ask / Mid
    haveBA1 = ~isempty(bid1Col) && ~isempty(ask1Col) && ...
              any(data.(bid1Col)~=0) && any(data.(ask1Col)~=0);
    
    if haveBA1
        Bid1 = data.(bid1Col) * conv1;
        Ask1 = data.(ask1Col) * conv1;
        Mid1 = (Bid1 + Ask1) / 2;
    else
        if isempty(mid1Col)
            error('Need either Bid+Ask or Mid for product 1.');
        end
        Mid1 = data.(mid1Col) * conv1;
        Bid1 = Mid1 - tick1/2;
        Ask1 = Mid1 + tick1/2;
    end
    
    % Product 2
    haveBA2 = ~isempty(bid2Col) && ~isempty(ask2Col) && ...
              any(data.(bid2Col)~=0) && any(data.(ask2Col)~=0);
    
    if haveBA2
        Bid2 = data.(bid2Col) * conv2;
        Ask2 = data.(ask2Col) * conv2;
        Mid2 = (Bid2 + Ask2) / 2;
    else
        if isempty(mid2Col)
            error('Need either Bid+Ask or Mid for product 2.');
        end
        Mid2 = data.(mid2Col) * conv2;
        Bid2 = Mid2 - tick2/2;
        Ask2 = Mid2 + tick2/2;
    end
    
    % Log-spread
    Rt = log(Mid1 ./ Mid2);
    
    % Assemble output
    T = table(Time, Bid1, Ask1, Mid1, Bid2, Ask2, Mid2, Rt, ...
              'VariableNames', {'Time','Bid1','Ask1','Mid1', ...
                                       'Bid2','Ask2','Mid2','Rt'});
end
