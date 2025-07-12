function T = buildCryptoPriceTable(file1, file2, ...
                                    tick1, conv1, ...
                                    tick2, conv2, ...
                                    startDate, endDate)
% Robust version that renames 'Open time' to 'Time'

    % Read and rename variables in file1
    data1 = readtable(file1, 'PreserveVariableNames', true);
    if ismember('Open time', data1.Properties.VariableNames)
        data1.Properties.VariableNames{strcmp(data1.Properties.VariableNames, 'Open time')} = 'Time';
    end

    % Read and rename variables in file2
    data2 = readtable(file2, 'PreserveVariableNames', true);
    if ismember('Open time', data2.Properties.VariableNames)
        data2.Properties.VariableNames{strcmp(data2.Properties.VariableNames, 'Open time')} = 'Time';
    end

    % Convert Time column to datetime if it's not already
    if ~isdatetime(data1.Time)
        data1.Time = datetime(data1.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    end
    if ~isdatetime(data2.Time)
        data2.Time = datetime(data2.Time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    end

    % Extract time and prices
    Mid1  = data1.Close * conv1;
    Bid1  = Mid1 - tick1/2;
    Ask1  = Mid1 + tick1/2;
    T1 = table(data1.Time, Bid1, Ask1, Mid1, 'VariableNames', {'Time', 'Bid1', 'Ask1', 'Mid1'});

    Mid2  = data2.Close * conv2;
    Bid2  = Mid2 - tick2/2;
    Ask2  = Mid2 + tick2/2;
    T2 = table(data2.Time, Bid2, Ask2, Mid2, 'VariableNames', {'Time', 'Bid2', 'Ask2', 'Mid2'});

    % Inner join on Time
    T = innerjoin(T1, T2, 'Keys', 'Time');

    % Optional time filtering
    if nargin >= 7 && ~isempty(startDate)
        T = T(T.Time >= startDate, :);
    end
    if nargin >= 8 && ~isempty(endDate)
        T = T(T.Time < endDate, :);
    end

    % Compute log-price spread
    T.Rt = log(T.Mid1 ./ T.Mid2);
end
