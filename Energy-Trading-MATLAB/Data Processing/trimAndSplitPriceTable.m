function [data_IS, data_OS] = trimAndSplitPriceTable(data, ...
                                   IS_startHour, IS_endHour, ...
                                   OS_startHour, OS_endHour, ...
                                   splitMonths)
%TRIMANDSPLITPRICETABLE
% Filters data by IS hour window and OS hour exclusion, then splits by months.
%
% INPUTS:
%   data          - Table with datetime column 'Time' (sorted or unsorted)
%   IS_startHour  - Start hour for in-sample (inclusive) or [] to keep all
%   IS_endHour    - End hour for in-sample (inclusive) or [] to keep all
%   OS_startHour  - Start hour for OOS exclusion window or [] to keep all
%   OS_endHour    - End hour for OOS exclusion window or [] to keep all
%   splitMonths   - Number of months for IS period (remaining is OS)
%
% OUTPUTS:
%   data_IS       - In-sample data filtered within IS hour window (or all)
%   data_OS       - Out-of-sample data filtered to exclude OS hour window (or all)

    % Sort by time (safety)
    data = sortrows(data, 'Time');

    % Split date
    splitDate = data.Time(1) + calmonths(splitMonths);

    % Split into raw IS and OS
    raw_IS = data(data.Time < splitDate, :);
    raw_OS = data(data.Time >= splitDate, :);

    % Filter IS if needed
    if ~isempty(IS_startHour) && ~isempty(IS_endHour)
        time_IS = hour(raw_IS.Time) + minute(raw_IS.Time)/60;
        mask_IS = (time_IS >= IS_startHour) & (time_IS <= IS_endHour);
        data_IS = raw_IS(mask_IS, :);
    else
        data_IS = raw_IS;  % keep all
    end

    % Filter OS if needed
    if ~isempty(OS_startHour) && ~isempty(OS_endHour)
        time_OS = hour(raw_OS.Time) + minute(raw_OS.Time)/60;
        mask_OS = (time_OS <= OS_startHour) | (time_OS >= OS_endHour);
        data_OS = raw_OS(mask_OS, :);
    else
        data_OS = raw_OS;  % keep all
    end
end
