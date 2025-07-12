
import pandas as pd
import numpy as np
from scipy.stats import iqr
import itertools
from statsmodels.tsa.stattools import adfuller


def build_multiindex_header(df_raw, row1=0, row2=1):
    # Extract the two rows for the header
    first_row = df_raw.iloc[row1].fillna(method='ffill').astype(str)
    second_row = df_raw.iloc[row2].fillna('').astype(str)
    combined = first_row + '_' + second_row
    return combined


def clean_numeric_columns(df, cols):
    for col in cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    # Optionally drop rows with NaNs in any of these columns
    df.dropna(subset=cols, inplace=True)
    return df


def build_price_table(
    data,
    time_col,
    bid1_col=None, ask1_col=None, mid1_col=None, tick1=None, conv1=1,
    bid2_col=None, ask2_col=None, mid2_col=None, tick2=None, conv2=1,
    start_date=None, end_date=None,
):
    # Optional date filter
    if start_date is not None and end_date is not None:
        data = data[(data[time_col] >= pd.to_datetime(start_date)) & 
                    (data[time_col] < pd.to_datetime(end_date))]

    time = data[time_col]

    # PRODUCT 1
    have_bid_ask1 = bid1_col in data.columns and ask1_col in data.columns \
                    and data[bid1_col].notna().any() and data[ask1_col].notna().any()

    if have_bid_ask1:
        bid1 = pd.to_numeric(data[bid1_col], errors='coerce') * conv1
        ask1 = pd.to_numeric(data[ask1_col], errors='coerce') * conv1
        mid1 = (bid1 + ask1) / 2
        print(type(mid1))
    else:
        if mid1_col is None or tick1 is None:
            raise ValueError("If Bid/Ask missing for product 1, mid1_col and tick1 must be provided")
        mid1 = pd.to_numeric(data[mid1_col], errors='coerce') * conv1
        bid1 = mid1 - tick1 / 2
        ask1 = mid1 + tick1 / 2

    # PRODUCT 2
    have_bid_ask2 = bid2_col in data.columns and ask2_col in data.columns \
                    and data[bid2_col].notna().any() and data[ask2_col].notna().any()

    if have_bid_ask2:
        bid2 = pd.to_numeric(data[bid2_col], errors='coerce') * conv2
        ask2 = pd.to_numeric(data[ask2_col], errors='coerce') * conv2
        mid2 = (bid2 + ask2) / 2
    else:
        if mid2_col is None or tick2 is None:
            raise ValueError("If Bid/Ask missing for product 2, mid2_col and tick2 must be provided")
        mid2 = pd.to_numeric(data[mid2_col], errors='coerce') * conv2
        bid2 = mid2 - tick2 / 2
        ask2 = mid2 + tick2 / 2

    # Compute log spread (Rt)
    rt = np.log(mid1 / mid2)

    # Assemble output DataFrame
    out = pd.DataFrame({
        'Time': time,
        'Bid1': bid1, 'Ask1': ask1, 'Mid1': mid1,
        'Bid2': bid2, 'Ask2': ask2, 'Mid2': mid2,
        'Rt': rt
    })

    return out


def load_and_process_price_data(filepath, time_col, bid_ask_cols, mid_cols=None, ticks=None, convs=None, start_date=None, end_date=None):
    import pandas as pd

    df_raw = pd.read_excel(filepath)
    header = build_multiindex_header(df_raw)
    
    df = df_raw.iloc[2:].copy()
    df.columns = header

    df.reset_index(drop=True, inplace=True)

    # Clean numeric columns (bid/ask columns)
    all_cols = list(bid_ask_cols)
    if mid_cols:
        all_cols += list(mid_cols)
    df = clean_numeric_columns(df, all_cols)

    # Prepare arguments for build_price_table, assuming function signature:
    # build_price_table(data, time_col, bid1_col, ask1_col, mid1_col, tick1, conv1, bid2_col, ask2_col, mid2_col, tick2, conv2, start_date, end_date)

    args = dict(
        data=df,
        time_col=time_col,
        bid1_col=bid_ask_cols[0], ask1_col=bid_ask_cols[1],
        mid1_col=None if mid_cols is None else mid_cols[0],
        tick1=0 if ticks is None else ticks[0],
        conv1=1 if convs is None else convs[0],
        bid2_col=bid_ask_cols[2], ask2_col=bid_ask_cols[3],
        mid2_col=None if mid_cols is None else mid_cols[1],
        tick2=0 if ticks is None else ticks[1],
        conv2=1 if convs is None else convs[1],
        start_date=start_date,
        end_date=end_date
    )

    return build_price_table(**args)



def trim_and_split_price_table(
    data: pd.DataFrame,
    IS_start_hour: float | None = None,
    IS_end_hour: float | None = None,
    OS_start_hour: float | None = None,
    OS_end_hour: float | None = None,
    split_months: int = 0
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filters data by IS hour window and OS hour exclusion, then splits by months.

    Parameters:
    - data: pd.DataFrame with a datetime column 'Time'
    - IS_start_hour: start hour for in-sample (inclusive), None to keep all
    - IS_end_hour: end hour for in-sample (inclusive), None to keep all
    - OS_start_hour: start hour for out-of-sample exclusion window, None to keep all
    - OS_end_hour: end hour for out-of-sample exclusion window, None to keep all
    - split_months: number of months for in-sample period (remaining is out-of-sample)

    Returns:
    - data_IS: filtered in-sample DataFrame
    - data_OS: filtered out-of-sample DataFrame
    """

    data['Time'] = pd.to_datetime(data['Time'])

    # Sort by time
    data = data.sort_values('Time').reset_index(drop=True)

    # Compute split date
    split_date = data['Time'].iloc[0] + pd.DateOffset(months=split_months)

    # Split data by date
    raw_IS = data[data['Time'] < split_date]
    raw_OS = data[data['Time'] >= split_date]

    # Filter IS if needed
    if IS_start_hour is not None and IS_end_hour is not None:
        time_IS = raw_IS['Time'].dt.hour + raw_IS['Time'].dt.minute / 60
        mask_IS = (time_IS >= IS_start_hour) & (time_IS <= IS_end_hour)
        data_IS = raw_IS[mask_IS]
    else:
        data_IS = raw_IS

    # Filter OS if needed
    if OS_start_hour is not None and OS_end_hour is not None:
        time_OS = raw_OS['Time'].dt.hour + raw_OS['Time'].dt.minute / 60
        # keep rows outside the exclusion window
        mask_OS = (time_OS <= OS_start_hour) | (time_OS >= OS_end_hour)
        data_OS = raw_OS[mask_OS]
    else:
        data_OS = raw_OS

    return data_IS.reset_index(drop=True), data_OS.reset_index(drop=True)



def split_price_table_by_months(
    data: pd.DataFrame,
    split_months: int = 0
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits the data into in-sample (IS) and out-of-sample (OS) based on a month offset.

    Parameters:
    - data: pd.DataFrame with a datetime column 'Time'
    - split_months: number of months for the in-sample period; remaining is out-of-sample

    Returns:
    - data_IS: in-sample DataFrame (first `split_months` of data)
    - data_OS: out-of-sample DataFrame (remaining data)
    """

    data['Time'] = pd.to_datetime(data['Time'])
    data = data.sort_values('Time').reset_index(drop=True)

    # Compute the split point
    split_date = data['Time'].iloc[0] + pd.DateOffset(months=split_months)

    # Perform the split
    data_IS = data[data['Time'] < split_date].reset_index(drop=True)
    data_OS = data[data['Time'] >= split_date].reset_index(drop=True)

    return data_IS, data_OS


def filter_log_spread_outliers(data):
    """
    Remove log-spread (Rt) outliers from a price DataFrame.

    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with at least the 'Rt' column.

    Returns:
    --------
    clean_data : pd.DataFrame
        DataFrame filtered to exclude Rt outliers outside 3*IQR.
    outlier_mask : pd.Series (bool)
        Boolean mask identifying outliers in the original data.
    outlier_values : pd.DataFrame
        Rows corresponding to outliers.
    """
    Rt = data['Rt']

    Q1 = Rt.quantile(0.25)
    Q3 = Rt.quantile(0.75)
    iqr_val = Q3 - Q1
    lo_bound = Q1 - 3 * iqr_val
    hi_bound = Q3 + 3 * iqr_val

    outlier_mask = (Rt < lo_bound) | (Rt > hi_bound)
    outlier_values = data.loc[outlier_mask]
    clean_data = data.loc[~outlier_mask]

    return clean_data, outlier_mask, outlier_values



def filter_antipersistent_outliers(data):
    """
    Removes antipersistent outliers from Rt time series.

    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame with a column 'Rt'.

    Returns:
    --------
    clean_data : pd.DataFrame
        Data with outliers removed.
    is_outlier : pd.Series
        Boolean Series marking outliers (with same index as input).
    outlier_values : pd.DataFrame
        Rows in original data that were outliers.
    """
    Rt = data['Rt']
    n = len(data)

    # Compute IQR
    IQR_val = Rt.quantile(0.75) - Rt.quantile(0.25)

    # Initialize logical index for outliers
    is_outlier = np.zeros(n, dtype=bool)

    # Apply the 3-point rule
    for t in range(1, n - 1):
        delta_prev = abs(Rt.iloc[t] - Rt.iloc[t - 1])
        delta_next = abs(Rt.iloc[t + 1] - Rt.iloc[t])

        if delta_prev > IQR_val and delta_next > 0.95 * IQR_val:
            is_outlier[t] = True

    # Convert to Series with original index
    is_outlier = pd.Series(is_outlier, index=data.index)

    # Extract outlier rows
    outlier_values = data[is_outlier]

    # Remove outliers
    clean_data = data[~is_outlier]

    return clean_data, is_outlier, outlier_values






def remove_outliers(data):
    """
    Removes log-spread and antipersistent outliers from IS data.

    Parameters:
    -----------
    data_IS : pd.DataFrame
        DataFrame with a column 'Rt' (log-spread).

    Returns:
    --------
    clean_IS : pd.DataFrame
        Data after removing outliers.
    outliers : pd.DataFrame
        Data rows identified as outliers.
    is_outlier : pd.Series or np.ndarray (boolean)
        Boolean mask of outliers for original data_IS.
    """

    # Step 1: Remove log spread IQR outliers
    data_no_spread_outliers, is_outlier_log, _ = filter_log_spread_outliers(data)
    is_outlier_log = pd.Series(is_outlier_log, index=data.index)

    # Step 2: Remove antipersistent outliers
    clean, is_outlier_anti, _ = filter_antipersistent_outliers(data_no_spread_outliers)
    is_outlier_anti = pd.Series(is_outlier_anti, index=data_no_spread_outliers.index)


    # Combine outlier masks
    # Note: Since clean_IS is subset of data_no_spread_outliers which is subset of data_IS,
    # we reconstruct mask for full data_IS.

    # Initialize full mask with False
    is_outlier = pd.Series(False, index=data.index)

    # Mark log spread outliers True
    is_outlier.loc[is_outlier_log.index] = is_outlier_log

    # Mark antipersistent outliers True
    is_outlier.loc[is_outlier_anti.index] = is_outlier_anti

    # Get outlier rows from original data_IS
    outliers = data[is_outlier]

    return clean, outliers, is_outlier




def load_and_process_crypto_data(file1, file2, time_col='Open time', ticks=None, convs=None,
                                 start_date=None, end_date=None, add_spread_col=True):
    """
    Loads and preprocesses two crypto CSVs with OHLCV 1-minute data.
    Returns a DataFrame ready for OU spread analysis.
    """

    def preprocess_crypto_file(filepath, label):
        df = pd.read_csv(filepath, parse_dates=[time_col])
        df = df[[time_col, 'Close']].copy()
        df.columns = [time_col, f'{label}_mid']
        df[f'{label}_mid'] = df[f'{label}_mid'].astype(float)
        df[f'{label}_bid'] = df[f'{label}_mid'] * (1 - 0.0001)
        df[f'{label}_ask'] = df[f'{label}_mid'] * (1 + 0.0001)
        return df

    # Preprocess both files
    df1 = preprocess_crypto_file(file1, 'A')
    df2 = preprocess_crypto_file(file2, 'B')

    # Merge on time
    df = pd.merge(df1, df2, on=time_col, how='inner')

    # Apply date filtering
    if start_date:
        df = df[df[time_col] >= pd.to_datetime(start_date)]
    if end_date:
        df = df[df[time_col] <= pd.to_datetime(end_date)]

    # Build final price table
    args = dict(
        data=df,
        time_col=time_col,
        bid1_col='A_bid', ask1_col='A_ask', mid1_col='A_mid',
        tick1=0 if ticks is None else ticks[0],
        conv1=1 if convs is None else convs[0],
        bid2_col='B_bid', ask2_col='B_ask', mid2_col='B_mid',
        tick2=0 if ticks is None else ticks[1],
        conv2=1 if convs is None else convs[1],
        start_date=start_date,
        end_date=end_date
    )

    df_output = build_price_table(**args)

    # Add spread (log-price ratio) if needed
    if add_spread_col:
        df_output['Rt'] = np.log(df_output['Mid1'] / df_output['Mid2'])

    return df_output


def analyze_all_pairs(prices_df, spread_type='log', adf_cutoff=0.05):
    """
    Analyzes all 45 unique pairs in a price DataFrame.
    Returns a DataFrame with ADF p-values and whether each pair is potentially stationary.
    """
    tickers = prices_df.columns
    results = []

    for t1, t2 in itertools.combinations(tickers, 2):
        df_pair = prices_df[[t1, t2]].dropna()
        if df_pair.empty:
            continue

        # Compute spread
        if spread_type == 'log':
            spread = np.log(df_pair[t1]) - np.log(df_pair[t2])
        else:
            spread = df_pair[t1] - df_pair[t2]  # simple price difference

        # ADF test
        try:
            adf_stat, p_value, *_ = adfuller(spread)
        except Exception as e:
            p_value = np.nan  # fallback in case of error

        results.append({
            'Ticker1': t1,
            'Ticker2': t2,
            'ADF_pvalue': p_value,
            'Stationary': p_value < adf_cutoff if pd.notnull(p_value) else False
        })

    return pd.DataFrame(results).sort_values('ADF_pvalue')



def process_adf_equity_pairs(prices_df, adf_results_df,
                             spread_type='log', tick_size=0.005,
                             start_date=None, end_date=None,
                             convs=None):
    """
    Processes ADF-true equity pairs and returns a dictionary of spread dataframes.
    
    Each key is a 'TICKER1_TICKER2' string and value is a processed DataFrame
    ready for OU analysis (like the crypto pipeline).
    """
    output_dict = {}

    # Loop through stationary pairs
    true_pairs = adf_results_df[adf_results_df['Stationary']]
    for _, row in true_pairs.iterrows():
        t1, t2 = row['Ticker1'], row['Ticker2']
        df = prices_df[[t1, t2]].dropna().copy()

        # Filter dates
        if start_date:
            df = df[df.index >= pd.to_datetime(start_date)]
        if end_date:
            df = df[df.index <= pd.to_datetime(end_date)]
        if df.empty:
            continue

        # Mid prices
        df_output = pd.DataFrame({
            'Time': df.index,
            'Mid1': df[t1],
            'Mid2': df[t2]
        })

        # Apply conversion if provided
        if convs:
            df_output['Mid1'] *= convs.get(t1, 1)
            df_output['Mid2'] *= convs.get(t2, 1)

        # Bid/Ask (symmetric spread around mid)
        df_output['Bid1'] = df_output['Mid1'] * (1 - tick_size)
        df_output['Ask1'] = df_output['Mid1'] * (1 + tick_size)
        df_output['Bid2'] = df_output['Mid2'] * (1 - tick_size)
        df_output['Ask2'] = df_output['Mid2'] * (1 + tick_size)

        # Add spread
        if spread_type == 'log':
            df_output['Rt'] = np.log(df_output['Mid1'] / df_output['Mid2'])
        else:
            df_output['Rt'] = df_output['Mid1'] - df_output['Mid2']

        # Store
        pair_label = f'{t1}_{t2}'
        output_dict[pair_label] = df_output

    return output_dict
