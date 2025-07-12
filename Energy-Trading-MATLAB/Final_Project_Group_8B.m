clc
close all
clear 
warning('off', 'all')
format long
rng(42)

%% Add required folders to path

addpath("Data/");
addpath("Data Processing/");
addpath("Prints and Plots/");
addpath("Trading Strategy & Band Optimization/");
addpath("OU model Estimation & Simulation/");
addpath("Cost & Analysis Utilities/");

%% Load data and build price table

file = "HO-LGO.xlsm";   % Data file
data = readtable(file); % Read Excel data into table

startDate = datetime(2015,04,22);
endDate = datetime(2016,04,22);

T = buildPriceTable(data, 'Var2', ...
        'HOc2','HOc2_1','', '', 42, ...       % product 1 with correction
        'LGOc6','LGOc6_1','', '', 1/7.44, ... % product 2 with correction 
        startDate, endDate);                  % Date range

%% Question a: Compute optimal trading bands for varying cost parameter c

% Setup parameters
M_a = 1;            % No bootstrap needed
l = -1.96;          % Stop-Loss in sigma units
c_sigma_max = 0.76; % Max sigma for the plot
f = 1;              % Leverage

G.parameters.sigma = 1; % Arbitrary for f=1
G.parameters.k = 1;     % Arbitrary for f=1
C_max = c_sigma_max*G.parameters.sigma/sqrt(2*G.parameters.k);

C_values = linspace(0, C_max, 100); % Corrisponding C for a c=(0,0.76) in sigma units

u = zeros(length(C_values),1); % Pre-Allocation
d = zeros(length(C_values),1); % Pre-Allocation 

for i = 1:100
    X = optimal_trading_bands(M_a, l, f, G, C_values(i)); % Compute bands for each c
    u(i) = X.u_estimated;
    d(i) = X.d_estimated;
end

c_values = C_values .* sqrt(2*G.parameters.k) / G.parameters.sigma; % c = (0,0.76) in sigma units
plot_bands(d, u, c_values);

%% Question b
% - Filter IS (9-16), Remove outliers, Statistical bootstrap for OU parameters estimation

% Data Cleaning

startHour_IS = 9;  % Start time window for IS data
endHour_IS = 16;   % End time window for IS data
startHour_OS = 17; % Start time window for OS data
endHour_OS = 20;   % End time window for OS data
splitMonths = 9;   % Month division between IS/OS data
M = 10000;         % Bootstrap samples
alpha = 0.05;      % Confidence level

% Creation of IS and OS dataset
[data_IS_pre_out, data_OS_pre_out] = trimAndSplitPriceTable(T,  ...
                                      startHour_IS, endHour_IS, ...
                                      startHour_OS, endHour_OS, splitMonths);

% Outliers removal and print
[data_IS_9_16, outliers_9_16, isOutlierMask_9_16] = removeOutliersIS(data_IS_pre_out);
plotandprintOutliers(data_IS_pre_out, isOutlierMask_9_16, 'Outliers from IS 9-16 Dataset');

% Statistical Bootstrap for 9-16 time window

% Setup parameters
N_9_16 = height(data_IS_9_16);                                 % Number of sample
t = yearfrac(data_IS_9_16.Time(1), data_IS_9_16.Time(end), 1); % Time 
dt = t / N_9_16;                                               % Delta time

% Bootstrap OU parameters
P_9_16 = ouBootstrap(data_IS_9_16.Rt, dt, M, alpha);
ou_parameters_print(P_9_16);
plot_bootstrap_distributions(P_9_16, M);

% - Filter IS (8-16), Remove outliers, Statistical bootstrap for OU parameters estimation

% Data Cleaning

startHour_IS = 8; % Start time window for IS data

[data_IS_pre_out, ~] = trimAndSplitPriceTable(T, startHour_IS, ...
                        endHour_IS, startHour_OS, endHour_OS, splitMonths);
[data_IS_8_16, outliers_8_16, isOutlierMask_8_16] = removeOutliersIS(data_IS_pre_out);
plotandprintOutliers(data_IS_pre_out, isOutlierMask_8_16, 'Outliers from IS 8-16 Dataset');

% Statistical Bootstrap for 8-16 time window

% Setup parameters
N_8_16 = height(data_IS_8_16);                                 % Number of sample
t = yearfrac(data_IS_8_16.Time(1), data_IS_8_16.Time(end), 1); % Time
dt = t / N_8_16;                                               % Delta time
% Bootstrap OU parameters
P_8_16 = ouBootstrap(data_IS_8_16.Rt, dt, M, alpha);
ou_parameters_print(P_8_16);
plot_bootstrap_distributions(P_8_16, M);

% Removing outliers in OS (excluded 17-20) dataset

[data_OS, outliers_OS, isOutlierMask_OS] = removeOutliersIS(data_OS_pre_out);
plotandprintOutliers(data_OS_pre_out, isOutlierMask_OS, 'Outliers from OS excluded 17-20 Dataset');

%% Question c-e: Compute transaction costs and apply trading strategy for every leverage value and stop-loss

C = transaction_cost(data_IS_9_16);

% Setup parameters
l_list = {-1.282, -1.645, -1.960, -2.236}; % Stop-Loss level
f_list = {1, 2, 5, 'opt'};                 % Leverage values ('opt' computes optimal value of leverage)

% Computing optimal bands and optimal OS return 
[optimal_bands, mu_OS] = trading_strategy(P_9_16, f_list, l_list, M, alpha, C, data_OS);

%Plot the return with respect to different leverage for every stop-loss in l_list
mu_leverage_plot(l_list, 40, P_9_16, optimal_bands, C, alpha); 

%% Question g: Computing same analysis of before for Governative-Futures pairs

file_g = "GovernativeFutures.xlsx";

data_g = readtable(file_g);

% Setup parameters
splitMonths = 4; % For dataset division
M = 10000;       % Bootstrap samples
alpha = 0.05;    % Confidence level

l_list = {-1.645};
f_list = {'opt'};

%% Analyzing a couple IKA-OATA  

product1 = 'IKA';
product2 = 'OATA';

fprintf('Analyzing optimal trading strategy for %s vs %s\n',product1, product2)

T_IKA_OATA = buildPriceTable(data_g, 'timestamp', ...
        '','',product1,0.01 ,1, ...  % product 1 IKA with just Mid and tick and conversion=1
        '','',product2,0.01, 1, ...  % product 2 OATA with just Mid and tick and conversion=1
        '','');                      % Date range (every date)

plot_mid_and_log_ratio(T_IKA_OATA);

% Computing optimal analysis 
IKA_OATA = pair_trading_analysis(T_IKA_OATA, splitMonths, M, alpha, l_list, f_list, 100);

%% Analyzing a couple OEA-OATA 

product1 = 'OEA';
product2 = 'OATA';

fprintf('Analyzing optimal trading strategy for %s vs %s\n',product1, product2)

T_OEA_OATA = buildPriceTable(data_g, 'timestamp', ...
        '','',product1, 0.01 ,1, ... % product 1 IKA with just Mid and tick and conversion=1
        '','',product2, 0.01, 1, ... % product 2 OATA with just Mid and tick and conversion=1
        '','');                      % Date range (every date)

plot_mid_and_log_ratio(T_OEA_OATA);

% Computing optimal analysis 
OEA_OATA = pair_trading_analysis(T_OEA_OATA, splitMonths, M, alpha, l_list, f_list, 100);

%% Analyzing a couple IKA-RXA

product1 = 'IKA';
product2 = 'RXA';

fprintf('Analyzing optimal trading strategy for %s vs %s\n',product1, product2)

T_IKA_RXA = buildPriceTable(data_g, 'timestamp', ...
        '','',product1,0.01 ,1, ...  % product 1 IKA with just Mid and tick and conversion=1
        '','',product2,0.01, 1, ...  % product 2 OATA with just Mid and tick and conversion=1
        '','');                      % Date range (every date)

plot_mid_and_log_ratio(T_IKA_RXA);

% Computing optimal analysis 
IKA_RXA = pair_trading_analysis(T_IKA_RXA, splitMonths, M, alpha, l_list, f_list, 150);