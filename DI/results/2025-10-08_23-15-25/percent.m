% =========================================================================
% MATLAB Script to Calculate Safe Set Percentages (Robust Version)
% =========================================================================
clear;
clc;
close all;

%% --- 1. Define Files and Parameters ---
results_dir = './'; % Assumes script is in the same directory as the CSVs

% --- AVR Data ---
avr_file = 'AVR_gain_map.csv';

% --- MDR Data ---
mdr_files = {
    'MDR_Z_map_lambda_0.0.csv', ...
    'MDR_Z_map_lambda_0.01.csv', ...
    'MDR_Z_map_lambda_0.1.csv', ...
    'MDR_Z_map_lambda_0.2.csv'
};
lambda_values = [0.0, 0.01, 0.015, 0.02];

%% --- 2. Perform Calculations ---

% --- AVR Calculation ---
try
    % readtable is more robust to formatting errors like trailing commas
    avr_table = readtable(fullfile(results_dir, avr_file), 'ReadVariableNames', false);
    data_vector = table2array(avr_table); % Convert table to a matrix
    data_vector = data_vector(:); % Flatten into a single column vector
    
    % The safe set is where the gain g(s) is approximately 1
    safe_count = sum(data_vector >=1.0)
    total_count = numel(data_vector);
    avr_percentage = (safe_count / total_count) * 100;
catch ME
    warning('Could not process %s. Setting percentage to NaN.', avr_file);
    avr_percentage = NaN;
end

% --- MDR Calculation ---
mdr_percentages = NaN(size(lambda_values)); % Pre-allocate results array
for i = 1:length(mdr_files)
    current_file = mdr_files{i};
    try
        mdr_table = readtable(fullfile(results_dir, current_file), 'ReadVariableNames', false);
        data_vector = table2array(mdr_table);
        data_vector = data_vector(:);

        % The safe set is where the value Z(x) is non-negative
        safe_count = sum(data_vector >=0.);
        total_count = numel(data_vector);
        mdr_percentages(i) = (safe_count / total_count) * 100;
    catch ME
        warning('Could not process %s. Skipping.', current_file);
        % The value will remain NaN
    end
end

%% --- 3. Display Results ---
fprintf('\n======================================================\n');
fprintf('      Safe Set Percentage of Total State Space      \n');
fprintf('======================================================\n');
fprintf('Method                     | Safe Set Percentage\n');
fprintf('---------------------------|----------------------\n');
if ~isnan(avr_percentage)
    fprintf('AVR (g >= 1)               | %.3f%%\n', avr_percentage);
else
    fprintf('AVR (g >= 1)               | Data Not Found\n');
end
fprintf('---------------------------|----------------------\n');
for i = 1:length(lambda_values)
    lambda = lambda_values(i);
    if ~isnan(mdr_percentages(i))
        fprintf('MDR (Z >= 0, lambda=%.3f)    | %.3f%%\n', lambda, mdr_percentages(i));
    else
        fprintf('MDR (Z >= 0, lambda=%.3f)    | Data Not Found\n', lambda);
    end
end
fprintf('======================================================\n');