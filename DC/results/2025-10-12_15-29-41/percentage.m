% =========================================================================
% MATLAB Script to Calculate and Compare Safe Set Percentages
% =========================================================================
% Description:
% This script loads the pre-computed AVR and MDR data from CSV files
% and calculates the percentage of the state space that is considered
% safe according to each method and discount factor.
%
% Author: Gemini
% Date: October 12, 2025
% =========================================================================

clear all;
clc;
close all;

%% --- 1. Load and Prepare Data ---
results_dir = './'; % Assumes script is in the same folder as the CSVs

% --- AVR Calculation ---
try
    G_flat = readmatrix(fullfile(results_dir, 'AVR_gain_map.csv'));
    % For AVR, the safe set is where the gain g(s) is approximately 1
    count_G = sum(G_flat >= 0.999); % Use a tolerance for floating point
    G_per = (count_G / numel(G_flat)) * 100;
catch ME
    error('Failed to load AVR_gain_map.csv. Error: %s', ME.message);
end

% --- MDR Calculation (Undiscounted and Discounted) ---
lambda_values = [0.0, 0.01, 0.1, 0.2];
Z_per = NaN(size(lambda_values)); % Pre-allocate array for percentages

for i = 1:length(lambda_values)
    lambda = lambda_values(i);
    mdr_filename = sprintf('MDR_Z_map_lambda_%.1f.csv', lambda);
    
    try
        Z_flat = readmatrix(fullfile(results_dir, mdr_filename));
        % For MDR, the safe set is where the value Z(x) is non-negative
        count_Z = sum(Z_flat >= 0);
        Z_per(i) = (count_Z / numel(Z_flat)) * 100;
        
    catch ME
        warning('Could not load %s. Skipping calculation for this lambda. Error: %s', mdr_filename, ME.message);
    end
end

%% --- 2. Display Results in a Formatted Table ---

fprintf('\n======================================================\n');
fprintf('      Safe Set Percentage of Total State Space      \n');
fprintf('======================================================\n');
fprintf('Method                     | Safe Set Percentage\n');
fprintf('---------------------------|----------------------\n');
fprintf('AVR (g >= 1)               | %.3f%%\n', G_per);
fprintf('---------------------------|----------------------\n');

for i = 1:length(lambda_values)
    lambda = lambda_values(i);
    if ~isnan(Z_per(i))
        fprintf('MDR (Z >= 0, lambda=%.2f)    | %.3f%%\n', lambda, Z_per(i));
    else
        fprintf('MDR (Z >= 0, lambda=%.2f)    | Data Not Found\n', lambda);
    end
end
fprintf('======================================================\n');