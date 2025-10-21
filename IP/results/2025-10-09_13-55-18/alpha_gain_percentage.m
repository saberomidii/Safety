% =========================================================================
% MATLAB Script to Calculate AVR Level Set Percentages
% =========================================================================
clear;
clc;
close all;

%% --- 1. Define Files and Parameters ---
results_dir = './'; % Assumes script is in the same directory as the CSVs

% --- AVR Data ---
avr_file = 'AVR_gain_map.csv';

% --- Define the gain level set percentages to analyze ---
% These represent the threshold P for calculating P% of the data where g >= P/100
level_set_percentages = [20, 30, 40, 50, 60, 70, 80, 90]; 

%% --- 2. Perform Calculations ---
% Initialize a structure to store results
avr_level_set_results = struct('threshold', num2cell(level_set_percentages), 'percentage', NaN);
data_vector = [];

% --- AVR Data Loading ---
try
    % readtable is more robust to formatting errors like trailing commas
    avr_table = readtable(fullfile(results_dir, avr_file), 'ReadVariableNames', false);
    data_vector = table2array(avr_table); % Convert table to a matrix
    data_vector = data_vector(:); % Flatten into a single column vector
    total_count = numel(data_vector);
    
catch ME
    warning('Could not process %s. Setting all level set percentages to NaN.', avr_file);
    total_count = 0; % Ensure division by zero doesn't happen if no data
end

% --- AVR Level Set Calculation ---
if total_count > 0 % Only proceed if data was successfully loaded
    for i = 1:length(level_set_percentages)
        P = level_set_percentages(i);
        threshold_g = P / 100; % e.g., 20% corresponds to g >= 0.2
        
        % Count points where the gain g is greater than or equal to the threshold
        level_set_count = sum(data_vector >= threshold_g);
        
        % Calculate the percentage
        avr_level_set_results(i).percentage = (level_set_count / total_count) * 100;
    end
end

%% --- 3. Display Results ---
fprintf('\n======================================================\n');
fprintf('     AVR Level Set Percentage of Total State Space    \n');
fprintf('======================================================\n');
fprintf('Gain Threshold (g >= P/100) | State Space Percentage\n');
fprintf('----------------------------|------------------------\n');

if isempty(data_vector)
    fprintf('Data File Not Found or Empty\n');
else
    for i = 1:length(avr_level_set_results)
        P = avr_level_set_results(i).threshold;
        percentage = avr_level_set_results(i).percentage;
        
        if ~isnan(percentage)
            % Format the output string: P/100 corresponds to P%
            fprintf('AVR (g >= %.2f)             | %.3f%%\n', P/100, percentage); 
        else
            fprintf('AVR (g >= %.2f)             | Calculation Error\n', P/100);
        end
    end
end
fprintf('======================================================\n');