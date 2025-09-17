clear
clc
close all

% Specify the filenames of the CSV files
file_sqrt_5 = 'approximation_results_sqrt(5).csv';
file_2_sqrt_5 = 'approximation_results_2sqrt(5).csv';
file_3_sqrt_5 = 'approximation_results.csv';

% Read the data from the CSV files into matrices
% This assumes the files are in the current MATLAB directory.
% If not, provide the full path to the files.
try
    sqrt_5 = readmatrix(file_sqrt_5);
    sqrt_2_5 = readmatrix(file_2_sqrt_5);
    sqrt_3_5 = readmatrix(file_3_sqrt_5);
catch ME
    disp('Error reading files. Make sure the CSV files are in the correct path.');
    disp(ME.message);
    return;
end

% --- Scale the approximation data to percentages ---
% Multiply the over-approximation (col 2) and under-approximation (col 3) by 100
sqrt_5(:,2:3) = sqrt_5(:,2:3);
sqrt_2_5(:,2:3) = sqrt_2_5(:,2:3);
sqrt_3_5(:,2:3) = sqrt_3_5(:,2:3);

% Create a new figure window for the combined plot
figure;

% --- Plot 1: Results for sqrt(5) ---
plot(sqrt_5(:,1), sqrt_5(:,2), '-', 'LineWidth', 2, 'DisplayName', 'Over-approx for $\sqrt{5}$');
hold on; % Keep the plot active to add more lines
plot(sqrt_5(:,1), sqrt_5(:,3), '--', 'LineWidth', 2, 'DisplayName', 'Under-approx for $\sqrt{5}$');

% --- Plot 2: Results for 2*sqrt(5) ---
plot(sqrt_2_5(:,1), sqrt_2_5(:,2), '-', 'LineWidth', 2, 'DisplayName', 'Over-approx for $2\sqrt{5}$');
plot(sqrt_2_5(:,1), sqrt_2_5(:,3), '--', 'LineWidth', 2, 'DisplayName', 'Under-approx for $2\sqrt{5}$');

% --- Plot 3: Results for 3*sqrt(5) ---
plot(sqrt_3_5(:,1), sqrt_3_5(:,2), '-', 'LineWidth', 2, 'DisplayName', 'Over-approx for $3\sqrt{5}$');
plot(sqrt_3_5(:,1), sqrt_3_5(:,3), '--', 'LineWidth', 2, 'DisplayName', 'Under-approx for $3\sqrt{5}$');

% --- Add the constant line for Average Reward ---
% The value is 29, which is now consistent with the scaled percentage data.
yline(29, 'k:', 'LineWidth', 4, 'DisplayName', 'Average Reward');

% Release the plot hold
hold off;

% --- Add Plot Enhancements ---
grid on;
ax = gca;
ax.FontSize = 12;

% Add title and labels using the LaTeX interpreter
title('Effect of L and $\lambda$ on the Discounted Minimum Reward Approach for Double Integrator', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);

% --- MODIFIED LINE ---
% Updated y-axis label to include the percentage symbol
ylabel('Approximated Proportion of Safe Set (\%)', 'Interpreter', 'latex', 'FontSize', 14);

% Create a legend to identify all lines
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);