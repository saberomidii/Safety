% This script visualizes results in two separate figures, one for each
% disturbance case (sigma=0 and sigma=1), and saves both as PDF files.
clear; clc;
close all;
fprintf('Plotting results from Julia simulation...\n');

% ---------------------------
% Load Data
% ---------------------------
try
    Z_lambda_0_sigma_0 = readmatrix('Value_function_lambda_0_sigma_0.txt');
    Z_lambda_2_sigma_0 = readmatrix('Value_function_lambda_2_sigma_0.txt');
    Z_lambda_2_sigma_0_L_sqrt5 = readmatrix('Value_function_lambda_2_sigma_0_lsqrt(5).txt');
    Z_lambda_2_sigma_1_L_sqrt5 = readmatrix('Value_function_lambda_2_sigma_1_lsqrt(5).txt');
    Z_lambda_0_sigma_1 = readmatrix('Value_function_lambda_0_sigma_1.txt');
    Z_lambda_2_sigma_1 = readmatrix('Value_function_lambda_2_sigma_1.txt');
    Z_interpolation_path_2 = "/Users/saber/Desktop/Safety/simulations/double_integrator/hji/interpolation_results/Value_function_lambda_0.2.txt";
    Z_Int_2 = readmatrix(Z_interpolation_path_2);
    Z_interpolation_path_0 = "/Users/saber/Desktop/Safety/simulations/double_integrator/hji/interpolation_results/Value_function_lambda_0.0.txt";
    Z_Int_0 = readmatrix(Z_interpolation_path_0);
    Z_interpolation_path_sigma_1 = '/Users/saber/Desktop/Safety/simulations/double_integrator/hji/interpolation_results_2player/Value_function_lambda_0.0_2player.txt';
    Z_Int_sigma_1= readmatrix(Z_interpolation_path_sigma_1);
    Path_G_sigma_0 = "/Users/saber/Desktop/Safety/simulations/double_integrator/average_reward/results_primal/0.0/G.csv";
    Path_G_sigma_1 = "/Users/saber/Desktop/Safety/simulations/double_integrator/average_reward/results_primal/1.0/G.csv";
    G_average_sigma_0 = readmatrix(Path_G_sigma_0)';
    G_average_sigma_1 = readmatrix(Path_G_sigma_1)';
    fprintf('Successfully loaded data files.\n');
catch ME
    error('Could not find or read a required data file. Please check file paths. Error: %s', ME.message);
end

% ---------------------------
% Reconstruct Grid
% ---------------------------
NUM_POINTS_STATE_1 = 161;
NUM_POINTS_STATE_2 = size(G_average_sigma_0, 2);
if size(Z_lambda_0_sigma_0, 2) ~= NUM_POINTS_STATE_2
    NUM_POINTS_STATE_2 = size(Z_lambda_0_sigma_0, 2);
end
x1 = linspace(-1, 5, NUM_POINTS_STATE_1);
x2 = linspace(-5, 5, NUM_POINTS_STATE_2);

% ###########################################################################
% ## FIGURE 1: RESULTS FOR SIGMA = 0
% ###########################################################################

fig1 = figure('Position', [100, 100, 600, 500]);
hold on;
plotBoxAndCorners();

% Plot all sigma=0 contours
contour(x1, x2, Z_lambda_0_sigma_0', [0.0, 0.0], "LineWidth", 2.5, "Color", "#0000CD", 'LineStyle', '-', 'DisplayName', 'MDR: $\lambda=0.0, L=2\sqrt{5}$');
contour(x1, x2, Z_lambda_2_sigma_0', [0.0, 0.0], "LineWidth", 2.5, "Color", "#1E90FF", 'LineStyle', '--', 'DisplayName', 'MDR: $\lambda=0.2, L=2\sqrt{5}$');
contour(x1, x2, Z_lambda_2_sigma_0_L_sqrt5', [0.0, 0.0], "LineWidth", 2.5, "Color", "#483D8B", 'LineStyle', ':', 'DisplayName', 'MDR: $\lambda=0.2, L=\sqrt{5}$');
contour(x1, x2, Z_Int_0', [0.0, 0.0], "LineWidth", 3, "Color", "#FF0000", 'LineStyle', '-', 'DisplayName', 'MDR (Interp): $\lambda=0.0, L=\sqrt{5}$');
contour(x1, x2, G_average_sigma_0', [1.0, 1.0], "LineWidth", 2.5, "Color", "#006400", 'LineStyle', '-', 'DisplayName', 'Average Reward');

hold off;
setupAxes('', true); % Show legend

% --- Save the first figure ---
fprintf('\nSaving figure 1 to PDF...\n');
exportgraphics(fig1, 'comparison_sigma0.pdf', 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Figure saved as comparison_sigma0.pdf\n');


% ###########################################################################
% ## FIGURE 2: RESULTS FOR SIGMA = 1
% ###########################################################################

fig2 = figure('Position', [750, 100, 600, 500]);
hold on;
plotBoxAndCorners();

% Plot all sigma=1 contours
contour(x1, x2, Z_lambda_0_sigma_1', [0.0, 0.0], "LineWidth", 2.5, "Color", "#FF4500", 'LineStyle', '-', 'DisplayName', 'MDR: $\lambda=0.0, L=2\sqrt{5}$');
contour(x1, x2, Z_lambda_2_sigma_1', [0.0, 0.0], "LineWidth", 2.5, "Color", "#FFA500", 'LineStyle', '--', 'DisplayName', 'MDR: $\lambda=0.2, L=2\sqrt{5}$');
contour(x1, x2, Z_lambda_2_sigma_1_L_sqrt5', [0.0, 0.0], "LineWidth", 2.5, "Color", "#9400D3", 'LineStyle', ':', 'DisplayName', 'MDR: $\lambda=0.2, L=\sqrt{5}$');
contour(x1, x2, Z_Int_sigma_1', [0.0, 0.0], "LineWidth", 3, "Color", "blue", 'LineStyle', '-', 'DisplayName', 'MDR (Interp): $\lambda=0.0, L=\sqrt{5}$');
contour(x1, x2, G_average_sigma_1', [1.0, 1.0], "LineWidth", 2.5, "Color", "#228B22", 'LineStyle', '-', 'DisplayName', 'Average Reward');

hold off;
setupAxes('', true); % Show legend

% --- Save the second figure ---
fprintf('\nSaving figure 2 to PDF...\n');
exportgraphics(fig2, 'comparison_sigma1.pdf', 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Figure saved as comparison_sigma1.pdf\n');


% ---------------------------
% Helper Functions
% ---------------------------
function plotBoxAndCorners()
    K1_MIN = 0.0; K1_MAX = 4.0; K2_MIN = -3.0; K2_MAX = 3.0;
    box_x = [K1_MIN, K1_MAX, K1_MAX, K1_MIN, K1_MIN];
    box_y = [K2_MIN, K2_MIN, K2_MAX, K2_MAX, K2_MIN];
    plot(box_x, box_y, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Safe Set Boundary');
    corners_x = [4, 0, 4, 0]; corners_y = [-3, -3, 3, 3];
    plot(corners_x, corners_y, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'HandleVisibility', 'off');
end

function setupAxes(title_text, show_legend)
    grid on; 
    axis equal;
    
    title(title_text, 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('State $x_1$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('State $x_2$', 'Interpreter', 'latex', 'FontSize', 14);
    
    if show_legend
        lgd = legend('show');
        set(lgd, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 8);
    end
end