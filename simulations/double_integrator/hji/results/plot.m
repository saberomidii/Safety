% This script visualizes the results of the MDR value iteration for the
% double-integrator system. It loads the saved value function and creates
% several plots to analyze the results.
clear; clc;

fprintf('Plotting results from Julia simulation...\n');

% ---------------------------
% Parameters (must match the Julia script)
% ---------------------------
X1_MIN = -1.0; X1_MAX = 5.0;
X2_MIN = -5.0; X2_MAX = 5.0;
NUM_POINTS_STATE_1 = 161;
NUM_POINTS_STATE_2 = 161;

% Safe Box K
K1_MIN = 0.0; K1_MAX = 4.0;
K2_MIN = -3.0; K2_MAX = 3.0;

% ---------------------------
% Load Data and Reconstruct Grid
% ---------------------------
try
    Z = readmatrix('results/Value_function.txt')';
    fprintf('Successfully loaded value_function.csv.\n');
catch ME
    error('Could not find or read "results/value_function_and_levels_sto.csv". Please run the Julia script first to generate the file.');
end

x1_grid = linspace(X1_MIN, X1_MAX, NUM_POINTS_STATE_1);
x2_grid = linspace(X2_MIN, X2_MAX, NUM_POINTS_STATE_2);

if ~isequal(size(Z), [NUM_POINTS_STATE_1, NUM_POINTS_STATE_2])
    error('The dimensions of the loaded value function do not match the grid parameters.');
end

% Z is upper bound 
    L = sqrt(5);
    TAU_BAR = 2;
    DISCOUNT_RATE= 0.1;
    coef_lower_bound = exp(-DISCOUNT_RATE*TAU_BAR);
    Lower_bound = (Z - L * (1 - coef_lower_bound)) / coef_lower_bound;
% ---------------------------
% Plot 1: Filled Contour with Safe Box
% ---------------------------
box_width = K1_MAX - K1_MIN;
box_height = K2_MAX - K2_MIN;
figure('Name', 'Value Function Contour');
contourf(x1_grid, x2_grid, Z, 30); % Transpose Z to match axes
colorbar;
xlabel('x_1 (position)');
ylabel('x_2 (velocity)');
title('Value Function Contour');
hold on;
rectangle('Position',[K1_MIN, K2_MIN, box_width, box_height], 'EdgeColor','k','LineWidth',2);
hold off;

% ---------------------------
% Plot 2: Zero Level Set with Safe Box
% ---------------------------

figure('Name', 'Zero Level Set');
rectangle('Position',[K1_MIN, K2_MIN, box_width, box_height], 'EdgeColor','k','LineWidth',2);
xlabel('x_1 (position)');
ylabel('x_2 (velocity)');
title('Zero Level Set of Value Function');
hold on;
contour(x1_grid, x2_grid, Z, [0.0, 0.0], 'green', 'LineWidth', 2);  % Plot only the zero level set
contour(x1_grid, x2_grid, Lower_bound, [0.0, 0.0], 'green','LineStyle','--', 'LineWidth', 2);  % Plot only the zero level set
legend('Upper bound', 'Lower bound');
hold off;

% ---------------------------
% Plot 3: 3D Surface
% ---------------------------
figure('Name', 'Value Function Surface');
surf(x1_grid, x2_grid, Z); % Transpose Z to match axes
shading interp;
colorbar;
xlabel('x_1 (position)');
ylabel('x_2 (velocity)');
zlabel('Value Z(x)');
title('Value Function Surface');

% ----------------------------------------------------
% Count Nodes Inside Each Safe Set Approximation
% ----------------------------------------------------
fprintf('\n--- Node Count & Percentage Results ---\n');

% Total number of nodes in the state space grid
total_nodes = NUM_POINTS_STATE_1 * NUM_POINTS_STATE_2;

% --- Count for the Upper Bound (Z >= 0) ---
% Create a logical matrix: 1 where Z>=0, 0 otherwise
upper_bound_mask = (Z >= 0);
% Count the number of 'true' (1) elements
nodes_in_upper_bound = sum(upper_bound_mask, 'all');
percent_upper = (nodes_in_upper_bound / total_nodes) * 100;

fprintf('Upper Bound (Z>=0): %d / %d nodes (%.2f%%)\n', ...
        nodes_in_upper_bound, total_nodes, percent_upper);


% --- Count for the Lower Bound (V >= 0) ---
% Create a logical matrix for the lower bound
lower_bound_mask = (Lower_bound >= 0);
% Count the number of 'true' (1) elements
nodes_in_lower_bound = sum(lower_bound_mask, 'all');
percent_lower = (nodes_in_lower_bound / total_nodes) * 100;

fprintf('Lower Bound (V>=0): %d / %d nodes (%.2f%%)\n', ...
        nodes_in_lower_bound, total_nodes, percent_lower);
