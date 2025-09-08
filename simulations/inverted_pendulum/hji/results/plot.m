% This script visualizes the results of the MDR value iteration for the
% double-integrator system. It loads the saved value function and creates
% several plots to analyze the results.

clear; clc;

fprintf('Plotting results from Julia simulation...\n');

% ---------------------------
% Parameters (must match the Julia script)
% ---------------------------
X1_MIN = -0.5; X1_MAX = 0.5;
X2_MIN = -1.0; X2_MAX = 1.0;
NUM_POINTS_STATE_1 = 201;
NUM_POINTS_STATE_2 = 201;

% Safe Box K
K1_MIN = -0.3; K1_MAX = 0.3;
K2_MIN = -0.6; K2_MAX = 0.6;

% ---------------------------
% Load Data and Reconstruct Grid
% ---------------------------
try
    Z = readmatrix('results/value_function_and_levels_deter.txt')';
    fprintf('Successfully loaded value_function.csv.\n');
catch ME
    error('Could not find or read "results/value_function_and_levels_sto.csv". Please run the Julia script first to generate the file.');
end

x1_grid = linspace(X1_MIN, X1_MAX, NUM_POINTS_STATE_1);
x2_grid = linspace(X2_MIN, X2_MAX, NUM_POINTS_STATE_2);

if ~isequal(size(Z), [NUM_POINTS_STATE_1, NUM_POINTS_STATE_2])
    error('The dimensions of the loaded value function do not match the grid parameters.');
end

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
DISCOUNT_RATE =0.1;
TAU_BAR = 2;
L = 0.44;
under_approx_level = L * (1 - exp(-DISCOUNT_RATE * TAU_BAR));


figure('Name', 'Zero Level Set');
rectangle('Position',[K1_MIN, K2_MIN, box_width, box_height], 'EdgeColor','k','LineWidth',2);
xlabel('x_1 (position)');
ylabel('x_2 (velocity)');
title('Zero Level Set of Value Function');
hold on;
contour(x1_grid, x2_grid, Z, [0.0, 0.0], 'green', 'LineWidth', 2);  % Plot only the zero level set
contour(x1_grid, x2_grid, Z, [under_approx_level, under_approx_level], 'green','LineStyle','--', 'LineWidth', 2);  % Plot only the zero level set
legend('Over-Approximation', 'Under-Approximation');
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

% ---------------------------
% Extract Zero-Level Set Contour and Compute Area
% ---------------------------
C = contourc(x1_grid, x2_grid, Z, [0.0 0.0]); % Get contour matrix for zero level set

total_area = 0;
k = 1;
while k < size(C,2)
    n = C(2,k);               % Number of points in this contour segment
    x_contour = C(1, k+1:k+n);
    y_contour = C(2, k+1:k+n);
    
    % Area of polygon defined by this contour segment
    area_segment = polyarea(x_contour, y_contour);
    total_area = total_area + area_segment;
    
    k = k + n + 1;            % Move to the start of the next contour segment
end

fprintf('\nTotal area enclosed by the zero-level set: %f\n', total_area);

% Total area of the state space grid
state_space_area = (X1_MAX - X1_MIN) * (X2_MAX - X2_MIN);
safe_set_percentage = (total_area / state_space_area) * 100;

fprintf('The safe set represents %.2f%% of the state space.\n', safe_set_percentage);