clear all
clc
close all

% Set LaTeX interpreters for professional-looking plots
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');

%% --- Configuration (Should match your Julia script) ---
% Path to the results folder from the Julia simulation
results_folder = 'results_dubins_MDR_sigma_0.0';

% Grid parameters used in the Julia script
num_points_state_1 = 201; % x-dimension
num_points_state_2 = 201; % y-dimension
num_points_state_3 =61; % theta-dimension

% Domain limits
x1_min = -2.0; x1_max = 2.0;
x2_min = -2.0; x2_max = 2.0;

%[cite_start]% Parameters for calculating over-approximation bound [cite: 233]
discount_rate = 0.1;
tau_bar = 2.0;

%% --- Load and Reshape Data ---
fprintf('Loading data from %s...\n', results_folder);

% Load the flattened Z value function and the 3D state grid
Z_val_flat = readmatrix(fullfile(results_folder, 'Z_value_function.csv'));
state_3d   = readmatrix(fullfile(results_folder, 'state_3d.csv'));

% Reshape the flat Z vector into a 3D grid corresponding to (x, y, theta)
Z_grid_3D = reshape(Z_val_flat, num_points_state_1, num_points_state_2, num_points_state_3);

% For visualization in 2D, we project the 3D value function.
% We take the minimum value across the theta dimension. This represents the
% worst-case safety value for any orientation at a given (x,y) position.
Z_grid_2D = min(Z_grid_3D, [], 3);

% Create the x and y grid vectors for plotting
x1_grid = linspace(x1_min, x1_max, num_points_state_1);
x2_grid = linspace(x2_min, x2_max, num_points_state_2);

%% --- Calculate Bounds and Thresholds ---
% Estimate L from the maximum value of Z, as Z = U + L and max(U) is approx. [cite_start]0 [cite: 201]
L = max(Z_val_flat);

%[cite_start]% Threshold for the under-approximation (viable set) [cite: 224]
under_approx_threshold = 0.0;

%[cite_start]% Threshold for the over-approximation (reachable set) [cite: 233]
over_approx_threshold = L * (1 - exp(-discount_rate * tau_bar));

%% --- Plotting ---
figure;
set(gcf, 'Position', [100, 100, 600, 550], 'Color', 'w');
hold on;

% -- 1. Plot the Under-Approximation (Viable Set) --
% This is the region where Z(x) <= 0. We fill this area.
[C_under, ~] = contour(x1_grid, x2_grid, Z_grid_2D', [under_approx_threshold, under_approx_threshold]);
fill_contour(C_under, [0.3 0.7 0.3]); % Fill with a light green color

% -- 2. Plot the Over-Approximation Contour --
% This is the line where Z(x) = L*(1 - exp(-λ*τ̄))
contour(x1_grid, x2_grid, Z_grid_2D', [over_approx_threshold, over_approx_threshold], ...
    'Color', [0.8 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');

% -- 3. Draw the Original Target Set (Elliptical Ring) --
t = linspace(0, 2*pi, 200);
% Ellipse parameters from the Julia script's s_d function
a_outer = 0.4; b_outer = 0.9;
a_inner = 0.2; b_inner = 0.7;
plot(a_outer*cos(t), b_outer*sin(t), 'k-', 'LineWidth', 2);
plot(a_inner*cos(t), b_inner*sin(t), 'k-', 'LineWidth', 2);

% -- 4. Configure Plot Appearance --
axis equal;
grid on;
box on;

xlabel('$x_1$ (Position)', 'FontSize', 16);
ylabel('$x_2$ (Position)', 'FontSize', 16);
title('MDR Viable Set for Dubins Car', 'FontSize', 18);
xlim([x1_min, x1_max]);
ylim([x2_min, x2_max]);
set(gca, 'FontSize', 12);

% -- 5. Create a Legend --
h_under = patch(NaN, NaN, [0.3 0.7 0.3], 'EdgeColor', 'k');
h_over  = plot(NaN, NaN, 'Color', [0.8 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
h_target = plot(NaN, NaN, 'k-', 'LineWidth', 2);
legend([h_under, h_over, h_target], ...
    {'Viable Set (Under-Approx.)', 'Reachable Set (Over-Approx.)', 'Target Set'}, ...
    'Location', 'best', 'FontSize', 12);

hold off;

%% --- Helper Function to Fill Contours ---
function fill_contour(C, color)
    % Helper function to parse contour matrix C and fill the areas
    idx = 1;
    while idx < size(C, 2)
        level = C(1, idx);
        n_points = C(2, idx);
        x_coords = C(1, idx+1 : idx+n_points);
        y_coords = C(2, idx+1 : idx+n_points);
        fill(x_coords, y_coords, color, 'EdgeColor', 'k', 'LineWidth', 2);
        idx = idx + n_points + 1;
    end
end