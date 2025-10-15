% --- MATLAB script to load data and generate plots with new style ---
clear; clc; close all;

% --- Load Data ---
fprintf("Loading data from 'figure_data' directory...\n");
data_dir = 'figure_data';

trajectory   = readmatrix(fullfile(data_dir, 'trajectory.csv'));
control_input= readmatrix(fullfile(data_dir, 'control_input.csv'));
outer_track  = readmatrix(fullfile(data_dir, 'outer_track.csv'));
inner_track  = readmatrix(fullfile(data_dir, 'inner_track.csv'));
start_point  = readmatrix(fullfile(data_dir, 'start_point.csv'));
longest_time = readmatrix(fullfile(data_dir, 'longest_time.txt'));

% Extract data into variables for clarity
best_X_traj = trajectory(:, 1);
best_Y_traj = trajectory(:, 2);
best_x = start_point(1);
best_y = start_point(2);

% --- Plot 1: Trajectory with new style ---
fprintf("Generating trajectory plot with updated style...\n");
fig1 = figure('Position', [100, 100, 800, 600], 'Color', 'w'); % Set figure background to white
hold on;

% Plot the constraint set boundary as a solid black line
plot(outer_track(:,1), outer_track(:,2), 'k-', 'LineWidth', 2);
if ~isempty(inner_track)
    plot(inner_track(:,1), inner_track(:,2), 'k-', 'LineWidth', 2);
end

% Plot the trajectory (blue), start point (green), and end point (black star)
plot(best_X_traj, best_Y_traj, 'b-', 'LineWidth', 2.5);
scatter(best_x, best_y, 200, 'g', 'filled');
scatter(best_X_traj(end), best_Y_traj(end), 200, 'red', '^','filled');

% Setup plot with LaTeX interpreter and updated terminology
% title(sprintf('AVR Policy (Time: %d steps)', longest_time), 'Interpreter', 'latex', 'FontSize', 20);
xlabel(' $\mathrm{X}$(meter)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel(' $\mathrm{Y}$(meter)', 'Interpreter', 'latex', 'FontSize', 18);
legend({'Constraint Set', 'AVR Policy', 'Start Point', 'End Point'}, 'Interpreter', 'latex', 'Location', 'northeast');
axis equal;
box on;
grid on;
hold off;

% Save the figure as a new PDF
fprintf("Saving updated trajectory plot as PDF...\n");
exportgraphics(fig1, 'AVR_policy_trajectory.pdf', 'ContentType', 'vector');
