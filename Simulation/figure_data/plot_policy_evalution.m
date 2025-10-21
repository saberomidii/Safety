% --- MATLAB script to load data and generate plots with new style ---
clear; clc; close all;

% --- Load Data ---
fprintf("Loading data from 'figure_data' directory...\n");



trajectory_1   = readmatrix(fullfile('trajectory_1.csv'));
trajectory_2   = readmatrix(fullfile('trajectory_2.csv'));
trajectory_3   = readmatrix(fullfile('trajectory_3.csv'));
trajectory_4   = readmatrix(fullfile('trajectory_4.csv'));
outer_track    = readmatrix(fullfile('outer_track.csv'));
inner_track    = readmatrix(fullfile('inner_track.csv'));

% Extract data into variables for clarity
X_1 = trajectory_1(:, 1);
Y_1 = trajectory_1(:, 2);
X_2 = trajectory_2(:, 1);
Y_2 = trajectory_2(:, 2);
X_3 = trajectory_3(:, 1);
Y_3 = trajectory_3(:, 2);
X_4 = trajectory_4(:, 1);
Y_4 = trajectory_4(:, 2);



% --- Plot: Trajectories with new style ---
fprintf("Generating trajectory plot with updated style...\n");
fig1 = figure('Position', [100, 100, 800, 600], 'Color', 'w'); % Set figure background to white
hold on;

% Plot the track boundaries (constraint set)
% We only need one legend entry for the track, so we use 'HandleVisibility','off' for the second part
p_track = plot(outer_track(:,1), outer_track(:,2), 'k-', 'LineWidth', 2);
if ~isempty(inner_track)
    plot(inner_track(:,1), inner_track(:,2), 'k-', 'LineWidth', 2, 'HandleVisibility','off');
end

% Plot the trajectories
p1 = plot(X_1, Y_1, '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5); % Blue
p2 = plot(X_2, Y_2, '-', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 2.5); % Green
p3 = plot(X_3, Y_3, '-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 2.5); % Red/Maroon
p4 = plot(X_4, Y_4, '-', 'Color', [0.6350, 0.4660, 0.1840], 'LineWidth', 2.5); % Red/Maroon

% Plot start points (black triangles)
% Only create a legend entry for the first one to avoid clutter
p_start = scatter(X_1(1), Y_1(1), 150, 'k^', 'filled', 'LineWidth', 1.5);
scatter(X_2(1), Y_2(1), 150, 'k^', 'filled', 'LineWidth', 1.5, 'HandleVisibility', 'off');
scatter(X_3(1), Y_3(1), 150, 'k^', 'filled', 'LineWidth', 1.5, 'HandleVisibility', 'off');
scatter(X_4(1), Y_4(1), 150, 'k^', 'filled', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Plot end points (red stars)
% Only create a legend entry for the first one
p_end = scatter(X_1(end), Y_1(end), 250, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5); % Orange pentagram
scatter(X_2(end), Y_2(end), 250, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5, 'HandleVisibility', 'off');
scatter(X_3(end), Y_3(end), 250, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5, 'HandleVisibility', 'off');
scatter(X_4(end), Y_4(end), 250, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5, 'HandleVisibility', 'off');

% --- Setup plot with LaTeX interpreter and updated terminology ---
xlim([-3, 3]);
ylim([-1.75, 1.75]);
grid on;
grid minor;
box on;
axis equal; % Ensure aspect ratio is correct

xlabel('$\mathrm{X}$ (meters)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Y}$ (meters)', 'Interpreter', 'latex', 'FontSize', 18);

% Create a more descriptive legend using the plot handles
legend([p_track, p1, p2, p3,p4, p_start, p_end], ...
       {'Constraint Set', 'Trajectory 1', 'Trajectory 2', 'Trajectory 3','Trajectory 4', 'Start Points', 'End Points'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12);

hold off;

% --- Save the figure as a new PDF ---
fprintf("Saving updated trajectory plot as PDF...\n");
exportgraphics(fig1, 'AVR_policy_trajectory.pdf', 'ContentType', 'vector');

fprintf("Script finished successfully.\n");