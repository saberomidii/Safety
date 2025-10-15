%% --- Main Script ---
clear; clc; close all;

% --- Setup: Ensure a dummy policy file exists ---
NUM_STATES_TOTAL = 101 * 61 * 15;
NUM_ACTIONS = 21;
policy_filename = "optimal_policy_avr.csv";
create_dummy_policy_file(policy_filename, NUM_STATES_TOTAL, NUM_ACTIONS);
% -----------------------------------------------------------------------

fprintf("Loading policy and setting up state space...\n");
policy_average = readmatrix(policy_filename);

% --- ENVIRONMENT AND DYNAMICS DEFINITIONS ---
L_straight = 4.0; % Length of the straight sections
R_outer = 3.0;    % Radius of the outer semicircles
R_inner = 0.0;    % Radius of the inner semicircles
xc_left = -L_straight / 2.0;
xc_right = L_straight / 2.0;
y_c = 0.0;

% --- STATE SPACE DISCRETIZATION ---
fprintf("Discretizing state space and building KD-Tree...\n");
x1_range = linspace(-5.0, 5.0, 101);
x2_range = linspace(-3.0, 3.0, 61);
x3_range = linspace(0, 2*pi, 15);

[X1, X2, X3] = meshgrid(x1_range, x2_range, x3_range);
states_3d = [X1(:), X2(:), X3(:)]; % N x 3 matrix of states
tree = KDTreeSearcher(states_3d);

% --- SIMULATION PARAMETERS ---
u_actions = linspace(-2.0, 2.0, NUM_ACTIONS);
max_sim_time = 2000;

% --- MAIN EXECUTION: SEARCH FOR THE BEST STARTING POINT ---
longest_time = 0;
best_starting_point = [0.0, 0.0, 0.0];
total_states_to_check = size(states_3d, 1);

fprintf("\n--- Searching for the best starting point out of %d candidates ---\n", total_states_to_check);

% Iterate through every possible starting state from our grid
for i = 1:total_states_to_check
    start_state = states_3d(i, :);
    sx = start_state(1);
    sy = start_state(2);
    st = start_state(3);

    % Provide progress updates
    if mod(i, 20000) == 0 % Adjusted for faster MATLAB loops
        progress = 100 * i / total_states_to_check;
        fprintf("... Progress: %.1f%% (%d / %d states checked)\n", progress, i, total_states_to_check);
    end

    % Optimization: skip starting points that are already outside the track
    if ~is_point_safe(sx, sy, L_straight, R_outer, R_inner)
        continue;
    end

    % Run simulation (only need trajectory length for search)
    [X_traj, ~, ~] = run_simulation(sx, sy, st, tree, policy_average, u_actions, max_sim_time, L_straight, R_outer, R_inner);
    current_time_safe = length(X_traj);

    % If this run is the best so far, save it
    if current_time_safe > longest_time
        longest_time = current_time_safe;
        best_starting_point = start_state;
        fprintf("✅ New best found! Time: %d steps, Start: (x=%.2f, y=%.2f, theta=%.2f)\n", longest_time, sx, sy, st);
    end
end

fprintf("\n--- 🏆 Search Complete! ---\n");
fprintf("Best Starting Point: (%.2f, %.2f, %.2f)\n", best_starting_point(1), best_starting_point(2), best_starting_point(3));
fprintf("Longest Time Safe: %d steps\n", longest_time);

% --- VISUALIZATION OF THE BEST TRAJECTORY ---
fprintf("Visualizing the best trajectory...\n");

% Run the simulation one last time with the best starting point
best_x = best_starting_point(1);
best_y = best_starting_point(2);
best_theta = best_starting_point(3);
[best_X_traj, best_Y_traj, best_U_traj] = run_simulation(best_x, best_y, best_theta, tree, policy_average, u_actions, max_sim_time, L_straight, R_outer, R_inner);

% --- Plot 1: Trajectory ---
fig1 = figure('Position', [100, 100, 800, 600]);
hold on;

% Plot the racetrack shape
[outer_track_x, outer_track_y] = generate_racetrack_points(R_outer, L_straight);
[inner_track_x, inner_track_y] = generate_racetrack_points(R_inner, L_straight);

fill(outer_track_x, outer_track_y, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
if R_inner > 0
    fill(inner_track_x, inner_track_y, [1 1 1], 'EdgeColor', 'r');
end

% Plot the trajectory and start/end points
plot(best_X_traj, best_Y_traj, 'm', 'LineWidth', 2.5);
scatter(best_x, best_y, 50, 'g', 'filled');
scatter(best_X_traj(end), best_Y_traj(end), 50, 'r', 'filled');

% Setup plot with LaTeX interpreter
title(sprintf('Longest Safe Trajectory (Time: %d steps)', longest_time), 'Interpreter', 'latex', 'FontSize', 14);
xlabel('x position ($x$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('y position ($y$)', 'Interpreter', 'latex', 'FontSize', 12);
legend({'Safe Region', 'Longest Trajectory', 'Start Point', 'End Point'}, 'Interpreter', 'latex', 'Location', 'northeast');
axis equal;
box on;
grid on;
hold off;

% Save the figure as a PDF
fprintf("Saving trajectory plot as PDF...\n");
exportgraphics(fig1, 'trajectory_plot.pdf', 'ContentType', 'vector');

% --- Plot 2: Control Input ---
fprintf("Visualizing the control input over time...\n");
fig2 = figure('Position', [900, 100, 800, 600]);
plot(1:length(best_U_traj), best_U_traj, 'b-', 'LineWidth', 1.5);

% Setup plot with LaTeX interpreter
title('Control Input ''$u$'' During Simulation', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Simulation Step', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Input Value ($u$)', 'Interpreter', 'latex', 'FontSize', 12);
legend({'Control Input $u$'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid on;
box on;

% Save the figure as a PDF
fprintf("Saving control input plot as PDF...\n");
exportgraphics(fig2, 'control_input_plot.pdf', 'ContentType', 'vector');
fprintf("Done.\n");


%% --- Helper Functions ---

function create_dummy_policy_file(filename, num_states, num_actions)
    if ~isfile(filename)
        fprintf("Creating dummy policy file: %s\n", filename);
        dummy_policy = round(num_actions / 2) * ones(num_states, 1);
        writematrix(dummy_policy, filename, 'Delimiter', ';');
    end
end

function safe = is_point_safe(x, y, L_straight, R_outer, R_inner)
    xc_left = -L_straight / 2.0;
    xc_right = L_straight / 2.0;
    y_c = 0.0;

    is_in_straights = (x >= xc_left && x <= xc_right) && (abs(y - y_c) <= R_outer);
    is_in_left_curve = x < xc_left && ((x - xc_left)^2 + (y - y_c)^2 <= R_outer^2);
    is_in_right_curve = x > xc_right && ((x - xc_right)^2 + (y - y_c)^2 <= R_outer^2);
    inside_outer = is_in_straights || is_in_left_curve || is_in_right_curve;
    
    if R_inner > 0
        is_in_inner_straights = (x >= xc_left && x <= xc_right) && (abs(y - y_c) < R_inner);
        is_in_inner_left_curve = x < xc_left && ((x - xc_left)^2 + (y - y_c)^2 < R_inner^2);
        is_in_inner_right_curve = x > xc_right && ((x - xc_right)^2 + (y - y_c)^2 < R_inner^2);
        inside_inner = is_in_inner_straights || is_in_inner_left_curve || is_in_inner_right_curve;
    else
        inside_inner = false;
    end
    
    safe = inside_outer && ~inside_inner;
end

function [x1_next, x2_next, x3_next] = di_dynamics(x1, x2, x3, u, d)
    dt = 0.2;
    V = 0.25;
    x1_next = x1 + V * cos(x3) * dt;
    x2_next = x2 + V * sin(x3) * dt;
    x3_next = x3 + (u + d) * dt;
end

function [X_traj, Y_traj, U_traj] = run_simulation(x_start, y_start, theta_start, tree, policy, u_actions, max_time, L, R_out, R_in)
    X_traj = [];
    Y_traj = [];
    U_traj = [];
    xn = x_start;
    yn = y_start;
    thetan = theta_start;

    for i = 1:max_time
        if ~is_point_safe(xn, yn, L, R_out, R_in)
            break;
        end

        X_traj = [X_traj; xn];
        Y_traj = [Y_traj; yn];
        
        % Find nearest state in the grid
        idx = knnsearch(tree, [xn, yn, thetan], 'K', 1);
        policy_index = round(policy(idx));
        
        input_u = u_actions(policy_index);
        U_traj = [U_traj; input_u];
        
        disturbance = 0.0;
        [xn, yn, thetan] = di_dynamics(xn, yn, thetan, input_u, disturbance);
    end
end

function [x_points, y_points] = generate_racetrack_points(R, L)
    if R <= 0, x_points = []; y_points = []; return; end
    n_points_curve = 100;
    xc_l = -L / 2.0;
    xc_r = L / 2.0;
    y_c = 0.0;

    t_right = linspace(pi/2, -pi/2, n_points_curve);
    x_r_curve = xc_r + R * cos(t_right);
    y_r_curve = y_c + R * sin(t_right);

    t_left = linspace(3*pi/2, pi/2, n_points_curve);
    x_l_curve = xc_l + R * cos(t_left);
    y_l_curve = y_c + R * sin(t_left);
    
    % Connect the curves with straight lines for a closed shape
    x_points = [x_r_curve, xc_l, x_l_curve, xc_r];
    y_points = [y_r_curve, -R, y_l_curve, R];
end