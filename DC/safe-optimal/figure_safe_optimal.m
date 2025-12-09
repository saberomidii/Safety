clear; clc; close all;

% --- CONFIGURATION ---
result_dir = 'safe-optimal';
fixed_theta_idx = 16; % Choose a theta index to slice (e.g., facing Right/Up)

% UPDATED GRID DIMENSIONS (Must match Julia: num_points_x/y/th)
nx = 76; 
ny = 76;
nth = 31; 

% GEOMETRY PARAMETERS (Must match Julia config)
L_straight = 2.0;
R_outer = 1.80;
R_inner = 0.90;
Obs_x = 0.0;
Obs_y = 1.1;
Obs_r = 0.4; % Based on your code snippet
y_c = 0.0;

% Targets and Start (Must match Julia config)
start_x = -2.5;
start_y = -0.15;
target_x = 1.5;
target_y = 1.35;

% Grid Definitions (Must match Julia config)
x_grid = linspace(-3.0, 3.0, nx);
y_grid = linspace(-3.0, 3.0, ny);
[X, Y] = meshgrid(x_grid, y_grid);

% --- 1. LOAD DATA ---

% Load Safety Gain (g)
g_flat = load(fullfile(result_dir, 'safety_values_g.csv'));

% Reshape: Julia saves as [theta, x, y] order typically or linear index
% Adjust reshaping based on Julia's linear indexing (Column-major: Y then X then Th)
% We use [] for the 3rd dimension to auto-calculate size and prevent errors
g_3d = reshape(g_flat, ny, nx, []);

% Load Value Functions
V_risky_flat = load(fullfile(result_dir, 'value_function_risky.csv'));
V_safe_flat  = load(fullfile(result_dir, 'value_function_safe.csv'));

V_risky_3d = reshape(V_risky_flat, ny, nx, []);
V_safe_3d  = reshape(V_safe_flat,  ny, nx, []);

% Load Trajectories (CSV)
% Format: RunID, Step, X, Y, Status
traj_risky_data = readtable(fullfile(result_dir, 'risky_trajectories.csv'));
traj_safe_data  = readtable(fullfile(result_dir, 'safe_trajectories.csv'));


% --- 3. PLOT SIMULATION TRAJECTORIES WITH CONSTRAINT SET ---

figure('Name', '100 Montro Carlo Sampling Simulation Results', 'Color', 'w');
hold on; axis equal;
xlim([-3, 3]); ylim([-2.25, 2.25]);
xlabel('X (m)'); ylabel('Y (m)');

% --- DRAW CONSTRAINT SET (TRACK) ---
% Re-calculate constraint set locally for plotting boundary
dx_rect = max(abs(X) - L_straight/2.0, 0.0);
dy_rect = abs(Y - y_c);
d_center = sqrt(dx_rect.^2 + dy_rect.^2);
in_track = (d_center <= R_outer) & (d_center >= R_inner);
in_obs = (X - Obs_x).^2 + (Y - Obs_y).^2 <= Obs_r^2;
constraint_mask = in_track & ~in_obs;

% Draw Outer Wall (Red)
contour(X, Y, d_center, [R_outer R_outer], 'black', 'LineWidth', 4);
% Draw Inner Wall (Red)
contour(X, Y, d_center, [R_inner R_inner], 'black', 'LineWidth', 4);
% Draw Obstacle (Black)
th_obs = linspace(0, 2*pi, 100);
fill(Obs_x + Obs_r*cos(th_obs), Obs_y + Obs_r*sin(th_obs), 'black');

% --- PLOT TRAJECTORIES ---
% Plot Risky Trajectories (Orange)
ids = unique(traj_risky_data.RunID);
for i = 1:min(length(ids), 30) 
    idx = traj_risky_data.RunID == ids(i);
    path_x = traj_risky_data.X(idx);
    path_y = traj_risky_data.Y(idx);
    status = string(traj_risky_data.Status(find(idx, 1, 'last')));
    plot(path_x, path_y, 'Color', [1, 0.5, 0, 0.4], 'LineWidth', 3);
    if status == "crash"
        plot(path_x(end), path_y(end), 'rx', 'MarkerSize', 12, 'LineWidth', 3,'Color','r');
    end
end

% Plot Safe Trajectories (Magenta)
ids = unique(traj_safe_data.RunID);
for i = 1:min(length(ids), 30)
    idx = traj_safe_data.RunID == ids(i);
    path_x = traj_safe_data.X(idx);
    path_y = traj_safe_data.Y(idx);
    status = string(traj_safe_data.Status(find(idx, 1, 'last')));
    plot(path_x, path_y, 'Color', [1, 0, 1, 0.4], 'LineWidth', 3);
    if status == "crash"
        plot(path_x(end), path_y(end), 'rx', 'MarkerSize', 12, 'LineWidth', 3,'Color','r');
    end
end

% Plot Start and Target Points
plot(start_x, start_y, 'ob', 'MarkerSize', 25, 'MarkerFaceColor', 'blue', 'LineWidth', 2.5); % Start (Cyan)
plot(target_x, target_y, 'og', 'MarkerSize', 30, 'MarkerFaceColor', 'green', 'LineWidth', 3.5); % Target (Red)

% Fake Legend for Policies
h1 = plot(NaN,NaN,'Color',[1, 0.5, 0], 'LineWidth', 2);
h2 = plot(NaN,NaN,'Color',[1, 0, 1], 'LineWidth', 2);
h3 = plot(NaN,NaN,'black', 'LineWidth', 2); % Wall
h4 = fill(NaN,NaN, 'k'); % Obstacle
h5 = plot(NaN,NaN, 'ob', 'MarkerFaceColor', 'blue'); % Start
h6 = plot(NaN,NaN, 'og', 'MarkerFaceColor', 'green'); % Target
h7 = plot(NaN,NaN, 'xr', 'MarkerFaceColor', 'red'); % Target

legend([h1, h2, h3, h4, h5, h6,h7], {'Risky Policy', 'Safe Policy', 'Constraint Boundary', 'Obstacle', 'Start', 'Target','Crash'}, 'Location', 'best');

hold off;

