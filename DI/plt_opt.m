clear; clc;

% --- Parameters (MUST match Julia Script) ---
num_points_x1 = 101; % Position
num_points_x2 = 101; % Velocity

% Define axes vectors
x1_vec = linspace(-1.0, 5.0, num_points_x1); % Position (Y-axis)
x2_vec = linspace(-5.0, 5.0, num_points_x2); % Velocity (X-axis)

% Create Meshgrid for Contour Plotting
% X1 = Position (Cols), X2 = Velocity (Rows)
[X1, X2] = meshgrid(x1_vec, x2_vec);

% ==========================================
% 1. LOAD DATA
% ==========================================
disp('Loading data...');
try
    % Value Iteration Data
    VI_V_data = readmatrix('value_function.csv');
    VI_P_data = readmatrix('optimal_policy.csv');
    
    % Policy Iteration Data
    PI_V_data = readmatrix('value_function_pi.csv');
    PI_P_data = readmatrix('optimal_policy_pi.csv');
catch
    error('Error loading files. Ensure all 4 CSV files (VI and PI) are in the current directory.');
end

% ==========================================
% 2. RESHAPE & CLEAN
% ==========================================
% Reshape (Rows=Velocity, Cols=Position)
VI_V_grid = reshape(VI_V_data, num_points_x2, num_points_x1);
VI_P_grid = reshape(VI_P_data, num_points_x2, num_points_x1);
PI_V_grid = reshape(PI_V_data, num_points_x2, num_points_x1);
PI_P_grid = reshape(PI_P_data, num_points_x2, num_points_x1);

% Clean Unsafe Values (Thresholding)
unsafe_thresh = 100; 
VI_V_grid(VI_V_grid > unsafe_thresh) = NaN;
PI_V_grid(PI_V_grid > unsafe_thresh) = NaN;

% Transparency Maps for Policies
alpha_VI = ~isnan(VI_P_grid');
alpha_PI = ~isnan(PI_P_grid');

% ==========================================
% 3. PLOTTING (2x2 Grid)
% ==========================================
figure('Name', 'Comparison: Value Iteration vs Policy Iteration', ...
       'Color', 'w', 'Position', [50, 50, 1200, 900]);

% ---------------------------------------------------------
% ROW 1: VALUE ITERATION
% ---------------------------------------------------------

% Plot 1: VI Value Function
subplot(2, 2, 1);
contourf(X2, X1, VI_V_grid, 20, 'LineStyle', 'none');
colormap(gca, 'parula'); colorbar;
title('VI: Value Function');
ylabel('Position (x_1)'); xlabel('Velocity (x_2)');
axis square; grid on;

% Plot 2: VI Optimal Policy
subplot(2, 2, 2);
imagesc(x2_vec, x1_vec, VI_P_grid'); 
set(gca, 'YDir', 'normal');
colormap(gca, 'jet'); colorbar; caxis([-2, 2]);
title('VI: Optimal Policy (u)');
ylabel('Position (x_1)'); xlabel('Velocity (x_2)');
axis square;
% Apply transparency
set(gca, 'Color', [0.8 0.8 0.8]);
set(findobj(gca, 'Type', 'Image'), 'AlphaData', alpha_VI);

% ---------------------------------------------------------
% ROW 2: POLICY ITERATION
% ---------------------------------------------------------

% Plot 3: PI Value Function
subplot(2, 2, 3);
contourf(X2, X1, PI_V_grid, 20, 'LineStyle', 'none');
colormap(gca, 'parula'); colorbar;
title('PI: Value Function');
ylabel('Position (x_1)'); xlabel('Velocity (x_2)');
axis square; grid on;

% Plot 4: PI Optimal Policy
subplot(2, 2, 4);
imagesc(x2_vec, x1_vec, PI_P_grid'); 
set(gca, 'YDir', 'normal');
colormap(gca, 'jet'); colorbar; caxis([-2, 2]);
title('PI: Optimal Policy (u)');
ylabel('Position (x_1)'); xlabel('Velocity (x_2)');
axis square;
% Apply transparency
set(gca, 'Color', [0.8 0.8 0.8]);
set(findobj(gca, 'Type', 'Image'), 'AlphaData', alpha_PI);