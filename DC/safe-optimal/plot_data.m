clear; clc; close all;

% =========================================================================
% PART 1: SETUP & CONFIGURATION
% =========================================================================
result_dir = pwd; 
fprintf('Processing data in: %s\n', result_dir);

% Grid & Geometry Config
nx = 55; ny = 55; nth = 32;
L_straight = 2.0;
R_outer    = 1.80;
R_inner    = 0.90;
Obs_x      = 0.0; Obs_y = 1.1; Obs_r = 0.4;
y_c        = 0.0;
start_x = -2.5; start_y = -0.15;

% Finish Line Geometry
finish_line_x = -1.3;
finish_dx = abs(finish_line_x) - L_straight/2.0;
finish_y_inner = -sqrt(R_inner^2 - finish_dx^2);
finish_y_outer = -sqrt(R_outer^2 - finish_dx^2);

% --- File Finding ---
files = dir(fullfile(result_dir, '*.csv'));
file_names = {files.name};
find_file = @(kwd) file_names(contains(lower(file_names), lower(kwd)));

f_risky = find_file('risky_trajectories');
f_safe  = find_file('safe_trajectories');
f_stats = find_file('stats');

if isempty(f_risky) || isempty(f_safe)
    error('Trajectory CSV files not found in directory.');
end

risky_file = fullfile(result_dir, f_risky{1});
safe_file  = fullfile(result_dir, f_safe{1});

% =========================================================================
% PART 2: TRAJECTORY PLOT (FIGURE 1)
% =========================================================================
fprintf('Generating Trajectory Plot...\n');

% Load Trajectory Data
opts = detectImportOptions(risky_file);
opts.VariableTypes = {'double', 'double', 'double', 'double', 'string'}; 
risky_data = readtable(risky_file, opts);
risky_data.Policy = repmat("Risky", height(risky_data), 1);

safe_data = readtable(safe_file, opts);
safe_data.Policy = repmat("Safe", height(safe_data), 1);

all_trajectories = [risky_data; safe_data];

% --- CALCULATE STATS (RESTORED LOGIC) ---
if ~isempty(f_stats)
    stats_file = fullfile(result_dir, f_stats{1});
    opts_stats = detectImportOptions(stats_file);
    opts_stats.VariableTypes = {'string', 'double', 'double', 'double'};
    sim_stats = readtable(stats_file, opts_stats);
else
    % Manual Calculation (Guarantees text box appears)
    fprintf('Stats file not found. Calculating manually...\n');
    policies = unique(all_trajectories.Policy);
    stats_cells = cell(length(policies), 4);
    for i = 1:length(policies)
        p = policies(i);
        sub_data = all_trajectories(all_trajectories.Policy == p, :);
        ids = unique(sub_data.RunID);
        n_succ = 0; n_stuck = 0; n_crash = 0;
        for j = 1:length(ids)
            run = sub_data(sub_data.RunID == ids(j), :);
            if isempty(run.Status), stat = "unknown"; else, stat = string(run.Status(end)); end
            if strcmpi(stat, 'success'), n_succ = n_succ+1; 
            elseif strcmpi(stat, 'stuck'), n_stuck = n_stuck+1; 
            elseif strcmpi(stat, 'crash'), n_crash = n_crash+1; end
        end
        stats_cells(i,:) = {p, n_succ, n_stuck, n_crash};
    end
    sim_stats = cell2table(stats_cells, 'VariableNames', {'Policy', 'Success', 'No_Feasible', 'Crash'});
end

% Prepare Grid for Geometry Drawing
x_grid = linspace(-3.0, 3.0, nx);
y_grid = linspace(-3.0, 3.0, ny);
[X, Y] = meshgrid(x_grid, y_grid);
dx_rect = max(abs(X) - L_straight/2.0, 0.0);
dy_rect = abs(Y - y_c);
d_center = sqrt(dx_rect.^2 + dy_rect.^2);

% --- FIGURE 1 SETUP ---
fig1 = figure('Name', 'Trajectory Simulation Results', 'Color', 'w', 'Position', [50, 50, 1000, 700]);
hold on; axis equal;
xlim([-3.45, 3.2]); ylim([-2., 2.5]);
xlabel('X (m)'); ylabel('Y (m)');
grid on; box on;

% Styling
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 20, ...
    'FontWeight', 'bold', 'LineWidth', 2.0, 'TickLength', [0.02 0.02]);

% Draw Environment
contour(X, Y, d_center, [R_outer R_outer], 'k', 'LineWidth', 3);
contour(X, Y, d_center, [R_inner R_inner], 'k', 'LineWidth', 3);
th_circ = linspace(0, 2*pi, 100);
fill(Obs_x + Obs_r*cos(th_circ), Obs_y + Obs_r*sin(th_circ), [0.3, 0.3, 0.3], 'EdgeColor', 'none');
plot([finish_line_x, finish_line_x], [finish_y_outer, finish_y_inner], 'g-', 'LineWidth', 5);
plot(start_x, start_y, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');

% Draw Trajectories
policies = unique(all_trajectories.Policy);
colors_map = containers.Map({'Risky', 'Safe'}, {[1, 0, 0], [0, 0, 1]}); 

for p = 1:length(policies)
    pol_name = string(policies(p));
    pol_data = all_trajectories(all_trajectories.Policy == pol_name, :);
    if isKey(colors_map, pol_name), base_color = colors_map(pol_name); else, base_color = [0,0,0]; end
    plot_color = [base_color, 0.3]; 
    
    ids = unique(pol_data.RunID);
    for i = 1:length(ids)
        run_data = pol_data(pol_data.RunID == ids(i), :);
        x_path = run_data.X; y_path = run_data.Y;
        if isempty(run_data.Status), final_status = "unknown"; else, final_status = string(run_data.Status(end)); end
        
        plot(x_path, y_path, 'Color', plot_color, 'LineWidth', 2);
        if strcmpi(final_status, "crash")
            plot(x_path(end), y_path(end), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
        elseif strcmpi(final_status, "stuck")
            plot(x_path(end), y_path(end), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        end
    end
end

% Legend & Annotation
h_risky = plot(NaN, NaN, 'Color', [1, 0, 0], 'LineWidth', 2);
h_safe  = plot(NaN, NaN, 'Color', [0, 0, 1], 'LineWidth', 2);
h_crash = plot(NaN, NaN, 'kx', 'MarkerSize', 15, 'LineWidth', 2);
h_stuck = plot(NaN, NaN, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
h_finish = plot(NaN, NaN, 'g-', 'LineWidth', 5);
h_start = plot(NaN, NaN, 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');

legend([h_risky, h_safe, h_crash, h_stuck, h_finish, h_start], ...
    {'\alpha \geq 0.70 (Red)', '\alpha \geq 0.95 (Blue)', 'Crash', 'Stuck', 'Finish Line', ...
    'Initial State: [-2.5, -0.15, \pi/4]'}, ...
    'Location', 'southwest', 'Interpreter', 'tex','FontSize', 14);

if ~isempty(sim_stats)
    stats_str = { '\bf{Simulation Results (100 Monte Carlo)}' };
    for i = 1:height(sim_stats)
        p_raw = string(sim_stats.Policy(i));
        succ = sim_stats.Success(i); stuck = sim_stats.No_Feasible(i); crash = sim_stats.Crash(i);
        if contains(p_raw, "Risky", 'IgnoreCase', true)
            c = '\color{red}'; p_name = "\alpha \geq 0.70";
        else
            c = '\color{blue}'; p_name = "\alpha \geq 0.95";
        end
        stats_str{end+1} = sprintf('%s%s: Success: %d | Stuck: %d | Crash: %d', c, p_name, succ, stuck, crash); %#ok<AGROW>
    end
    annotation('textbox', [0.52, 0.76, 0.35, 0.12], 'String', stats_str, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'w', 'FaceAlpha', 0.9, ...
        'Interpreter', 'tex', 'FontSize', 15, 'EdgeColor', 'k');
end
hold off;

% Save Figure 1
print(fig1, 'Simulation_Trajectories.pdf', '-dpdf', '-r0', '-painters');
% =========================================================================
% PART 3: 3x3 HEATMAP GRID (FIGURE 2)
% =========================================================================
fprintf('Generating 3x3 Heatmap Grid...\n');

% Config
trunc_val = 10^7;

target_thetas = [1, 9, 17]; 
theta_labels  = {'0', '\pi/2', '\pi'};

% Load Heatmap Data
gain_raw      = readmatrix(fullfile(result_dir, 'safety_values_g.csv'));
val_risky_raw = readmatrix(fullfile(result_dir, 'value_function_risky.csv'));
val_safe_raw  = readmatrix(fullfile(result_dir, 'value_function_safe.csv'));

G_3D       = reshape(gain_raw, nth, nx, ny);
V_risky_3D = reshape(val_risky_raw, nth, nx, ny);
V_safe_3D  = reshape(val_safe_raw, nth, nx, ny);

% Safety Mask
[Xm, Ym] = meshgrid(x_grid, y_grid); 
dx_m = max(abs(Xm) - L_straight/2.0, 0.0);
dy_m = abs(Ym - y_c);
dm = sqrt(dx_m.^2 + dy_m.^2);
in_track_m = (dm <= R_outer) & (dm >= R_inner);
in_obs_m   = (Xm - Obs_x).^2 + (Ym - Obs_y).^2 <= Obs_r^2;
safe_mask  = in_track_m & ~in_obs_m;

% --- FIGURE 2 SETUP ---
fig2 = figure('Name', 'Heatmap Comparison 3x3', 'Color', 'w', 'Position', [100, 100, 1200, 1000]);

for row = 1:3
    t_idx = target_thetas(row);
    t_lbl = theta_labels{row};
    
    % Extract Slices
    slice_G  = squeeze(G_3D(t_idx, :, :))';
    slice_Vr = squeeze(V_risky_3D(t_idx, :, :))';
    slice_Vs = squeeze(V_safe_3D(t_idx, :, :))';
    
    % Mask & Truncate
    slice_G(~safe_mask)  = NaN;
    slice_Vr(~safe_mask) = NaN;
    slice_Vs(~safe_mask) = NaN;
    
    slice_Vr(slice_Vr > trunc_val) = trunc_val;
    slice_Vs(slice_Vs > trunc_val) = trunc_val;
    
    % --- PLOT COLUMN 1: GAIN ---
    sub_idx = (row-1)*3 + 1;
    subplot(3, 3, sub_idx);
    im = imagesc(x_grid, y_grid, slice_G);
    set(im, 'AlphaData', ~isnan(slice_G));
    axis xy equal tight;
    colormap(gca, 'parula'); colorbar;
    
    title(['Safety Gain ($g$) $\theta \approx ' t_lbl '$'], 'Interpreter', 'latex', 'FontSize', 16);
    if row == 3; xlabel('$X$ (m)', 'Interpreter', 'latex', 'FontSize', 14); end
    ylabel('$Y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
    
    % --- PLOT COLUMN 2: RISKY VALUE ---
    sub_idx = (row-1)*3 + 2;
    subplot(3, 3, sub_idx);
    im = imagesc(x_grid, y_grid, slice_Vr);
    set(im, 'AlphaData', ~isnan(slice_Vr));
    axis xy equal tight;
    colormap(gca, 'turbo'); colorbar;
    
    title(['Risky Value $\theta \approx ' t_lbl '$'], 'Interpreter', 'latex', 'FontSize', 16);
    if row == 3; xlabel('$X$ (m)', 'Interpreter', 'latex', 'FontSize', 14); end
    set(gca, 'YTickLabel', [], 'FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
    
    % --- PLOT COLUMN 3: SAFE VALUE ---
    sub_idx = (row-1)*3 + 3;
    subplot(3, 3, sub_idx);
    im = imagesc(x_grid, y_grid, slice_Vs);
    set(im, 'AlphaData', ~isnan(slice_Vs));
    axis xy equal tight;
    colormap(gca, 'turbo'); colorbar;
    
    title(['Safe Value $\theta \approx ' t_lbl '$'], 'Interpreter', 'latex', 'FontSize', 16);
    if row == 3; xlabel('$X$ (m)', 'Interpreter', 'latex', 'FontSize', 14); end
    set(gca, 'YTickLabel', [], 'FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
end

sgtitle('\textbf{Safety and Value Function at Different Orientations}', 'Interpreter', 'latex', 'FontSize', 22);

% Save Figure 2
print(fig2, 'Heatmaps_3x3.pdf', '-dpdf', '-r0', '-painters');

fprintf('Done. Saved Simulation_Trajectories.pdf and Heatmaps_3x3.pdf\n');