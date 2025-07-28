clc; clear; close all;

% Ellipse parameters
x_c = 0; v_c = 0;

a_outer = 0.4; b_outer = 0.9;
a_inner = 0.2;  b_inner = 0.7;

theta = linspace(0, 2*pi, 500);

% Outer ellipse points
x_outer = x_c + a_outer * cos(theta);
v_outer = v_c + b_outer * sin(theta);

% Inner ellipse points
x_inner = x_c + a_inner * cos(theta);
v_inner = v_c + b_inner * sin(theta);

figure; hold on; grid on; axis equal

% Plot outer ellipse
plot(x_outer, v_outer, 'b-', 'LineWidth', 2, 'DisplayName', 'Outer Ellipse');

% Plot inner ellipse
plot(x_inner, v_inner, 'r-', 'LineWidth', 2, 'DisplayName', 'Inner Ellipse');

xlabel('x');
ylabel('v');
title('Ring between Two Ellipses (Safe Set)');
legend('Location', 'best');

% Fix axis limits
xlim([-0.5, 0.5]);
ylim([-1, 1]);
