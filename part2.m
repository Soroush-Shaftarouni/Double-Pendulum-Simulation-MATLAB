clc; clear; close all;

p.m1 = 1; p.m2 = 1; p.L1 = 1; p.L2 = 1; p.g = 9.81;

z0 = [pi/2; 0; pi/2; 0]; 
tbaze = linspace(0, 20, 2000);

[t, z] = ode45(@(t, z) pendulum_equations(t, z, p), tbaze, z0);

teta1 = z(:, 1); omega1 = z(:, 2);
teta2 = z(:, 3); omega2 = z(:, 4);

x2 = p.L1 * sin(teta1) + p.L2 * sin(teta2);
y2 = -p.L1 * cos(teta1) - p.L2 * cos(teta2);

figure('Name', 'مسیر حرکت');
plot(x2, y2, 'k', 'LineWidth', 1);
xlabel('x_2'); ylabel('y_2'); title('مسیر حرکت جرم دوم'); grid on;

figure('Name', 'زوایا');
plot(t, teta1, 'r', t, teta2, 'b', 'LineWidth', 1.5);
xlabel('زمان'); ylabel('زاویه'); title('زوایا بر حسب زمان'); legend('teta_1', 'teta_2'); grid on;

figure('Name', 'سرعت های زاویه ای');
plot(t, omega1, 'r', t, omega2, 'b', 'LineWidth', 1.5);
xlabel('زمان'); ylabel('سرعت زاویه ای'); title('سرعت زاویه ای بر حسب زمان'); legend('\omega_1', '\omega_2'); grid on;
