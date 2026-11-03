clc; clear; close all;

p.m1 = 1; p.m2 = 1; p.L1 = 1; p.L2 = 1; p.g = 9.81;
z0 = [pi/2; 0; pi/2; 0];
target_angle = pi/4;

t_nr = 0.5; 
tol = 1e-4;
iter_nr = 0;
for i = 1:50
    [f_val, f_der] = eval_teta(t_nr, p, z0, target_angle);
    if abs(f_val) < tol
        break;
    end
    t_nr = t_nr - f_val / f_der;
    iter_nr = iter_nr + 1;
end
fprintf('--- نیوتن-رافسون ---\nزمان: %f | تکرار: %d\n\n', t_nr, iter_nr);

t_fp = 0.5; 
k = 0.1; 
iter_fp = 0;
for i = 1:200
    [f_val, ~] = eval_teta(t_fp, p, z0, target_angle);
    if abs(f_val) < tol
        break;
    end
    t_fp = t_fp - k * f_val; 
    iter_fp = iter_fp + 1;
end
fprintf('--- نقطه ثابت ---\nزمان: %f | تکرار: %d\n\n', t_fp, iter_fp);

a = 0.1; b = 1.0; 
iter_bi = 0;
for i = 1:50
    c = (a + b) / 2;
    [f_c, ~] = eval_teta(c, p, z0, target_angle);
    if abs(f_c) < tol || (b - a)/2 < tol
        break;
    end
    [f_a, ~] = eval_teta(a, p, z0, target_angle);
    if sign(f_c) == sign(f_a)
        a = c;
    else
        b = c;
    end
    iter_bi = iter_bi + 1;
end
fprintf('--- نصف کردن (Bisection) ---\nزمان: %f | تکرار: %d\n', c, iter_bi);


function [f_val, f_der] = eval_teta(t_end, p, z0, target)
    if t_end <= 0
        f_val = z0(1) - target;
        f_der = z0(2);
        return;
    end
    [~, z] = ode45(@(t, z) pendulum_equations(t, z, p), [0, t_end], z0);
    f_val = z(end, 1) - target; 
    f_der = z(end, 2);          
end
